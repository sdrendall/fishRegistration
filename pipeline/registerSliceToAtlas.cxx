#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"

#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkBSplineTransform.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"

#include "itkDemonsRegistrationFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"

#include "itkPermuteAxesImageFilter.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkCommand.h"

const char * atlasReferencePath = "/home/sam/Pictures/allenReferenceAtlas_mouseCoronal/atlasVolume/atlasVolume.mhd";
const char * atlasLabelsPath = "/home/sam/Pictures/allenReferenceAtlas_mouseCoronal/annotation.mhd";
const unsigned int BSplineOrder = 3;

typedef itk::Image<unsigned char, 2> ImageType;
typedef itk::CenteredRigid2DTransform<double> RigidTransformType;
typedef itk::BSplineTransform<double, 2, BSplineOrder> BSplineTransformType; // <CoordinateRepType, Dims, BSplineOrder>

// Types for the demons registration
typedef itk::Vector<float, 2> DisplacementPixelType;
typedef itk::Image<DisplacementPixelType, 2> DisplacementFieldType;
typedef itk::WarpImageFilter<ImageType, ImageType, DisplacementFieldType> WarperType;

// Function Declarations
ImageType::Pointer getCoronalAtlasSlice(int, const char *);
ImageType::Pointer rotateImage(ImageType::Pointer);
RigidTransformType::Pointer getRigidRegistrationTransform(ImageType::Pointer, ImageType::Pointer);
BSplineTransformType::Pointer getBSPlineRegistrationResults(ImageType::Pointer, ImageType::Pointer);
DisplacementFieldType::Pointer getDemonsDisplacementField(ImageType::Pointer, ImageType::Pointer);
WarperType::Pointer createAndConfigureDemonsWarper(DisplacementFieldType::Pointer, ImageType::Pointer);

double degreesToRadians(double);
void writeImage(ImageType::Pointer, const char *);
void displayRegistrationResults(itk::OptimizerParameters<double>, const unsigned int, const double);
void debugOut(const char *);

// Observer class, to output the progress of the registration
// Both observers are copied directly from examples, without alteration
class CommandIterationUpdate : public itk::Command {
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};

public:
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef   const OptimizerType *                  OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    OptimizerPointer optimizer =
      dynamic_cast< OptimizerPointer >( object );
    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};

class BSplineTransformIterationUpdate : public itk::Command
{
public:
  typedef  BSplineTransformIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );

protected:
  BSplineTransformIterationUpdate() {};

public:
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef   const OptimizerType *                  OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    OptimizerPointer optimizer =
      dynamic_cast< OptimizerPointer >( object );
    if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
      return;
      }
    std::cout << "Iteration : ";
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << std::endl;
    }
};

  class DemonsIterationUpdate : public itk::Command
  {
  public:
    typedef  DemonsIterationUpdate                     Self;
    typedef  itk::Command                              Superclass;
    typedef  itk::SmartPointer<DemonsIterationUpdate>  Pointer;
    itkNewMacro( DemonsIterationUpdate );
  protected:
    DemonsIterationUpdate() {};

    typedef itk::Image< float, 2 >            InternalImageType;
    typedef itk::Vector< float, 2 >           VectorPixelType;
    typedef itk::Image<  VectorPixelType, 2 > DisplacementFieldType;

    typedef itk::DemonsRegistrationFilter<
                                InternalImageType,
                                InternalImageType,
                                DisplacementFieldType>   RegistrationFilterType;

  public:

    void Execute(itk::Object *caller, const itk::EventObject & event)
      {
        Execute( (const itk::Object *)caller, event);
      }

    void Execute(const itk::Object * object, const itk::EventObject & event)
      {
         const RegistrationFilterType * filter =
          dynamic_cast< const RegistrationFilterType * >( object );
        if( !(itk::IterationEvent().CheckEvent( &event )) )
          {
          return;
          }
        std::cout << filter->GetMetric() << std::endl;
      }
  };

int main(int argc, char *argv[]){
    //Check args
    if (argc < 5) {
            std::cerr << "Missing Parameters " << std::endl;
            std::cerr << "Usage: " << argv[0];
            std::cerr << " fixedImage sliceDepth";
            std::cerr << " outputAtlasImage outputAtlasLabels";
            return EXIT_FAILURE;
    }

    // Declare Image, filter types ect
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::ImageFileWriter<ImageType> WriterType;

    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;

    // The annotated images must be resampled with a nearest neighbor interpolator
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NearestNeighborInterpolatorType;


    // Get corresponding atlas reference slice
    ImageType::Pointer atlasSlice = getCoronalAtlasSlice(atoi(argv[2]), atlasReferencePath);
    ImageType::Pointer annotationSlice = getCoronalAtlasSlice(atoi(argv[2]), atlasLabelsPath);

    // Load the input image
    ReaderType::Pointer inputReader = ReaderType::New();
    inputReader->SetFileName(argv[1]);
    inputReader->Update();
    ImageType::Pointer inputImage = inputReader->GetOutput();

    // Set the input spacing and atlas direction
    ImageType::SpacingType spacing;
    spacing[0] = 25;
    spacing[1] = 25;
    inputImage->SetSpacing(spacing);
    atlasSlice->SetDirection(inputImage->GetDirection());
    annotationSlice->SetDirection(inputImage->GetDirection());

    //FOR DEBUGGING - Save the slice locally, until I feel more confident
    writeImage(atlasSlice, "/home/sam/Desktop/extractedSlice.jpg");
    writeImage(annotationSlice, "/home/sam/Desktop/extAnnoSlice.jpg");
    writeImage(inputImage, "/home/sam/Desktop/inputImage.jpg");

    RigidTransformType::Pointer rigidTransform = getRigidRegistrationTransform(inputImage, atlasSlice);
    
    // Resample the atlas using the computed rigid transform
    ResampleFilterType::Pointer rigidResampler = ResampleFilterType::New();
    rigidResampler->SetTransform(rigidTransform);
    rigidResampler->SetInput(atlasSlice);

    rigidResampler->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    rigidResampler->SetOutputOrigin(inputImage->GetOrigin());
    rigidResampler->SetOutputSpacing(inputImage->GetSpacing());
    rigidResampler->SetOutputDirection(inputImage->GetDirection());
    rigidResampler->SetDefaultPixelValue(0);

    // Update, to resample the image so it can be used for the deformable registration
    rigidResampler->Update();
    writeImage(rigidResampler->GetOutput(), "/home/sam/Desktop/afterRigidRegistration.jpg");

    // Compute a deformable registration transform using the resampled atlas slice and the input image
    BSplineTransformType::Pointer deformableTransform = getBSPlineRegistrationResults(rigidResampler->GetOutput(), inputImage); // (movingImage, fixedImage)

    ResampleFilterType::Pointer deformableResampler = ResampleFilterType::New();
    deformableResampler->SetTransform(deformableTransform);
    deformableResampler->SetInput(rigidResampler->GetOutput());

    // This should really be done with a function...
    deformableResampler->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    deformableResampler->SetOutputOrigin(inputImage->GetOrigin());
    deformableResampler->SetOutputSpacing(inputImage->GetSpacing());
    deformableResampler->SetOutputDirection(inputImage->GetDirection());
    deformableResampler->SetDefaultPixelValue(0);

    // Compute the demons registration displacement field and create a warper to manipulate images with
    DisplacementFieldType::Pointer displacementField = getDemonsDisplacementField(atlasSlice, inputImage); // (movingImage, fixedImage)
    WarperType::Pointer demonsWarper = createAndConfigureDemonsWarper(displacementField, inputImage); // (displacementField, targetImage)

    // Add the warper to the pipeline
    demonsWarper->SetInput(deformableResampler->GetOutput());

    // Write the registered reference image to the specified filepath
    WriterType::Pointer outputWriter = WriterType::New();
    outputWriter->SetFileName(argv[3]);
    outputWriter->SetInput(demonsWarper->GetOutput());
    outputWriter->Update();

    // Transform the annotated atlas slice and write the output to the specified filepath
    // Create a nearest neighbor interpolator to use for the interpolation
    NearestNeighborInterpolatorType::Pointer nnInterpolator = NearestNeighborInterpolatorType::New();
    // Use the nnInterpolator for each resampler in the pipeline
    rigidResampler->SetInterpolator(nnInterpolator);
    deformableResampler->SetInterpolator(nnInterpolator);
    demonsWarper->SetInterpolator(nnInterpolator);

    // The rigid resampler is the start of the pipeline
    rigidResampler->SetInput(annotationSlice);
    // Updating the output writer will pull the annotated slice through
    outputWriter->SetFileName(argv[4]);
    outputWriter->Update();    

    return EXIT_SUCCESS;
}

ImageType::Pointer getCoronalAtlasSlice(int sliceDepth, const char * atlasPath) {
    // Declare 3D and 2D image types
    typedef unsigned char AtlasPixelType;
    typedef unsigned char SlicePixelType;
    
    typedef itk::Image<AtlasPixelType, 3> AtlasImageType;
    typedef itk::Image<SlicePixelType, 2> SliceImageType;
    
    // Declare Readers and Writers for 3D images and 2D images respectively
    typedef itk::ImageFileReader<AtlasImageType> AtlasReaderType;
    typedef itk::ExtractImageFilter<AtlasImageType, SliceImageType> SliceAtlasFilterType;

    // Collapsing along the x axis will yield coronal sections
    const int axisToCollapse = 0;
    
    // Create a reader, to read the input image
    AtlasReaderType::Pointer reader = AtlasReaderType::New();
    reader->SetFileName(atlasPath);

    // Update, to load the image so that size information can be ascertained
    reader->Update();

    // Get the whole atlas region
    AtlasImageType::Pointer atlas = reader->GetOutput();
    AtlasImageType::RegionType entireAtlasRegion = atlas->GetLargestPossibleRegion();

    // Determine Slice Size
    AtlasImageType::SizeType inputSize = entireAtlasRegion.GetSize();
    AtlasImageType::SizeType sliceSize = inputSize;
    sliceSize[axisToCollapse] = 0;  // 0 tells the ExtractionFilter to return an Image without that dimension

    // Initialize a slice region
    AtlasImageType::RegionType sliceRegion;
    sliceRegion.SetSize(sliceSize);

    // Get Start Point
    AtlasImageType::IndexType sliceStartIndex = entireAtlasRegion.GetIndex();

    // Reference an extraction filter and connect it to the reader
    SliceAtlasFilterType::Pointer extFilter = SliceAtlasFilterType::New();
    extFilter->SetInput(reader->GetOutput());
    extFilter->SetDirectionCollapseToIdentity();

    // Find the corresponding atlas index to the point specified by the user
    typedef itk::Point <double, AtlasImageType::ImageDimension> PointType;
    PointType point;
    point.Fill(0.0);
    point[axisToCollapse] = sliceDepth;
    bool isInside = atlas->TransformPhysicalPointToIndex(point, sliceStartIndex);

    // Extract the slice, if it is in the image
    if (isInside) {
        sliceRegion.SetIndex(sliceStartIndex);
        extFilter->SetExtractionRegion(sliceRegion);
        extFilter->Update();
        SliceImageType::Pointer atlasSlice = extFilter->GetOutput();
        atlasSlice = rotateImage(atlasSlice);  // Think this should work...
        return atlasSlice;
    } else {
        std::cerr << "Plane outside of the atlas boundries!" << std::endl;
        std::cerr << "Please specify a point within the atlas..." << std::endl;
        throw 1337;
    }
}

ImageType::Pointer rotateImage(ImageType::Pointer inputImage) {
    typedef itk::PermuteAxesImageFilter<ImageType> PermuteAxesFilterType;

    // Create a permutation filter
    PermuteAxesFilterType::Pointer permutationFilter = PermuteAxesFilterType::New();

    // Define the axes permutation order
    itk::FixedArray<unsigned int, 2> order;
    order[0] = 1;
    order[1] = 0;
    
    // Permute the axes
    permutationFilter->SetInput(inputImage);
    permutationFilter->SetOrder(order);
    permutationFilter->Update();

    // Return the permuted image
    ImageType::Pointer outputImage = permutationFilter->GetOutput();
    return outputImage;
}

RigidTransformType::Pointer getRigidRegistrationTransform(ImageType::Pointer inputImage, ImageType::Pointer atlasSlice) {
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> MetricType;
    typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
    typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
    typedef itk::CenteredTransformInitializer<RigidTransformType, ImageType, ImageType> TransformInitializerType;

    // Instantiate the metric, optimizer, interpolator and registration objects
    MetricType::Pointer           metric        = MetricType::New();
    OptimizerType::Pointer        optimizer     = OptimizerType::New();
    InterpolatorType::Pointer     interpolator  = InterpolatorType::New();
    RigidTransformType::Pointer   transform     = RigidTransformType::New();
    RegistrationType::Pointer     registration  = RegistrationType::New();

    // Set up the registration
    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetInterpolator(interpolator);
    registration->SetTransform(transform);

    // Set the inputs
    registration->SetFixedImage(inputImage);
    registration->SetMovingImage(atlasSlice);

    registration->SetFixedImageRegion(inputImage->GetBufferedRegion());

    // Create an initializer and connect it to the input images and the transform
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetFixedImage(inputImage);
    initializer->SetMovingImage(atlasSlice);
    initializer->SetTransform(transform);

    // MomentsOn() sets the initializer to center of mass mode
    initializer->MomentsOn();

    // Compute the initial transform using the center of mass
    initializer->InitializeTransform();
    transform->SetAngle(0.0);

    // Get the computed parameters and initialize the transform with them
    registration->SetInitialTransformParameters(transform->GetParameters());

    // Set up the optimizer
    //  The Optimizer is responsible for sustaining and terminating registration iteration
    // I'm not entirely sure what the purpose of scales is
    OptimizerType::ScalesType optimizerScales(transform->GetNumberOfParameters());
    const double translationScale = 25.0/1000.0;

    optimizerScales[0] = 0.0001;
    for (int i = 1; i <= 4; i++) {
        optimizerScales[i] = translationScale;
    }
    optimizer->SetScales(optimizerScales);

    // These parameters determine when the optimizer will terminate registration
    optimizer->SetMaximumStepLength(0.1);
    optimizer->SetMinimumStepLength(0.001);
    optimizer->SetNumberOfIterations(200);

    // The command observer will report on the registrations progress at each iteration
    //  The class definition is below, and taken directly from the ImageRegistration6.cxx example
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    optimizer->AddObserver(itk::IterationEvent(), observer);

    // Begin Registration by calling Update()
    try {
        registration->Update();
        std::cout << "Optimizer stop condition: "
                  << registration->GetOptimizer()->GetStopConditionDescription()
                  << std::endl;
    } catch ( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        throw 1337;
    }

    // Display the final registration transformation parameters
    OptimizerType::ParametersType registrationParameters = registration->GetLastTransformParameters();
    displayRegistrationResults(registrationParameters, optimizer->GetCurrentIteration(), optimizer->GetValue());

    // return the transform computed by the registration
    RigidTransformType::Pointer registrationTransform = RigidTransformType::New();
    registrationTransform->SetParameters(registrationParameters);
    registrationTransform->SetFixedParameters(transform->GetFixedParameters());
    return registrationTransform;
}

BSplineTransformType::Pointer getBSPlineRegistrationResults(ImageType::Pointer movingImage, ImageType::Pointer fixedImage){
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType> MetricType;
    typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
    typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
    typedef BSplineTransformType::ParametersType ParametersType;

    // Instantiate the metric, optimizer, interpolator, transform and registration objects
    MetricType::Pointer metric = MetricType::New();
    OptimizerType::Pointer optimizer = OptimizerType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    RegistrationType::Pointer registration = RegistrationType::New();
    BSplineTransformType::Pointer transform = BSplineTransformType::New();

    // Connect everything to the registration object
    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetInterpolator(interpolator);
    registration->SetTransform(transform);

    // Connect the inputs
    registration->SetFixedImage(fixedImage);
    registration->SetMovingImage(movingImage);
    registration->SetFixedImageRegion(fixedImage->GetBufferedRegion());

    // Calculate image physical dimensions and meshSize
    BSplineTransformType::PhysicalDimensionsType fixedPhysicalDimensions;
    BSplineTransformType::MeshSizeType meshSize;
    BSplineTransformType::OriginType fixedOrigin = fixedImage->GetOrigin();
    ImageType::SizeType fixedImageSize = fixedImage->GetLargestPossibleRegion().GetSize();
    unsigned int numberOfGridNodesInOneDimension = 8;

    for (int i = 0; i < 2; i++) {
        fixedPhysicalDimensions[i] = (fixedImageSize[i] - 1) * fixedImage->GetSpacing()[i];
    }
    meshSize.Fill(numberOfGridNodesInOneDimension - BSplineOrder);

    // Define the transform domain
    transform->SetTransformDomainOrigin(fixedOrigin);
    transform->SetTransformDomainPhysicalDimensions(fixedPhysicalDimensions);
    transform->SetTransformDomainMeshSize(meshSize);
    transform->SetTransformDomainDirection(fixedImage->GetDirection());

    // Set the initial transform parameters
    const unsigned int numberOfParameters = transform->GetNumberOfParameters();
    ParametersType parameters(numberOfParameters);
    parameters.Fill(0.0);
    transform->SetParameters(parameters);
    registration->SetInitialTransformParameters(transform->GetParameters());

    // Specify the optimizer parameters
    optimizer->MaximizeOff();
    optimizer->SetMaximumStepLength(25.0);
    optimizer->SetMinimumStepLength(0.001);
    optimizer->SetRelaxationFactor(0.7);
    optimizer->SetNumberOfIterations(200);

    BSplineTransformIterationUpdate::Pointer observer = BSplineTransformIterationUpdate::New();
    optimizer->AddObserver(itk::IterationEvent(), observer);

    // Prepare time and memory probes
    itk::TimeProbesCollectorBase chronometer;
    itk::MemoryProbesCollectorBase memorymeter;
    
    // Start Registration
    std::cout << std::endl << "Starting Registration" << std::endl;
    try {
        memorymeter.Start("Registration");
        chronometer.Start("Registration");

        registration->Update();

        chronometer.Stop("Registration");
        memorymeter.Stop("Registration");

        std::cout << "Optimizer stop condition = "
                  << registration->GetOptimizer()->GetStopConditionDescription()
                  << std::endl;
    } catch ( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        throw 1337;
    }

    // Return the computed transform
    OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();

    std::cout << "Final BSpline Transform Parameters" << std::endl;
    std::cout << finalParameters << std::endl;

    transform->SetParameters(finalParameters);
    return transform;
}

DisplacementFieldType::Pointer getDemonsDisplacementField(ImageType::Pointer movingImage, ImageType::Pointer fixedImage) {
    typedef itk::Image<float, 2> InternalImageType;
    typedef itk::CastImageFilter<ImageType, InternalImageType> ImageCasterType;
    typedef itk::HistogramMatchingImageFilter<InternalImageType, InternalImageType> MatchingFilterType;
    typedef itk::DemonsRegistrationFilter<InternalImageType, InternalImageType, DisplacementFieldType> DemonsFilterType;

    // Instantiate casters to cast from input image types to the internal image type
    ImageCasterType::Pointer movingCaster = ImageCasterType::New();
    ImageCasterType::Pointer fixedCaster = ImageCasterType::New();

    movingCaster->SetInput(movingImage);
    fixedCaster->SetInput(fixedImage);

    // Create a histogram matcher, and feed it the caster outputs
    MatchingFilterType::Pointer matcher = MatchingFilterType::New();
    matcher->SetInput(movingCaster->GetOutput());
    matcher->SetReferenceImage(fixedCaster->GetOutput());

    // Set matcher parameters
    matcher->SetNumberOfHistogramLevels(1024);
    matcher->SetNumberOfMatchPoints(7);
    // Simple, built in segmentation
    matcher->ThresholdAtMeanIntensityOn();

    // Create the demons registration filter
    DemonsFilterType::Pointer filter = DemonsFilterType::New();

    // Set filter parameters
    filter->SetFixedImage(fixedCaster->GetOutput());
    filter->SetMovingImage(matcher->GetOutput());

    filter->SetNumberOfIterations(50);
    filter->SetStandardDeviations(1);

    // Add an observer
    DemonsIterationUpdate::Pointer observer = DemonsIterationUpdate::New();
    filter->AddObserver(itk::IterationEvent(), observer);

    // Start registration
    filter->Update();

    // Return output
    DisplacementFieldType::Pointer outputField = filter->GetOutput();
    return outputField;
}

WarperType::Pointer createAndConfigureDemonsWarper(DisplacementFieldType::Pointer displacementField, ImageType::Pointer targetImage) {
    // Define warper and interpolator types
    typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;

    // Instantiate a warper and an interpolator
    WarperType::Pointer warper = WarperType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();

    warper->SetInterpolator(interpolator);
    warper->SetOutputSpacing(targetImage->GetSpacing());
    warper->SetOutputOrigin(targetImage->GetOrigin());
    warper->SetOutputDirection(targetImage->GetDirection());
    warper->SetDisplacementField(displacementField);

    return warper;
}

void writeImage(ImageType::Pointer im, const char * path) {
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer fileWriter = WriterType::New();
    fileWriter->SetFileName(path);
    fileWriter->SetInput(im);
    fileWriter->Update();
}

double degreesToRadians(double degrees) {
    const double conversionFactor = vcl_atan(1.0)/45.0;
    double radians = degrees * conversionFactor;
    return radians;
}

void debugOut(const char * msg) {
    std::cout << "[Debug] " << msg << std::endl;
}

void displayRegistrationResults(itk::OptimizerParameters<double> finalParameters, const unsigned int numberOfIterations, const double bestValue) {
    const double finalAngle           = finalParameters[0];
    const double finalRotationCenterX = finalParameters[1];
    const double finalRotationCenterY = finalParameters[2];
    const double finalTranslationX    = finalParameters[3];
    const double finalTranslationY    = finalParameters[4];

    const double finalAngleInDegrees = finalAngle * 180.0 / vnl_math::pi;
    
    std::cout << "Result = " << std::endl;
    std::cout << " Angle (radians) " << finalAngle  << std::endl;
    std::cout << " Angle (degrees) " << finalAngleInDegrees  << std::endl;
    std::cout << " Center X      = " << finalRotationCenterX  << std::endl;
    std::cout << " Center Y      = " << finalRotationCenterY  << std::endl;
    std::cout << " Translation X = " << finalTranslationX  << std::endl;
    std::cout << " Translation Y = " << finalTranslationY  << std::endl;
    std::cout << " Iterations    = " << numberOfIterations << std::endl;
    std::cout << " Metric value  = " << bestValue          << std::endl;
}