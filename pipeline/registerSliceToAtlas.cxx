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

#include "itkBSplineTransform.h"
#include "itkCombinedImageToImageMetricAdaptor.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"

#include "itkPermuteAxesImageFilter.h"

#include "itkExceptionObject.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkCommand.h"

const char * atlasPath = "/home/sam/Pictures/allenReferenceAtlas_mouseCoronal/atlasVolume/atlasVolume.mhd";
const unsigned int BSplineOrder = 3;

typedef itk::Image<unsigned char, 2> ImageType;
typedef itk::CenteredRigid2DTransform<double> RigidTransformType;
typedef itk::BSplineTransform<double, 2, BSplineOrder> DeformableTransformType; // <CoordinateRepType, Dims, BSplineOrder>

ImageType::Pointer getCoronalAtlasSlice(int);
ImageType::Pointer rotateImage(ImageType::Pointer);
RigidTransformType::Pointer getRigidRegistrationTransform(ImageType::Pointer, ImageType::Pointer);
DeformableTransformType::Pointer getDeformableRegistrationTransform(ImageType::Pointer, ImageType::Pointer, ImageType::Pointer, ImageType::Pointer);
double degreesToRadians(double);
void writeImage(ImageType::Pointer, const char *);
void displayRegistrationResults(itk::OptimizerParameters<double>, const unsigned int, const double);
void debugOut(const char *);

// Observer class, to output the progress of the registration
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

int main(int argc, char *argv[]){
    //Check args
    if (argc < 6) {
            std::cerr << "Missing Parameters " << std::endl;
            std::cerr << "Usage: " << argv[0];
            std::cerr << " fixedImagePath movingImagePath";
            std::cerr << " fixedMaskPath movingMaskPath";
            std::cerr << " outputPath";
            return EXIT_FAILURE;
    }

    // Declare Image, filter types ect
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::ImageFileWriter<ImageType> WriterType;

    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
    

    // Get corresponding atlas slice
    //ImageType::Pointer atlasSlice = getCoronalAtlasSlice(atoi(argv[2]));

    ReaderType::Pointer atlasReader = ReaderType::New();
    atlasReader->SetFileName(argv[2]);
    atlasReader->Update();
    ImageType::Pointer atlasSlice = atlasReader->GetOutput();

    // Load the input image
    ReaderType::Pointer inputReader = ReaderType::New();
    inputReader->SetFileName(argv[1]);
    inputReader->Update();
    ImageType::Pointer inputImage = inputReader->GetOutput();

    // Load the masks
    ReaderType::Pointer fixedMaskReader = ReaderType::New();
    fixedMaskReader->SetFileName(argv[3]);
    fixedMaskReader->Update();
    ImageType::Pointer fixedMask = fixedMaskReader->GetOutput();

    ReaderType::Pointer movingMaskReader = ReaderType::New();
    movingMaskReader->SetFileName(argv[4]);
    movingMaskReader->Update();
    ImageType::Pointer movingMask = movingMaskReader->GetOutput();

    // Set the input spacing and atlas direction
    ImageType::SpacingType spacing;
    spacing[0] = 25.7764;
    spacing[1] = 25.7764;
    inputImage->SetSpacing(spacing);
    atlasSlice->SetSpacing(spacing);
    atlasSlice->SetDirection(inputImage->GetDirection());

    // Save the slice locally, until I feel more confident
    writeImage(atlasSlice, "/home/sam/Desktop/extractedSlice.jpg");
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
    writeImage(rigidResampler->GetOutput(), "/home/sam/afterRigidRegistration.jpg");

    // Compute a deformable registration transform using the resampled atlas slice and the input image
    DeformableTransformType::Pointer deformableTransform = getDeformableRegistrationTransform(rigidResampler->GetOutput(), inputImage, movingMask, fixedMask); // (movingImage, fixedImage)

    ResampleFilterType::Pointer deformableResampler = ResampleFilterType::New();
    deformableResampler->SetTransform(deformableTransform);
    deformableResampler->SetInput(rigidResampler->GetOutput());

    // This should really be done with a function...
    deformableResampler->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    deformableResampler->SetOutputOrigin(inputImage->GetOrigin());
    deformableResampler->SetOutputSpacing(inputImage->GetSpacing());
    deformableResampler->SetOutputDirection(inputImage->GetDirection());
    deformableResampler->SetDefaultPixelValue(0);

    // Finally, connect the file writer to write the output file
    WriterType::Pointer outputWriter = WriterType::New();
    outputWriter->SetFileName(argv[3]);
    outputWriter->SetInput(deformableResampler->GetOutput());
    outputWriter->Update();

    return EXIT_SUCCESS;
}

ImageType::Pointer getCoronalAtlasSlice(int sliceDepth) {
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

DeformableTransformType::Pointer getDeformableRegistrationTransform(ImageType::Pointer movingImage, ImageType::Pointer fixedImage, ImageType::Pointer movingMask, ImageType::Pointer fixedMask){
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType> NormalizedCorrelationMetricType;
    typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> MeansSquaredMetricType;
    typedef itk::CombinedImageToImageMetricAdaptor<ImageType, ImageType> CombinedMetricType;
    typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
    typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
    typedef DeformableTransformType::ParametersType ParametersType;

    // Instantiate the metric, optimizer, interpolator, transform and registration objects
    CombinedMetricType::Pointer combinedMetric = CombinedMetricType::New();
    OptimizerType::Pointer optimizer = OptimizerType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    RegistrationType::Pointer registration = RegistrationType::New();
    DeformableTransformType::Pointer transform = DeformableTransformType::New();

    // Connect everything to the registration object
    registration->SetMetric(combinedMetric);
    registration->SetOptimizer(optimizer);
    registration->SetInterpolator(interpolator);
    registration->SetTransform(transform);

    // Connect the inputs
    registration->SetFixedImage(fixedImage);
    registration->SetMovingImage(movingImage);
    registration->SetFixedImageRegion(fixedImage->GetBufferedRegion());

    // Calculate image physical dimensions and meshSize
    DeformableTransformType::PhysicalDimensionsType fixedPhysicalDimensions;
    DeformableTransformType::MeshSizeType meshSize;
    DeformableTransformType::OriginType fixedOrigin = fixedImage->GetOrigin();
    ImageType::SizeType fixedImageSize = fixedImage->GetLargestPossibleRegion().GetSize();
    unsigned int numberOfGridNodesInOneDimension = 12;

    for (int i = 0; i < 2; i++) {
        fixedPhysicalDimensions[i] = (fixedImageSize[i] - 2) * fixedImage->GetSpacing()[i];
    }
    meshSize.Fill(numberOfGridNodesInOneDimension - BSplineOrder);

    // Define the transform domain
    transform->SetTransformDomainOrigin(fixedOrigin);
    transform->SetTransformDomainPhysicalDimensions(fixedPhysicalDimensions);
    transform->SetTransformDomainMeshSize(meshSize);
    transform->SetTransformDomainDirection(fixedImage->GetDirection());

    // Configure the two image metrics
    // The normalized correlation metric inputs will be taken directly from the registration
    NormalizedCorrelationMetricType::Pointer normCorrMetric = NormalizedCorrelationMetricType::New();

    // The means squared metric inputs will come from segmented masks
    MeansSquaredMetricType::Pointer meanSqMetric = MeansSquaredMetricType::New();
    InterpolatorType::Pointer meanSqInterpolator = InterpolatorType::New(); // Not sure if I need a separate interpolator here
    meanSqMetric->SetFixedImage(fixedMask);
    meanSqMetric->SetMovingImage(movingMask);
    meanSqMetric->SetFixedImageRegion(fixedMask->GetBufferedRegion());
    meanSqMetric->SetInterpolator(meanSqInterpolator);

    combinedMetric->SetMetric(0, normCorrMetric, 1.0);
    combinedMetric->SetMetric(1, meanSqMetric, 0.001);

    // Set the initial transform parameters
    const unsigned int numberOfParameters = transform->GetNumberOfParameters();
    ParametersType parameters(numberOfParameters);
    parameters.Fill(0.0);
    transform->SetParameters(parameters);
    registration->SetInitialTransformParameters(transform->GetParameters());

    // Specify the optimizer parameters
    optimizer->MaximizeOff();
    optimizer->SetMaximumStepLength(10.0);
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