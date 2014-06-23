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

#include "itkPermuteAxesImageFilter.h"

#include "itkExceptionObject.h"

#include "itkCommand.h"

const char * atlasPath = "/home/sam/Pictures/allenReferenceAtlas_mouseCoronal/atlasVolume/atlasVolume.mhd";

itk::Image<unsigned char, 2>::Pointer getCoronalAtlasSlice(int);
itk::Image<unsigned char, 2>::Pointer rotateImage(itk::Image<unsigned char, 2>::Pointer);
double degreesToRadians(double);
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

int main(int argc, char *argv[]){
    //Check args
    if (argc < 4) {
            std::cerr << "Missing Parameters " << std::endl;
            std::cerr << "Usage: " << argv[0];
            std::cerr << " sliceToRegisterPath anteriorPosteriorCoorinate";
            std::cerr << " outputPath";
            return EXIT_FAILURE;
    }

    // Declare Image, filter types ect
    typedef itk::Image<unsigned char, 2> AtlasSliceType;
    typedef itk::Image<unsigned char, 2> InputImageType;

    typedef itk::ImageFileReader<InputImageType> InputReaderType;
    typedef itk::ImageFileWriter<AtlasSliceType> SliceWriterType;

    typedef itk::CenteredRigid2DTransform<double> TransformType;
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef itk::MeanSquaresImageToImageMetric<InputImageType, AtlasSliceType> MetricType;
    typedef itk::LinearInterpolateImageFunction<AtlasSliceType, double> InterpolatorType;
    typedef itk::ImageRegistrationMethod<InputImageType, AtlasSliceType> RegistrationType;
    typedef itk::ResampleImageFilter<AtlasSliceType, InputImageType> ResampleFilterType;
    
    typedef itk::CenteredTransformInitializer<TransformType, InputImageType, AtlasSliceType> TransformInitializerType;

    // Get corresponding atlas slice
    AtlasSliceType::Pointer atlasSlice = getCoronalAtlasSlice(atoi(argv[2]));

    // Save the slice locally, until I feel more confident
    SliceWriterType::Pointer debugWriter = SliceWriterType::New();
    debugWriter->SetFileName("extractedSlice.jpg");
    debugWriter->SetInput(atlasSlice);
    debugWriter->Update();
    
    // --- Rough Rigid Registration ---

    // Set up our input reader
    InputReaderType::Pointer inputReader = InputReaderType::New();
    inputReader->SetFileName(argv[1]);
    inputReader->Update();
    InputImageType::Pointer inputImage = inputReader->GetOutput();

    // Instantiate the metric, optimizer, interpolator and registration objects
    MetricType::Pointer         metric        = MetricType::New();
    OptimizerType::Pointer      optimizer     = OptimizerType::New();
    InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
    TransformType::Pointer      transform     = TransformType::New();
    RegistrationType::Pointer   registration  = RegistrationType::New();

    // Set up the registration
    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetInterpolator(interpolator);
    registration->SetTransform(transform);

    registration->SetInitialTransformParameters(transform->GetParameters());

    // Set the inputs
    registration->SetFixedImage(inputImage);
    registration->SetMovingImage(atlasSlice);

    registration->SetFixedImageRegion(atlasSlice->GetBufferedRegion());

    // Create an initializer and connect it to the input images and the transform
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetFixedImage(inputImage);
    initializer->SetMovingImage(atlasSlice);
    initializer->SetTransform(transform);

    // MomentsOn() sets the initializer to center of mass mode
    initializer->MomentsOn();
    initializer->InitializeTransform();

    // Set up the optimizer
    //  The Optimizer is responsible for sustaining and terminating registration iteration
    // I'm not entirely sure what the purpose of scales is
    OptimizerType::ScalesType optimizerScales(transform->GetNumberOfParameters());
    const double translationScale = 1.0/1000.0;

    optimizerScales[0] = 1.0;
    for (int i = 1; i <= 4; i++) {
        optimizerScales[i] = translationScale;
    }
    optimizer->SetScales(optimizerScales);

    // These metrics determine when the optimizer will terminate registration
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
    } catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    // Display the final registration transformation parameters
    OptimizerType::ParametersType registrationParameters = registration->GetLastTransformParameters();
    displayRegistrationResults(registrationParameters, optimizer->GetCurrentIteration(), optimizer->GetValue());

    // Apply the computed rigid registration to the slice
    TransformType::Pointer registrationTransform = TransformType::New();
    registrationTransform->SetParameters(registrationParameters);
    registrationTransform->SetFixedParameters(transform->GetFixedParameters());

    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetTransform(registrationTransform);
    resampler->SetInput(atlasSlice);

    resampler->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputOrigin(inputImage->GetOrigin());
    resampler->SetOutputSpacing(inputImage->GetSpacing());
    resampler->SetOutputDirection(inputImage->GetDirection());
    resampler->SetDefaultPixelValue(0);

    // Finally, connect the file writer to write the output file
    SliceWriterType::Pointer outputWriter = SliceWriterType::New();
    outputWriter->SetFileName(argv[3]);
    outputWriter->SetInput(resampler->GetOutput());
    outputWriter->Update();

    return EXIT_SUCCESS;
}

itk::Image<unsigned char, 2>::Pointer getCoronalAtlasSlice(int sliceDepth) {
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

itk::Image<unsigned char, 2>::Pointer rotateImage(itk::Image<unsigned char, 2>::Pointer inputImage) {
    typedef itk::Image<unsigned char, 2> ImageType;
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