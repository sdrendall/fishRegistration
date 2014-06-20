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
#include "itkAffineTransform.h"

#include "itkExceptionObject.h"

#include "itkCommand.h"

const char * atlasPath = "/home/sam/Pictures/allenReferenceAtlas_mouseCoronal/atlasVolume/atlasVolume.mhd";

itk::Image<unsigned char, 2>::Pointer getCoronalAtlasSlice(int);
itk::Image<unsigned char, 2>::Pointer rotateImage(itk::Image<unsigned char, 2>::Pointer);
double degreesToRadians(double);
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

    // Set the inputs
    registration->SetFixedImage(inputReader->GetOutput());
    registration->SetMovingImage(atlasSlice);

    registration->SetFixedImageRegion(atlasSlice->GetBufferedRegion());

    // Create an initializer and connect it to the input images and the transform
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetFixedImage(inputReader->GetOutput());
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

    // Finally, connect the file writer to write the output file
    SliceWriterType::Pointer outputWriter = SliceWriterType::New();
    outputWriter->SetFileName(argv[3]);
    outputWriter->SetInput(registration->GetOutput());

    // Begin Registration by calling Update() on the outputWriter object
    try {
    outputWriter->Update();
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    } catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    // Display final results
    OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();

    const double finalAngle           = finalParameters[0];
    const double finalRotationCenterX = finalParameters[1];
    const double finalRotationCenterY = finalParameters[2];
    const double finalTranslationX    = finalParameters[3];
    const double finalTranslationY    = finalParameters[4];
    
    const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
    const double bestValue = optimizer->GetValue();
    
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
    typedef itk::ResampleImageFilter<ImageType, ImageType> RotationFilterType;
    typedef itk::AffineTransform<double, 2> TransformType;
    typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;

    const float rotationAngle = 90.0;

    RotationFilterType::Pointer rotationFilter = RotationFilterType::New();
    rotationFilter->SetInterpolator(InterpolatorType::New());
    rotationFilter->SetDefaultPixelValue(0);

    // Copy parameters from the input image
    const ImageType::SpacingType & spacing = inputImage->GetSpacing();
    const ImageType::PointType & origin = inputImage->GetOrigin();
    ImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();

    // Set Output Parameters
    rotationFilter->SetOutputOrigin(origin);
    rotationFilter->SetOutputSpacing(spacing);
    rotationFilter->SetOutputDirection(inputImage->GetDirection());
    rotationFilter->SetSize(size);

    // Specify the set of tranformations to be carried out by the rotationFilter
    // Translate the image so that it's centered on it's physical origin
    TransformType::Pointer transform = TransformType::New();
    TransformType::OutputVectorType translation1;

    const double imageCenterX = origin[0] + spacing[0] * size[0] / 2.0;
    const double imageCenterY = origin[1] + spacing[1] * size[1] / 2.0;

    translation1[0] = -imageCenterX;
    translation1[1] = -imageCenterY;

    transform->Translate(translation1);

    // Rotate the image around it's physical origin
    const double angle = degreesToRadians(rotationAngle);
    transform->Rotate2D(-angle, false);

    // Translate back to original position
    TransformType::OutputVectorType translation2;

    translation2[0] = imageCenterX;
    translation2[1] = imageCenterY;

    transform->Translate(translation2, false);

    // Now that the transform is specified, connect it to the rotationFilter.
    rotationFilter->SetTransform(transform);

    // Perform the rotation and return the result
    rotationFilter->Update();
    ImageType::Pointer rotatedImage = rotationFilter->GetOutput();
    return rotatedImage;
}


double degreesToRadians(double degrees) {
    const double conversionFactor = vcl_atan(1.0)/45.0;
    double radians = degrees * conversionFactor;
    return radians;
}

void debugOut(const char * msg) {
    std::cout << "[Debug] " << msg << std::endl;
}