#include <stdio.h>
#include <iostream>
#include <fstream>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"

#include "itkImageRegistrationMethodv4.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkBSplineTransform.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkBSplineTransformInitializer.h"

#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkCastImageFilter.h"


#include "itkPermuteAxesImageFilter.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkCommand.h"

// Globals (TODO: Do this differently)
const char * atlasReferencePath = "/home/sr235/atlasVolume/atlasVolume.mhd";
const char * atlasAnnotationsPath = "/home/sr235/atlasVolume/annotation.mhd";
//const char * atlasReferencePath = "/home/sam/Dropbox/grayLab/allenReferenceAtlas_mouseCoronal/atlasVolume/atlasVolume.mhd";
//const char * atlasAnnotationsPath = "/home/sam/Dropbox/grayLab/allenReferenceAtlas_mouseCoronal/atlasVolume/annotation.mhd";
const unsigned int BSplineOrder = 3;

// Types
typedef float ImagePixelType;
typedef itk::Image<ImagePixelType, 2> ImageType;

typedef unsigned int AnnotationPixelType;
typedef itk::Image<AnnotationPixelType, 2> AnnotationImageType;

typedef itk::CenteredRigid2DTransform<double> RigidTransformType;
typedef itk::BSplineTransform<double, 2, BSplineOrder> BSplineTransformType; // <CoordinateRepType, Dims, BSplineOrder>

typedef itk::Image<unsigned char, 2> EightBitImageType;
typedef itk::CastImageFilter<ImageType, EightBitImageType> InternalToEightBitCasterType;

// Function Declarations
ImageType::Pointer getCoronalAtlasSlice(int, const char *);
AnnotationImageType::Pointer getCoronalAnnotationAtlasSlice(int, const char *);
ImageType::Pointer rotateImage(ImageType::Pointer);
AnnotationImageType::Pointer rotateAnnotationImage(AnnotationImageType::Pointer);
RigidTransformType::Pointer getRigidRegistrationTransform(ImageType::Pointer, ImageType::Pointer);
BSplineTransformType::Pointer computeBsplineTransform(ImageType::Pointer, ImageType::Pointer, const char *);

void applyTransformToAnnotations(RigidTransformType::Pointer, BSplineTransformType::Pointer, ImageType::Pointer, const char *, const char *, int);
void applyTransformToReference(RigidTransformType::Pointer, BSplineTransformType::Pointer, ImageType::Pointer, const char *, const char *, int);

void writeImage(ImageType::Pointer, const char *);
void displayRegistrationResults(itk::OptimizerParameters<double>, const unsigned int, const double);
void debugOut(const char *);

// Observer class, to output the progress of the registration
// This was copied directly from a rigid registration example, without modification
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

// Observer for the BSplineTransform Computation.
class BsplineTransformIterationUpdater : public itk::Command {
public:
    typedef BsplineTransformIterationUpdater Self;
    typedef itk::Command Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    itkNewMacro(Self);

private:
    std::fstream outputFile;
    const char * outputFilePath;

protected:
    BsplineTransformIterationUpdater() {};

    ~BsplineTransformIterationUpdater() {
        closeOutputFile();
    };

public:
    typedef itk::RegularStepGradientDescentOptimizerv4<double> OptimizerType;
    typedef const OptimizerType * OptimizerPointer;

    void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE {
        Execute( (const itk::Object *)caller, event);
    }
  
    void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE {
        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
        if( !(itk::IterationEvent().CheckEvent( &event )) ){
            return;
        }

        ensureOutputFile();
        outputFile << optimizer->GetCurrentIteration() << "\t";
        outputFile << optimizer->GetCurrentMetricValue() << std::endl;
    }

    void ensureOutputFile() {
        if (!outputFile.is_open()) {
            openOutputFile();
        }
    }
    void setOutputFilePath(const char * path) {
        outputFilePath = path;
    }

    void openOutputFile() {
        outputFile.open(outputFilePath, std::ios::out | std::ios::trunc);
    }

    void closeOutputFile() {
        outputFile.close();
    }

};

int main(int argc, char *argv[]){
    //Check args
    if (argc < 6) {
            std::cerr << "Missing Parameters " << std::endl;
            std::cerr << "Usage: " << argv[0];
            std::cerr << " fixedImage sliceIndex";
            std::cerr << " referenceImageOutputPath";
            std::cerr << " annotationImageOutputPath";
            std::cerr << " registrationLogPath";
            return EXIT_FAILURE;
    }

    // Declare Image, filter types ect
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::ImageFileWriter<EightBitImageType> WriterType;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;


    // Get corresponding atlas reference slice
    int sliceIndex = atoi(argv[2]);
    ImageType::Pointer atlasSlice = getCoronalAtlasSlice(sliceIndex, atlasReferencePath);
    
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

    RigidTransformType::Pointer rigidTransform = getRigidRegistrationTransform(inputImage, atlasSlice);
    
    // Apply the computed rigid transform to the atlas slice
    ResampleFilterType::Pointer rigidResampler = ResampleFilterType::New();
    rigidResampler->SetTransform(rigidTransform);
    rigidResampler->SetInput(atlasSlice);

    rigidResampler->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    rigidResampler->SetOutputOrigin(inputImage->GetOrigin());
    rigidResampler->SetOutputSpacing(inputImage->GetSpacing());
    rigidResampler->SetOutputDirection(inputImage->GetDirection());
    rigidResampler->SetDefaultPixelValue(0);

    // Compute the Bspline transform mapping the atlas slice to the input image
    BSplineTransformType::Pointer deformableTransform = computeBsplineTransform(rigidResampler->GetOutput(), inputImage, argv[5]); // (movingImage, fixedImage, logPath)

    // Apply the computed transforms to the corresponding reference and annotation images, and save the results
    // applyTransform(rigidTransform, deformableTransform, inputImage, outputName, atlasPath, sliceIndex)
    applyTransformToReference(rigidTransform, deformableTransform, inputImage, argv[3], atlasReferencePath, sliceIndex);
    applyTransformToAnnotations(rigidTransform, deformableTransform, inputImage, argv[4], atlasAnnotationsPath, sliceIndex);

    return EXIT_SUCCESS;
}


void applyTransformToReference(RigidTransformType::Pointer rigidTransform, 
        BSplineTransformType::Pointer deformableTransform, ImageType::Pointer inputImage, 
        const char * outputName, const char * atlasReferencePath, int sliceIndex) {

    typedef itk::ResampleImageFilter<ImageType, ImageType> ReferenceImageResamplerType;
    typedef itk::ImageFileWriter<EightBitImageType> EightBitImageWriterType;

    ImageType::Pointer referenceImage = getCoronalAtlasSlice(sliceIndex, atlasReferencePath);
    referenceImage->SetDirection(inputImage->GetDirection());

    // Apply the computed rigid transform to the atlas slice
    ReferenceImageResamplerType::Pointer rigidResampler = ReferenceImageResamplerType::New();
    rigidResampler->SetTransform(rigidTransform);
    rigidResampler->SetInput(referenceImage);

    rigidResampler->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    rigidResampler->SetOutputOrigin(inputImage->GetOrigin());
    rigidResampler->SetOutputSpacing(inputImage->GetSpacing());
    rigidResampler->SetOutputDirection(inputImage->GetDirection());
    rigidResampler->SetDefaultPixelValue(0);

    ReferenceImageResamplerType::Pointer deformableResampler = ReferenceImageResamplerType::New();
    deformableResampler->SetTransform(deformableTransform);
    deformableResampler->SetInput(rigidResampler->GetOutput());

    // Configure the resampler to apply the BSpline Transform
    deformableResampler->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    deformableResampler->SetOutputOrigin(inputImage->GetOrigin());
    deformableResampler->SetOutputSpacing(inputImage->GetSpacing());
    deformableResampler->SetOutputDirection(inputImage->GetDirection());
    deformableResampler->SetDefaultPixelValue(0);

    // Output images need to be cast from float to char (uint8)
    InternalToEightBitCasterType::Pointer caster = InternalToEightBitCasterType::New();
    caster->SetInput(deformableResampler->GetOutput());

    // Write the registered reference image to the specified filepath
    EightBitImageWriterType::Pointer outputWriter = EightBitImageWriterType::New();
    outputWriter->SetFileName(outputName);
    outputWriter->SetInput(caster->GetOutput());
    outputWriter->Update();
}

void applyTransformToAnnotations(RigidTransformType::Pointer rigidTransform, 
        BSplineTransformType::Pointer deformableTransform, ImageType::Pointer inputImage, 
        const char * outputName, const char * atlasAnnotationsPath, int sliceIndex) {

    typedef itk::ResampleImageFilter<AnnotationImageType, AnnotationImageType> AnnotationImageResamplerType;
    typedef itk::NearestNeighborInterpolateImageFunction<AnnotationImageType, double> NearestNeighborInterpolatorType;
    typedef itk::ImageFileWriter<AnnotationImageType> AnnotationImageWriterType;

    AnnotationImageType::Pointer annotationImage = getCoronalAnnotationAtlasSlice(sliceIndex, atlasAnnotationsPath);
    annotationImage->SetDirection(inputImage->GetDirection());

    // Apply the computed rigid transform to the atlas slice
    AnnotationImageResamplerType::Pointer rigidResampler = AnnotationImageResamplerType::New();
    rigidResampler->SetTransform(rigidTransform);
    rigidResampler->SetInput(annotationImage);

    rigidResampler->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    rigidResampler->SetOutputOrigin(inputImage->GetOrigin());
    rigidResampler->SetOutputSpacing(inputImage->GetSpacing());
    rigidResampler->SetOutputDirection(inputImage->GetDirection());
    rigidResampler->SetDefaultPixelValue(0);

    AnnotationImageResamplerType::Pointer deformableResampler = AnnotationImageResamplerType::New();
    deformableResampler->SetTransform(deformableTransform);
    deformableResampler->SetInput(rigidResampler->GetOutput());

    // Configure the resampler to apply the BSpline Transform
    deformableResampler->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    deformableResampler->SetOutputOrigin(inputImage->GetOrigin());
    deformableResampler->SetOutputSpacing(inputImage->GetSpacing());
    deformableResampler->SetOutputDirection(inputImage->GetDirection());
    deformableResampler->SetDefaultPixelValue(0);

    // Use nearest neighbor interpolation to resample the annotation images
    NearestNeighborInterpolatorType::Pointer nnInterpolator = NearestNeighborInterpolatorType::New();
    rigidResampler->SetInterpolator(nnInterpolator);
    deformableResampler->SetInterpolator(nnInterpolator);

    // Write the registered reference image to the specified filepath
    AnnotationImageWriterType::Pointer outputWriter = AnnotationImageWriterType::New();
    outputWriter->SetFileName(outputName);
    outputWriter->SetInput(deformableResampler->GetOutput());
    outputWriter->Update();
}


ImageType::Pointer getCoronalAtlasSlice(int sliceIndex, const char * atlasPath) {
    // Declare 3D and 2D image types
    typedef itk::Image<ImagePixelType, 3> AtlasImageType;
    typedef itk::Image<ImagePixelType, 2> SliceImageType;
    
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

    // Get Start Point
    AtlasImageType::IndexType sliceStartIndex = entireAtlasRegion.GetIndex();
    sliceStartIndex.Fill(0);
    sliceStartIndex[axisToCollapse] = sliceIndex;  // Switched to slice index, I'll change this later
    
    // Initialize a slice region
    AtlasImageType::RegionType sliceRegion(sliceStartIndex, sliceSize);

    // Reference an extraction filter and connect it to the reader
    SliceAtlasFilterType::Pointer extFilter = SliceAtlasFilterType::New();
    extFilter->SetInput(atlas);
    extFilter->SetDirectionCollapseToIdentity();

    // Extract the slice
    extFilter->SetExtractionRegion(sliceRegion);
    extFilter->Update();
    SliceImageType::Pointer atlasSlice = extFilter->GetOutput();
    //atlasSlice = rotateImage(atlasSlice);
    return atlasSlice;
}

// Identical to getCoronalAtlasSlice(), except that a 32-bit uint image is loaded
AnnotationImageType::Pointer getCoronalAnnotationAtlasSlice(int sliceIndex, const char * atlasPath) {
    // Declare 3D and 2D image types
    typedef itk::Image<AnnotationPixelType, 3> AtlasImageType;
    typedef itk::Image<AnnotationPixelType, 2> SliceImageType;
    
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

    // Get Start Point
    AtlasImageType::IndexType sliceStartIndex = entireAtlasRegion.GetIndex();
    sliceStartIndex.Fill(0);
    sliceStartIndex[axisToCollapse] = sliceIndex;  // Switched to slice index, I'll change this later
    
    // Initialize a slice region
    AtlasImageType::RegionType sliceRegion(sliceStartIndex, sliceSize);

    // Reference an extraction filter and connect it to the reader
    SliceAtlasFilterType::Pointer extFilter = SliceAtlasFilterType::New();
    extFilter->SetInput(atlas);
    extFilter->SetDirectionCollapseToIdentity();

    // Extract the slice
    extFilter->SetExtractionRegion(sliceRegion);
    extFilter->Update();
    SliceImageType::Pointer atlasSlice = extFilter->GetOutput();
    atlasSlice = rotateAnnotationImage(atlasSlice);
    return atlasSlice;
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
    optimizer->SetMaximumStepLength(0.04);
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


BSplineTransformType::Pointer computeBsplineTransform(ImageType::Pointer movingImage, ImageType::Pointer fixedImage, const char * logPath){
    typedef itk::Image<float, 2> InternalImageType;
    typedef itk::CastImageFilter<ImageType, InternalImageType> ExternalToInternalImageCasterType;
    typedef itk::RegularStepGradientDescentOptimizerv4<double> OptimizerType;
    typedef itk::MattesMutualInformationImageToImageMetricv4<InternalImageType, InternalImageType> MetricType;
    typedef itk::ImageRegistrationMethodv4<InternalImageType, InternalImageType, BSplineTransformType> RegistrationType;
    typedef BSplineTransformType::ParametersType ParametersType;
    typedef itk::BSplineTransformInitializer<BSplineTransformType, InternalImageType> BSplineTransformInitializerType;


    // Instantiate the metric, optimizer, interpolator, transform and registration objects
    MetricType::Pointer metric = MetricType::New();
    OptimizerType::Pointer optimizer = OptimizerType::New();
    RegistrationType::Pointer registration = RegistrationType::New();
    BSplineTransformType::Pointer transform = BSplineTransformType::New();
    BSplineTransformInitializerType::Pointer transformInitializer = BSplineTransformInitializerType::New();

    // itkRegistrationMethodv4 accepts images with floating point pixels
    // Casters must be created to convert the fixed and moving images to floating point
    ExternalToInternalImageCasterType::Pointer fixedCaster = ExternalToInternalImageCasterType::New();
    ExternalToInternalImageCasterType::Pointer movingCaster = ExternalToInternalImageCasterType::New();
    fixedCaster->SetInput(fixedImage);
    movingCaster->SetInput(movingImage);

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

    // Transform initial parameters are determined by the transformInitializer
    transformInitializer->SetTransform(transform);
    transformInitializer->SetImage(fixedImage);
    transformInitializer->SetTransformDomainMeshSize(meshSize);
    transformInitializer->InitializeTransform();

    transform->SetIdentity();

    // Set Metic Parameters
    metric->SetNumberOfHistogramBins(64);

    // Specify the optimizer parameters
    optimizer->SetMinimumStepLength(0.001);
    optimizer->SetRelaxationFactor(0.7);
    optimizer->SetNumberOfIterations(2000);
    optimizer->SetLearningRate(150);

    // Add an observer to the optimizer
    BsplineTransformIterationUpdater::Pointer observer = BsplineTransformIterationUpdater::New();
    observer->setOutputFilePath(logPath);
    optimizer->AddObserver(itk::IterationEvent(), observer);

    // Connect everything to the registration object
    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetInitialTransform(transform);

    // Connect the filters to the registration object
    registration->SetFixedImage(fixedImage);
    registration->SetMovingImage(movingImage);

    // Set Multi-Resolution Options
    // The shrink factor denotes to the factor by which the image will be downsized
    // The smoothing sigma determines the width of the gaussian kernel used to smooth the downsampled image
    const unsigned int numberOfLevels = 3;

    RegistrationType::ShrinkFactorsArrayType shrinkFactorPerLevel;
    shrinkFactorPerLevel.SetSize(numberOfLevels);
    shrinkFactorPerLevel[0] = 4;
    shrinkFactorPerLevel[1] = 2;
    shrinkFactorPerLevel[2] = 1;

    RegistrationType::SmoothingSigmasArrayType smoothingSigmaPerLevel;
    smoothingSigmaPerLevel.SetSize(numberOfLevels);
    smoothingSigmaPerLevel[0] = 4;
    smoothingSigmaPerLevel[1] = 2;
    smoothingSigmaPerLevel[2] = 0;

    registration->SetNumberOfLevels(numberOfLevels);
    registration->SetShrinkFactorsPerLevel(shrinkFactorPerLevel);
    registration->SetSmoothingSigmasPerLevel(smoothingSigmaPerLevel);

    // Prepare time and memory probes
    //itk::TimeProbesCollectorBase chronometer;
    //itk::MemoryProbesCollectorBase memorymeter;
    
    // Start Registration
    try {
        //memorymeter.Start("Registration");
        //chronometer.Start("Registration");

        registration->Update();

        //chronometer.Stop("Registration");
        //memorymeter.Stop("Registration");

        std::cout << "Optimizer stop condition = "
                  << registration->GetOptimizer()->GetStopConditionDescription()
                  << std::endl;
    } catch ( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        throw 1337;
    }

    return transform;
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


AnnotationImageType::Pointer rotateAnnotationImage(AnnotationImageType::Pointer inputImage) {
    typedef itk::PermuteAxesImageFilter<AnnotationImageType> PermuteAxesFilterType;

    // Create a permutation filter
    PermuteAxesFilterType::Pointer permutationFilter = PermuteAxesFilterType::New();

    // Define the axes permutation order
    itk::FixedArray<unsigned int, 2> order;
    order[0] = 1;
    order[1] = 0;
    
    // Reorder the axes
    permutationFilter->SetInput(inputImage);
    permutationFilter->SetOrder(order);
    permutationFilter->Update();

    // Return the permuted image
    AnnotationImageType::Pointer outputImage = permutationFilter->GetOutput();
    return outputImage;
}


void writeImage(ImageType::Pointer im, const char * path) {
    typedef itk::ImageFileWriter<EightBitImageType> WriterType;

    InternalToEightBitCasterType::Pointer caster = InternalToEightBitCasterType::New();
    caster->SetInput(im);

    WriterType::Pointer fileWriter = WriterType::New();
    fileWriter->SetFileName(path);
    fileWriter->SetInput(caster->GetOutput());
    fileWriter->Update();
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
