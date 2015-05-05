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
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkBSplineTransform.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkBSplineTransformInitializer.h"

#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkCastImageFilter.h"

// #include "itkGridImageSource.h"
#include "itkGridForwardWarpImageFilter.h"
#include "itkTransformToDisplacementFieldFilter.h"

#include "itkPermuteAxesImageFilter.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkCommand.h"

// Should I make this adjustable?
const unsigned int BSplineOrder = 3;

// Types
typedef float ImagePixelType;
typedef itk::Image<ImagePixelType, 2> ImageType;

typedef itk::CenteredRigid2DTransform<double> RigidTransformType;
typedef itk::BSplineTransform<double, 2, BSplineOrder> BSplineTransformType; // <CoordinateRepType, Dims, BSplineOrder>

typedef itk::Image<unsigned char, 2> EightBitImageType;
typedef itk::CastImageFilter<ImageType, EightBitImageType> InternalToEightBitCasterType;

// Function Declarations
RigidTransformType::Pointer compute_affine_transform(ImageType::Pointer, ImageType::Pointer);
BSplineTransformType::Pointer compute_bSpline_transform(ImageType::Pointer, ImageType::Pointer, const char *);

void write_image(ImageType::Pointer, const char *);
void display_registration_results(itk::OptimizerParameters<double>, const unsigned int, const double);
void debug_out(const char *);

// Function Templates for atlas manipulations

template<typename ATLAS_PIXEL_TYPE, typename INTERPOLATOR_TYPE>
void apply_transform_to_atlas(RigidTransformType::Pointer rigid_transform, 
                              BSplineTransformType::Pointer deformable_transform,
                              INTERPOLATOR_TYPE interpolator,
                              ImageType::Pointer input_image, 
                              const char * output_path, 
                              const char * atlas_path, 
                              int slice_index);

template<typename ATLAS_PIXEL_TYPE>
typename itk::Image<ATLAS_PIXEL_TYPE, 2>::Pointer get_atlas_slice(int slice_index, const char * atlas_path);

template<typename IMAGE_POINTER_TYPE>
IMAGE_POINTER_TYPE rotateImage(IMAGE_POINTER_TYPE input_image);


template<typename ATLAS_PIXEL_TYPE, typename INTERPOLATOR_TYPE>
void apply_transform_to_atlas(RigidTransformType::Pointer rigid_transform, 
        BSplineTransformType::Pointer deformable_transform, INTERPOLATOR_TYPE interpolator, ImageType::Pointer input_image, 
        const char * output_path, const char * atlas_path, int slice_index) {

    typedef itk::Image<ATLAS_PIXEL_TYPE, 2> AtlasSliceType;
    typedef itk::ResampleImageFilter<AtlasSliceType, AtlasSliceType> ImageResamplerType;
    typedef itk::ImageFileWriter<AtlasSliceType> ImageWriterType;

    typename AtlasSliceType::Pointer atlas_slice = get_atlas_slice<ATLAS_PIXEL_TYPE>(slice_index, atlas_path);
    atlas_slice->SetDirection(input_image->GetDirection());

    // Apply the computed rigid transform to the atlas slice
    typename ImageResamplerType::Pointer rigid_resampler = ImageResamplerType::New();
    rigid_resampler->SetTransform(rigid_transform);
    rigid_resampler->SetInput(atlas_slice);

    rigid_resampler->SetSize(input_image->GetLargestPossibleRegion().GetSize());
    rigid_resampler->SetOutputOrigin(input_image->GetOrigin());
    rigid_resampler->SetOutputSpacing(input_image->GetSpacing());
    rigid_resampler->SetOutputDirection(input_image->GetDirection());
    rigid_resampler->SetDefaultPixelValue(0);

    typename ImageResamplerType::Pointer deformable_resampler = ImageResamplerType::New();
    deformable_resampler->SetTransform(deformable_transform);
    deformable_resampler->SetInput(rigid_resampler->GetOutput());

    // Configure the resampler to apply the BSpline Transform
    deformable_resampler->SetSize(input_image->GetLargestPossibleRegion().GetSize());
    deformable_resampler->SetOutputOrigin(input_image->GetOrigin());
    deformable_resampler->SetOutputSpacing(input_image->GetSpacing());
    deformable_resampler->SetOutputDirection(input_image->GetDirection());
    deformable_resampler->SetDefaultPixelValue(0);

    // Link the interpolator to the resamplers
    rigid_resampler->SetInterpolator(interpolator);
    deformable_resampler->SetInterpolator(interpolator);

    // Write the registered reference image to the specified filepath
    typename ImageWriterType::Pointer output_writer = ImageWriterType::New();
    output_writer->SetFileName(output_path);
    output_writer->SetInput(deformable_resampler->GetOutput());
    output_writer->Update();
}


// template<typename ATLAS_PIXEL_TYPE, typename INTERPOLATOR_TYPE>


template<typename ATLAS_PIXEL_TYPE>
typename itk::Image<ATLAS_PIXEL_TYPE, 2>::Pointer get_atlas_slice(int slice_index, const char * atlas_path) {
    // Declare 3D and 2D image types
    typedef itk::Image<ATLAS_PIXEL_TYPE, 3> AtlasImageType;
    typedef itk::Image<ATLAS_PIXEL_TYPE, 2> SliceImageType;
    
    // Declare Readers and Writers for 3D images and 2D images respectively
    typedef itk::ImageFileReader<AtlasImageType> AtlasReaderType;
    typedef itk::ExtractImageFilter<AtlasImageType, SliceImageType> SliceAtlasFilterType;

    // Collapsing along the x axis will yield coronal sections
    const int axis_to_collapse = 0;
    
    // Create a reader, to read the input image
    typename AtlasReaderType::Pointer reader = AtlasReaderType::New();
    reader->SetFileName(atlas_path);

    // Update, to load the image so that size information can be ascertained
    reader->Update();

    // Get the whole atlas region
    typename AtlasImageType::Pointer atlas = reader->GetOutput();
    typename AtlasImageType::RegionType entire_atlas_region = atlas->GetLargestPossibleRegion();

    // Determine Slice Size
    typename AtlasImageType::SizeType input_size = entire_atlas_region.GetSize();
    typename AtlasImageType::SizeType slice_size = input_size;
    slice_size[axis_to_collapse] = 0;  // 0 tells the ExtractionFilter to return an Image without that dimension

    // Get Start Point
    typename AtlasImageType::IndexType slice_start_index = entire_atlas_region.GetIndex();
    slice_start_index.Fill(0);
    slice_start_index[axis_to_collapse] = slice_index;
    
    // Initialize a slice region
    typename AtlasImageType::RegionType sliceRegion(slice_start_index, slice_size);

    // Reference an extraction filter and connect it to the reader
    typename SliceAtlasFilterType::Pointer extraction_filter = SliceAtlasFilterType::New();
    extraction_filter->SetInput(atlas);
    extraction_filter->SetDirectionCollapseToIdentity();

    // Extract the slice
    extraction_filter->SetExtractionRegion(sliceRegion);
    extraction_filter->Update();
    typename SliceImageType::Pointer atlas_slice = extraction_filter->GetOutput();
    //atlas_slice = rotateImage(atlas_slice);
    return atlas_slice;
}


template<typename IMAGE_POINTER_TYPE>
IMAGE_POINTER_TYPE rotateImage(IMAGE_POINTER_TYPE input_image) {
    typedef itk::PermuteAxesImageFilter<IMAGE_POINTER_TYPE> PermuteAxesFilterType;

    // Create a permutation filter
    typename PermuteAxesFilterType::Pointer permute_image_filter = PermuteAxesFilterType::New();

    // Define the axes permutation order
    itk::FixedArray<unsigned int, 2> order;
    order[0] = 1;
    order[1] = 0;
    
    // Permute the axes
    permute_image_filter->SetInput(input_image);
    permute_image_filter->SetOrder(order);
    permute_image_filter->Update();

    // Return the permuted image
    IMAGE_POINTER_TYPE output_image = permute_image_filter->GetOutput();
    return output_image;
}


template<typename TRANSFORM_TYPE, typename IMAGE_TYPE>
int saveTransformAsDisplacementGrid(const char * save_path, TRANSFORM_TYPE transform, IMAGE_TYPE reference_image) {
    // Display or save the displacement field here
    // Create a displacement field based on the transform
    typedef itk::Vector<double, 2> DisplacementVectorType;
    typedef itk::Image<DisplacementVectorType, 2> DisplacementFieldType;
    typedef itk::TransformToDisplacementFieldFilter<DisplacementFieldType, double> DisplacementFieldGeneratorType;

    DisplacementFieldGeneratorType::Pointer displacement_field_generator = DisplacementFieldGeneratorType::New();
    displacement_field_generator->UseReferenceImageOn();
    displacement_field_generator->SetReferenceImage(reference_image);
    displacement_field_generator->SetTransform(transform);

    try {
        displacement_field_generator->Update();
    } catch ( itk::ExceptionObject & err ) {
        std::cout << "Could not generate deformation field from the given transform!" << std::endl;
        std::cout << "A displacement field grid will not be saved at " << save_path << std::endl;
        std::cout << err << std::endl;
        return 1;
    }

    typedef itk::GridForwardWarpImageFilter<DisplacementFieldType, EightBitImageType> GridFilterType;
    typename GridFilterType::Pointer grid_warper = GridFilterType::New();
    grid_warper->SetInput(displacement_field_generator->GetOutput());
    grid_warper->SetForegroundValue(itk::NumericTraits<unsigned char>::max());  // Foreground should be the maximum value for pixels in the grid image, which are uchars

    typedef itk::ImageFileWriter<EightBitImageType> GridImageWriterType;
    GridImageWriterType::Pointer grid_writer = GridImageWriterType::New();
    grid_writer->SetInput(grid_warper->GetOutput());
    grid_writer->SetFileName(save_path);

    try {
        grid_writer->Update();
    } catch (itk::ExceptionObject & err) {
        std::cout << "Failed to save deformation grid to " << save_path << std::endl;
        std::cout << err << std::endl;
        return 1;
    }

    return 0;
}


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
    if (argc < 8) {
            std::cerr << "Missing Parameters " << std::endl;
            std::cerr << "Usage: " << argv[0];
            std::cerr << " inputImagePath sliceIndex";
            std::cerr << " referenceImageOutputPath";
            std::cerr << " annotationImageOutputPath";
            std::cerr << " hemisphereImageOutputPath";
            std::cerr << " registrationLogPath";
            std::cerr << " atlasReferencePath";
            std::cerr << " atlasAnnotationsPath";
            std::cerr << " atlasHemispheresPath";
            return EXIT_FAILURE;
    }

    // Input paths
    const char * input_image_path = argv[1];
    const char * reference_image_output_path = argv[3];
    const char * annotation_image_output_path = argv[4];
    const char * hemisphere_image_output_path = argv[5];
    const char * registration_log_path = argv[6];
    const char * atlas_reference_path = argv[7]; 
    const char * atlas_annotations_path = argv[8]; 
    const char * atlas_hemispheres_path = argv[9]; 

    // Declare Image, filter types ect
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::ImageFileWriter<EightBitImageType> WriterType;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;


    // Get corresponding atlas reference slice
    int slice_index = atoi(argv[2]);
    ImageType::Pointer atlas_slice = get_atlas_slice<ImagePixelType>(slice_index, atlas_reference_path);
    
    // Load the input image
    ReaderType::Pointer input_reader = ReaderType::New();
    input_reader->SetFileName(input_image_path);
    input_reader->Update();
    ImageType::Pointer input_image = input_reader->GetOutput();

    // Set the input spacing and atlas direction
    ImageType::SpacingType spacing;
    spacing[0] = 25;
    spacing[1] = 25;
    input_image->SetSpacing(spacing);
    atlas_slice->SetDirection(input_image->GetDirection());

    RigidTransformType::Pointer rigid_transform = compute_affine_transform(input_image, atlas_slice);
    
    saveTransformAsDisplacementGrid("/home/sam/Desktop/rigid_displacement_grid.png", rigid_transform, input_image);

    // Apply the computed rigid transform to the atlas slice
    ResampleFilterType::Pointer rigidResampler = ResampleFilterType::New();
    rigidResampler->SetTransform(rigid_transform);
    rigidResampler->SetInput(atlas_slice);

    rigidResampler->SetSize(input_image->GetLargestPossibleRegion().GetSize());
    rigidResampler->SetOutputOrigin(input_image->GetOrigin());
    rigidResampler->SetOutputSpacing(input_image->GetSpacing());
    rigidResampler->SetOutputDirection(input_image->GetDirection());
    rigidResampler->SetDefaultPixelValue(0);

    // Compute the Bspline transform mapping the atlas slice to the input image
    BSplineTransformType::Pointer deformable_transform = compute_bSpline_transform(rigidResampler->GetOutput(), input_image, registration_log_path); // (movingImage, fixedImage, logPath)

    saveTransformAsDisplacementGrid("/home/sam/Desktop/deformable_displacement_grid.png", deformable_transform, input_image);

    // Apply the computed transforms to the corresponding reference and annotation images, and save the results
    // applyAtlas(rigid_transform, deformable_transform, interpolator, input_image, outputName, atlasPath, slice_index)

    // Different interpolators are needed for the different atlas images to account for the differences in bit-depth.
    typedef unsigned char ReferencePixelType;
    typedef itk::Image<ReferencePixelType, 2> ReferenceImageType;
    typedef itk::BSplineInterpolateImageFunction<ReferenceImageType, double, double> ReferenceInterpolatorType;
    ReferenceInterpolatorType::Pointer reference_interpolator = ReferenceInterpolatorType::New();
    reference_interpolator->SetSplineOrder(3);

    // Apply transforms to the atlas reference image
    apply_transform_to_atlas<ReferencePixelType, ReferenceInterpolatorType::Pointer>(rigid_transform,
                                                                            deformable_transform, 
                                                                            reference_interpolator,
                                                                            input_image,
                                                                            reference_image_output_path,
                                                                            atlas_reference_path,
                                                                            slice_index);

    typedef unsigned char HemispherePixelType;
    typedef itk::Image<HemispherePixelType, 2> HemisphereImageType;
    typedef itk::NearestNeighborInterpolateImageFunction<HemisphereImageType, double> HemisphereInterpolatorType;
    HemisphereInterpolatorType::Pointer hemisphere_interpolator = HemisphereInterpolatorType::New();

    // Apply transforms to the hemisphere atlas images
    apply_transform_to_atlas<HemispherePixelType, HemisphereInterpolatorType::Pointer>(rigid_transform,
                                                                              deformable_transform,
                                                                              hemisphere_interpolator,
                                                                              input_image,
                                                                              hemisphere_image_output_path,
                                                                              atlas_hemispheres_path,
                                                                              slice_index);

    typedef double AnnotationPixelType;
    typedef itk::Image<AnnotationPixelType, 2> AnnotationImageType;
    typedef itk::LinearInterpolateImageFunction<AnnotationImageType, double> AnnotationInterpolatorType;
    AnnotationInterpolatorType::Pointer annotation_interpolator = AnnotationInterpolatorType::New();

    // Apply transforms to the annotated atlas images
    apply_transform_to_atlas<AnnotationPixelType, AnnotationInterpolatorType::Pointer>(rigid_transform, 
                                                                              deformable_transform,
                                                                              annotation_interpolator,
                                                                              input_image,
                                                                              annotation_image_output_path,
                                                                              atlas_annotations_path,
                                                                              slice_index);

    return EXIT_SUCCESS;
}


RigidTransformType::Pointer compute_affine_transform(ImageType::Pointer inputImage, ImageType::Pointer atlasSlice) {
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
    display_registration_results(registrationParameters, optimizer->GetCurrentIteration(), optimizer->GetValue());

    // // return the transform computed by the registration
    // RigidTransformType::Pointer registrationTransform = RigidTransformType::New();
    // registrationTransform->SetParameters(registrationParameters);
    // registrationTransform->SetFixedParameters(transform->GetFixedParameters());
    transform->SetParameters(registration->GetLastTransformParameters());

    return transform;
}


BSplineTransformType::Pointer compute_bSpline_transform(ImageType::Pointer movingImage, ImageType::Pointer fixedImage, const char * logPath){
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
    unsigned int numberOfGridNodesInOneDimension = 5;

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


void write_image(ImageType::Pointer im, const char * path) {
    typedef itk::ImageFileWriter<EightBitImageType> WriterType;

    InternalToEightBitCasterType::Pointer caster = InternalToEightBitCasterType::New();
    caster->SetInput(im);

    WriterType::Pointer fileWriter = WriterType::New();
    fileWriter->SetFileName(path);
    fileWriter->SetInput(caster->GetOutput());
    fileWriter->Update();
}


void debug_out(const char * msg) {
    std::cout << "[Debug] " << msg << std::endl;
}


void display_registration_results(itk::OptimizerParameters<double> finalParameters, const unsigned int numberOfIterations, const double bestValue) {
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
