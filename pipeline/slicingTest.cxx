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

typedef itk::Image<unsigned char, 2> ImageType;

ImageType::Pointer getCoronalAtlasSlice(int, const char *);
ImageType::Pointer rotateImage(ImageType::Pointer);
void writeImage(ImageType::Pointer, const char *);
void debugOut(const char *);

int main(int argc, char *argv[]) {
    debugOut("Checking input arguments...");
    int axisToCollapse, coordinateOnAxis;
    // Check arguments
    if (argc < 4) {
            std::cerr << "Missing Parameters " << std::endl;
            std::cerr << "Usage: " << argv[0];
            std::cerr << " atlasOutputPath annotationOutputPath coordinateOnAxis";
            std::cerr << " [axisToCollapse] ";
            return EXIT_FAILURE;
    } else if (argc >= 5) {
        axisToCollapse = atoi(argv[4]);
    } else {
        axisToCollapse = 0;
    }

    debugOut("Slicing Atlas");
    ImageType::Pointer atlasSlice = getCoronalAtlasSlice(coordinateOnAxis, atlasReferencePath);
    debugOut("Slicing Annotations");
    ImageType::Pointer annotationSlice = getCoronalAtlasSlice(coordinateOnAxis, atlasLabelsPath);

    debugOut("Saving to output paths");
    writeImage(atlasSlice, argv[1]);
    writeImage(annotationSlice, argv[2]);
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

void debugOut(const char * msg) {
    std::cout << "[Debug] " << msg << std::endl;
}

void writeImage(ImageType::Pointer im, const char * path) {
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer fileWriter = WriterType::New();
    fileWriter->SetFileName(path);
    fileWriter->SetInput(im);
    fileWriter->Update();
}