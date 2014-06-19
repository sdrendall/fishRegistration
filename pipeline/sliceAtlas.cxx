#include <stdio.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"

const char atlasPath[80] = "/home/sam/Pictures/allenReferenceAtlas_mouseCoronal/atlasVolume/atlasVolume.mhd";
int axisToCollapse;

int main(int argc, char *argv[]){
    // Check arguments
    if (argc < 4) {
        axisToCollapse = atoi(argv[3]);
    } else if (argc < 3) {
            std::cerr << "Missing Parameters " << std::endl;
            std::cerr << "Usage: " << argv[0];
            std::cerr << " outputPath coordinateOnAxis";
            std::cerr << " [axisToCollapse] ";
            return EXIT_FAILURE;
    } 

    // Declare 3D and 2D image types
    typedef unsigned char AtlasPixelType;
    typedef unsigned char SlicePixelType;
    
    typedef itk::Image<AtlasPixelType, 3> AtlasImageType;
    typedef itk::Image<SlicePixelType, 2> SliceImageType;
    
    // Declare Readers and Writers for 3D images and 2D images respectively
    typedef itk::ImageFileReader<AtlasImageType> AtlasReaderType;
    typedef itk::ImageFileWriter<SliceImageType> SliceWriterType;
    typedef itk::ExtractImageFilter<AtlasImageType, SliceImageType> SliceAtlasFilterType;
    
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
    extFilter->SetDirectionCollapseToSubmatrix();

    // Create a FileWriter to save the images -- Connect it to the extraction filter
    SliceWriterType::Pointer writer = SliceWriterType::New();
    writer->SetInput(extFilter->GetOutput());

    // Find the corresponding atlas index to the point specified by the user
    typedef itk::Point <double, AtlasImageType::ImageDimension> PointType;
    PointType point; // Let's hope this initializes with 0s...
    point[axisToCollapse] = atoi(argv[2]);
    bool isInside = atlas->TransformPhysicalPointToIndex(point, sliceStartIndex);

    // Extract the slice, if it is in the image
    if (isInside) {
        sliceRegion.SetIndex(sliceStartIndex);
        extFilter->SetExtractionRegion(sliceRegion);
        // Generate Filename
        char outputPath[50];
        sprintf(outputPath, "%s.mhd", argv[1]);
        // Set new output path
        writer->SetFileName(outputPath);
        // Update writer to pull everything through
        writer->Update();
    } else {
        std::cerr << "Plane outside of the atlas boundries!" << std::endl;
        std::cerr << "Please specify a point within the atlas..." << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}