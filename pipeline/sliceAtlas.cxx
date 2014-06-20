#include <stdio.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"

const char * atlasPath = "/home/sam/Pictures/allenReferenceAtlas_mouseCoronal/atlasVolume/atlasVolume.mhd";

void debugOut(const char * msg) {
    std::cout << "[Debug] " << msg << std::endl;
}

int main(int argc, char *argv[]){
    debugOut("Checking input arguments...");
    int axisToCollapse;
    // Check arguments
    if (argc < 3) {
            std::cerr << "Missing Parameters " << std::endl;
            std::cerr << "Usage: " << argv[0];
            std::cerr << " outputPath coordinateOnAxis";
            std::cerr << " [axisToCollapse] ";
            return EXIT_FAILURE;
    } else if (argc >= 4) {
        axisToCollapse = atoi(argv[3]);
    } else {
        axisToCollapse = 1;
    }

    debugOut("Declaring Pixel, Image, Reader and Filter types...");
    // Declare 3D and 2D image types
    typedef unsigned char AtlasPixelType;
    typedef unsigned char SlicePixelType;
    
    typedef itk::Image<AtlasPixelType, 3> AtlasImageType;
    typedef itk::Image<SlicePixelType, 2> SliceImageType;
    
    // Declare Readers and Writers for 3D images and 2D images respectively
    typedef itk::ImageFileReader<AtlasImageType> AtlasReaderType;
    typedef itk::ImageFileWriter<SliceImageType> SliceWriterType;
    typedef itk::ExtractImageFilter<AtlasImageType, SliceImageType> SliceAtlasFilterType;
    
    debugOut("Instantiating Reader...");
    // Create a reader, to read the input image
    AtlasReaderType::Pointer reader = AtlasReaderType::New();
    reader->SetFileName(atlasPath);

    debugOut("Updating Reader...");
    // Update, to load the image so that size information can be ascertained
    reader->Update();

    debugOut("Referencing atlas...");
    // Get the whole atlas region
    AtlasImageType::Pointer atlas = reader->GetOutput();
    AtlasImageType::RegionType entireAtlasRegion = atlas->GetLargestPossibleRegion();

    debugOut("Computing slice size...");
    // Determine Slice Size
    AtlasImageType::SizeType inputSize = entireAtlasRegion.GetSize();
    AtlasImageType::SizeType sliceSize = inputSize;
    sliceSize[axisToCollapse] = 0;  // 0 tells the ExtractionFilter to return an Image without that dimension

    // Initialize a slice region
    AtlasImageType::RegionType sliceRegion;
    sliceRegion.SetSize(sliceSize);

    debugOut("Getting sliceStartIndex...");
    // Get Start Point
    AtlasImageType::IndexType sliceStartIndex = entireAtlasRegion.GetIndex();

    debugOut("Creating extraction filter...");
    // Reference an extraction filter and connect it to the reader
    SliceAtlasFilterType::Pointer extFilter = SliceAtlasFilterType::New();
    extFilter->SetInput(reader->GetOutput());
    extFilter->SetDirectionCollapseToIdentity();

    debugOut("Creating FileWriter...");
    // Create a FileWriter to save the images -- Connect it to the extraction filter
    SliceWriterType::Pointer writer = SliceWriterType::New();
    writer->SetInput(extFilter->GetOutput());

    debugOut("Computing array index of physical coordinate...");
    // Find the corresponding atlas index to the point specified by the user
    typedef itk::Point <double, AtlasImageType::ImageDimension> PointType;
    PointType point;
    point.Fill(0.0);
    point[axisToCollapse] = atoi(argv[2]);
    bool isInside = atlas->TransformPhysicalPointToIndex(point, sliceStartIndex);

    debugOut("Extracting Slice...");
    // Extract the slice, if it is in the image
    if (isInside) {
        debugOut("Configuring Extraction Filter...");
        sliceRegion.SetIndex(sliceStartIndex);
        extFilter->SetExtractionRegion(sliceRegion);
        // Generate Filename
        char outputPath[50];
        sprintf(outputPath, "%s.jpg", argv[1]);
        // Set new output path
        writer->SetFileName(outputPath);
        // Update writer to pull everything through
        debugOut("Updating Writer...");
        writer->Update();
    } else {
        std::cerr << "Plane outside of the atlas boundries!" << std::endl;
        std::cerr << "Please specify a point within the atlas..." << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}