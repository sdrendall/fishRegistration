#include <stdio.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"

int main(int argc, char *argv[]){
    // Check arguments
    if (argc < 4) {
            std::cerr << "Missing Parameters " << std::endl;
            std::cerr << "Usage: " << argv[0];
            std::cerr << " 3DImagePath 2DBaseOutputPath ";
            return EXIT_FAILURE;
        }

    // Declare 3D and 2D image types
    typedef unsigned char PixelType3D;
    typedef unsigned char PixelType2D;
    
    typedef itk::Image<PixelType3D, 3> ImageType3D;
    typedef itk::Image<PixelType2D, 2> ImageType2D;
    
    // Declare Readers and Writers for 3D images and 2D images respectively
    typedef itk::ImageFileReader<ImageType3D> ReaderType3D;
    typedef itk::ImageFileWriter<ImageType2D> WriterType2D;
    typedef itk::ExtractImageFilter<ImageType3D, ImageType2D> FilterType3Dto2D;
    
    // Create a reader, to read the input image
    ReaderType3D::Pointer reader = ReaderType3D::New();

    // The input filename is the first command line argument
    const char * filename = argv[1];
    reader->SetFileName(filename);

    // Update, to load the image so that size information can be ascertained
    reader->Update();

    // Get the whole image
    ImageType3D::RegionType inputRegion = reader->GetOutput()->GetLargestPossibleRegion();

    // Determine Slice Size
    ImageType3D::SizeType inputSize = inputRegion.GetSize();
    ImageType3D::SizeType sliceSize = inputSize;
    sliceSize[2] = 0; // Slices have no z dimension

    // Initialize a slice region
    ImageType3D::RegionType sliceRegion;
    sliceRegion.SetSize(sliceSize);

    // Get Start Point
    ImageType3D::IndexType sliceStart = inputRegion.GetIndex();

    // Reference an extraction filter and connect it to the reader
    FilterType3Dto2D::Pointer extFilter = FilterType3Dto2D::New();
    extFilter->SetInput(reader->GetOutput());
    extFilter->SetDirectionCollapseToSubmatrix();

    // Create a FileWriter to save the images -- Connect it to the extraction filter
    WriterType2D::Pointer writer = WriterType2D::New();
    writer->SetInput(extFilter->GetOutput());

    // Iterate over each slice along the z dimension
    const unsigned int numSlices = inputSize[2];
    printf("Number of Slices: %d", numSlices);
    for (unsigned int i = 0; i < numSlices; i++) {
        // Move the slice region to the next slice
        sliceStart[2] = i;
        sliceRegion.SetIndex(sliceStart);
        // Pass the extraction filter the new slice region
        extFilter->SetExtractionRegion(sliceRegion);
        // Generate Filename
        char outputPath[50];
        sprintf(outputPath, "%s%d.jpg", argv[2], i);
        // Set new output path
        writer->SetFileName(outputPath);
        // Update writer to pull everything through
        writer->Update();
    }
    return EXIT_SUCCESS;
}