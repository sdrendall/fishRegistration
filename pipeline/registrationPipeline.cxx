// Registration pipeline for registering FISH images to the 
//  Allen Brain Reference Atlas
//
// Heavily based off of the ImageRegistration6.cxx example
//  included with the itk 

// For rigid registration
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkCenteredTransformInitializer.h"

// File readers and writers
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// For data management
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSubtractImageFilter.h"

//Command observer -- monitors the evolution of registration
#include "itkCommand.h"

class CommandIterationUpdate : public itk::Comand {
public:
    typedef CommandIterationUpdate Self;
    typedef itk::Command Superclass;
    typedef itk::SmartPointer<Self>
}