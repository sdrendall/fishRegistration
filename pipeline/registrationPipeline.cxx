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
class CommandIterationUpdate : public itk::Command
{
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

int main(int argc, char *argv[]) {
    if(argc < 4) {
        std::cerr << "Missing Parameters "<< std::endl;
        std::cerr << "Usage: " << argv[0];
        std::cerr << " fixedImageFile movingImageFile ";
        std::cerr << " outputImagefile  [differenceBeforeRegistration] ";
        std::cerr << " [differenceAfterRegistration] "<< std::endl;
        return EXIT_FAILURE;
    }
}