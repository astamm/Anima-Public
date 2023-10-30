#include <animaMultiCompartmentModelCreator.h>
#include <animaMCMFileWriter.h>
#include <animaMCMImage.h>
#include <animaReadWriteFunctions.h>
#include <animaVectorOperations.h>

#include <itkImage.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkVectorImage.h>

#include <tclap/CmdLine.h>

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    // Required arguments
    TCLAP::ValueArg<std::string> globalWeightsArg(
        "w", "global-weights", 
        "A 3D vector image specifying the weights of free water, isotropic restricted and anisotropic compartments.", 
        true, "", "global weights image", cmd
    );
    TCLAP::ValueArg<std::string> irwDiffusivityArg(
        "d", "irw-diffusivity", 
        "A 3D scalar image specifying the diffusivity of the isotropic restricted compartment.", 
        true, "", "irw diffusivity image", cmd
    );
    TCLAP::MultiArg<std::string> anisoOrientationArg(
        "i", "aniso-orientation", 
        "A 3D vector image specifying the orientation of the anisotropic compartment (accepted mulitple times).", 
        true, "aniso orientation image", cmd
    );
    TCLAP::ValueArg<std::string> anisoWeightsArg(
        "a", "aniso-weights", 
        "A 3D vector image specifying the relative anisotropic weights.", 
        true, "", "aniso weights image", cmd
    );
    TCLAP::ValueArg<std::string> outputMCMArg(
        "o", "output-mcm", 
        "A vector image storing the output MCM image.", 
        true, "", "output mcm image", cmd
    );

    // Optional arguments
    TCLAP::ValueArg<std::string> irwRadiusArg(
        "r", "irw-radius", 
        "A 3D scalar image specifying the radius for the isotropic restricted compartment.", 
        false, "", "irw radius image", cmd
    );

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    using ScalarImageType = itk::Image<double,3>;
    using ScalarIteratorType = itk::ImageRegionConstIterator<ScalarImageType>;
    using VectorImageType = itk::VectorImage<double,3>;
    using VectorIteratorType = itk::ImageRegionConstIterator<VectorImageType>;
    using OutputImageType = anima::MCMImage<double,3>;
    using OutputIteratorType = itk::ImageRegionIterator<OutputImageType>;
    using MCMFileWriterType = anima::MCMFileWriter<double,3>;
    using MCMCreatorType = anima::MultiCompartmentModelCreator;
    using MCMPointer = MCMCreatorType::MCMPointer;

    // Load images
    VectorImageType::Pointer globalWeightsImage = anima::readImage<VectorImageType>(globalWeightsArg.getValue());
    ScalarImageType::Pointer irwDiffusivityImage = anima::readImage<ScalarImageType>(irwDiffusivityArg.getValue());
    std::vector <std::string> anisoOrientationFiles = anisoOrientationArg.getValue();
    unsigned int numberOfAnisotropicCompartments = anisoOrientationFiles.size();
    std::vector <VectorImageType::Pointer> anisoOrientationImages(numberOfAnisotropicCompartments);
    for (unsigned int i = 0;i < numberOfAnisotropicCompartments;++i)
        anisoOrientationImages[i] = anima::readImage<VectorImageType>(anisoOrientationFiles[i]);
    VectorImageType::Pointer anisoWeightsImage = anima::readImage<VectorImageType>(anisoWeightsArg.getValue());

    ScalarImageType::Pointer irwRadiusImage;

    // Generate output MCM image
    OutputImageType::Pointer outputImage = OutputImageType::New();
    outputImage->SetRegions(irwDiffusivityImage->GetLargestPossibleRegion());
    outputImage->CopyInformation(irwDiffusivityImage);

    // Create fake MCM output to get its length
    MCMCreatorType *tmpMCMCreator = new MCMCreatorType;
    tmpMCMCreator->SetModelWithFreeWaterComponent(true);
    tmpMCMCreator->SetModelWithStationaryWaterComponent(false);
    if (irwRadiusArg.getValue() == "")
    {
        tmpMCMCreator->SetModelWithRestrictedWaterComponent(true);
        tmpMCMCreator->SetModelWithStaniszComponent(false);
    }
    else
    {
        irwRadiusImage = anima::readImage<ScalarImageType>(irwRadiusArg.getValue());
        tmpMCMCreator->SetModelWithRestrictedWaterComponent(false);
        tmpMCMCreator->SetModelWithStaniszComponent(true);
    }
    tmpMCMCreator->SetCompartmentType(anima::Zeppelin);
    tmpMCMCreator->SetNumberOfCompartments(numberOfAnisotropicCompartments);

    MCMPointer tmpMCM = tmpMCMCreator->GetNewMultiCompartmentModel();

    outputImage->SetVectorLength(tmpMCM->GetSize());
    outputImage->SetDescriptionModel(tmpMCM);
    outputImage->Allocate();

    delete tmpMCMCreator;

    // Define global useful variables
    std::vector<double> compartmentWeights(numberOfAnisotropicCompartments + 2);
    VectorImageType::PixelType globalWeights(3), anisoWeights(numberOfAnisotropicCompartments), anisoOrientation(3), anisoAngles(3);
    OutputImageType::PixelType outputValue(tmpMCM->GetSize());
    outputValue.Fill(0.0);
    outputImage->FillBuffer(outputValue);

    // Define iterators
    VectorIteratorType globalWeightsItr(globalWeightsImage, irwDiffusivityImage->GetLargestPossibleRegion());
    ScalarIteratorType irwDiffusivityItr(irwDiffusivityImage, irwDiffusivityImage->GetLargestPossibleRegion());
    ScalarIteratorType irwRadiusItr;
    if (irwRadiusImage)
        irwRadiusItr = ScalarIteratorType(irwRadiusImage, irwDiffusivityImage->GetLargestPossibleRegion());
    std::vector<VectorIteratorType> anisoOrientationItrs(numberOfAnisotropicCompartments);
    for (unsigned int i = 0;i < numberOfAnisotropicCompartments;++i)
        anisoOrientationItrs[i] = VectorIteratorType(anisoOrientationImages[i], irwDiffusivityImage->GetLargestPossibleRegion());
    VectorIteratorType anisoWeightsItr(anisoWeightsImage, irwDiffusivityImage->GetLargestPossibleRegion());
    OutputIteratorType outputItr(outputImage, irwDiffusivityImage->GetLargestPossibleRegion());

	while (!outputItr.IsAtEnd())
	{
        if (irwDiffusivityItr.Get() == 0)
        {
            ++globalWeightsItr;
            ++irwDiffusivityItr;
            if (irwRadiusImage)
                ++irwRadiusItr;
            for (unsigned int i = 0;i < numberOfAnisotropicCompartments;++i)
                ++anisoOrientationItrs[i];
            ++anisoWeightsItr;
            ++outputItr;
            continue;
        }

        globalWeights = globalWeightsItr.Get();
        anisoWeights = anisoWeightsItr.Get();

        compartmentWeights[0] = globalWeights[0];
        compartmentWeights[1] = globalWeights[1];
        for (unsigned int i = 0;i < numberOfAnisotropicCompartments;++i)
        {
            compartmentWeights[i + 2] = globalWeights[2] * anisoWeights[i];
            anisoOrientation = anisoOrientationItrs[i].Get();
            anima::TransformCartesianToSphericalCoordinates(anisoOrientation, anisoAngles);
            tmpMCM->GetCompartment(i + 2)->SetOrientationTheta(anisoAngles[0]);
            tmpMCM->GetCompartment(i + 2)->SetOrientationPhi(anisoAngles[1]);
            tmpMCM->GetCompartment(i + 2)->SetAxialDiffusivity(1.71e-3);
            tmpMCM->GetCompartment(i + 2)->SetRadialDiffusivity1(1.71e-4);
        }
            
        tmpMCM->SetCompartmentWeights(compartmentWeights);

        tmpMCM->GetCompartment(1)->SetAxialDiffusivity(irwDiffusivityItr.Get());
        if (irwRadiusImage)
            tmpMCM->GetCompartment(1)->SetTissueRadius(irwRadiusItr.Get());

        outputValue = tmpMCM->GetModelVector();
        outputItr.Set(outputValue);

         ++globalWeightsItr;
        ++irwDiffusivityItr;
        if (irwRadiusImage)
            ++irwRadiusItr;
        for (unsigned int i = 0;i < numberOfAnisotropicCompartments;++i)
            ++anisoOrientationItrs[i];
        ++anisoWeightsItr;
        ++outputItr;
	}

    MCMFileWriterType writer;
    writer.SetInputImage(outputImage);
    writer.SetFileName(outputMCMArg.getValue());
    writer.Update();

    return EXIT_SUCCESS;
}
