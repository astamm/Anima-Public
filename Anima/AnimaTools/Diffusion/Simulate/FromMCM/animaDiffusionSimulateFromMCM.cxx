#include <animaGradientFileReader.h>
#include <animaMCMConstants.h>
#include <animaMCMFileReader.h>
#include <animaMCMImage.h>
#include <animaMultiCompartmentModel.h>
#include <animaReadWriteFunctions.h>

#include <boost/lexical_cast.hpp>

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

#include <cmath>
#include <fstream>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputmcm","mcm image",true,"","mcm image",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","Result DWI volume",true,"","result DWI volume",cmd);
    TCLAP::ValueArg<double> smallDeltaArg("", "small-delta", "Diffusion small delta (in seconds)", false, anima::DiffusionSmallDelta, "small delta", cmd);
    TCLAP::ValueArg<double> bigDeltaArg("", "big-delta", "Diffusion big delta (in seconds)", false, anima::DiffusionBigDelta, "big delta", cmd);

    TCLAP::ValueArg<std::string> s0ValueArg("s","s0","S0 of DWI (constant value or image)",false,"500","S0 of DWI",cmd);
    TCLAP::ValueArg<std::string> gradsArg("g","grad","Input gradients",true,"","Input gradients",cmd);
    TCLAP::ValueArg<std::string> bvalArg("b","bval","Input b-values",true,"","Input b-values",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    using ScalarType = double;
    using BValueListType = std::vector<ScalarType>;

    using GradientType = vnl_vector_fixed<double,3>;
    using GradientListType = std::vector<GradientType >;

    using GFReaderType = anima::GradientFileReader < GradientType, double >;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(gradsArg.getValue());
    gfReader.SetBValueBaseString(bvalArg.getValue());
    gfReader.SetSmallDelta(smallDeltaArg.getValue());
    gfReader.SetBigDelta(bigDeltaArg.getValue());
    gfReader.SetGradientIndependentNormalization(false);
    gfReader.Update();

    GradientListType gradients = gfReader.GetGradients();

    GFReaderType::BValueVectorType gradientStrengths = gfReader.GetGradientStrengths();

    int numberOfDirections = gradientStrengths.size();

    using TInputImageType = anima::MCMImage <ScalarType, 3>;
    using MCMPointer = TInputImageType::MCMPointer;
    using TOutputImageType = itk::VectorImage <ScalarType, 3>;

    /* Image typedef support */
    using InputImagePointer = TInputImageType::Pointer;
    using OutputImagePointer = TOutputImageType::Pointer;

    std::cout << "Loading MCM image " << inArg.getValue() << std::endl;

    anima::MCMFileReader <ScalarType,3> mcmReader;
    mcmReader.SetFileName(inArg.getValue());
    mcmReader.Update();

    InputImagePointer mcmImage = mcmReader.GetModelVectorImage();

    using InputIteratorType = itk::ImageRegionIterator <TInputImageType>;
    using OutputIteratorType = itk::ImageRegionIterator <TOutputImageType>;

    OutputImagePointer outImage = TOutputImageType::New();
    outImage->Initialize();
    outImage->SetRegions(mcmImage->GetLargestPossibleRegion());
    outImage->SetOrigin(mcmImage->GetOrigin());
    outImage->SetDirection(mcmImage->GetDirection());
    outImage->SetSpacing(mcmImage->GetSpacing());
    outImage->SetVectorLength(numberOfDirections);
    outImage->Allocate();

    using InputPixelType = TInputImageType::PixelType;
    InputPixelType mcmValue;

    using OutputPixelType = TOutputImageType::PixelType;

    InputIteratorType mcmIt (mcmImage, mcmImage->GetLargestPossibleRegion());
    OutputIteratorType outIt (outImage, outImage->GetLargestPossibleRegion());

    OutputPixelType voxelOutputValue(numberOfDirections);

    MCMPointer referenceModel = mcmImage->GetDescriptionModel();

    using S0ImageType = itk::Image <double, 3>;
    S0ImageType::Pointer s0Image;
    try
    {
        s0Image = anima::readImage <S0ImageType> (s0ValueArg.getValue());
    }
    catch(itk::ExceptionObject &)
    {
        s0Image = 0;
    }

    using S0ImageIteratorType = itk::ImageRegionConstIterator <S0ImageType>;
    S0ImageIteratorType s0Itr;

    if (s0Image)
        s0Itr = S0ImageIteratorType(s0Image, s0Image->GetLargestPossibleRegion());

    while(!outIt.IsAtEnd())
    {
        mcmValue = mcmIt.Get();
        referenceModel->SetModelVector(mcmValue);

        double b0Value = 1;
        if (s0Image)
            b0Value = s0Itr.Get();
        else
            b0Value = boost::lexical_cast<double> (s0ValueArg.getValue());

        for (int dir = 0;dir < numberOfDirections;++dir)
            voxelOutputValue[dir] = b0Value * referenceModel->GetPredictedSignal(smallDeltaArg.getValue(),bigDeltaArg.getValue(),
                                                                                 gradientStrengths[dir],gradients[dir]);

        outIt.Set(voxelOutputValue);
        ++outIt;
        ++mcmIt;
        if (s0Image)
            ++s0Itr;
    }

    std::cout << "Writing DWI Image : " << resArg.getValue() << std::endl;
    anima::writeImage <TOutputImageType> (resArg.getValue(),outImage);

    return 0;
}
