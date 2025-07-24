#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>

#include <animaReadWriteFunctions.h>

int main(int ac, const char** av)
{
    TCLAP::CmdLine cmd("Outputs the 3D physical coordinates of a 3D image.\n"
                       "INRIA / IRISA - VisAGeS/Empenn Team",' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inputArg("i",
            "input",
            "input_filename",
            true,
            "",
            "Input image.",
            cmd);

    TCLAP::ValueArg<std::string> outputArg("o",
            "output",
            "output_filename",
            false,
            "",
            "Output image.",
            cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    using ImageType = itk::Image<double, 3>;
    using InputIteratorType = itk::ImageRegionConstIteratorWithIndex<ImageType>;

    using VectorImageType = itk::VectorImage<double, 3>;
    using OutputIteratorType = itk::ImageRegionIterator<VectorImageType>;

    ImageType::Pointer inputImage = anima::readImage<ImageType>(inputArg.getValue());
    VectorImageType::Pointer outputImage = VectorImageType::New();
    outputImage->Initialize();
    outputImage->SetRegions(inputImage->GetLargestPossibleRegion());
    outputImage->SetSpacing(inputImage->GetSpacing());
    outputImage->SetOrigin(inputImage->GetOrigin());
    outputImage->SetDirection(inputImage->GetDirection());
    outputImage->SetVectorLength(3);
    outputImage->Allocate();

    InputIteratorType inItr(inputImage, inputImage->GetLargestPossibleRegion());
    OutputIteratorType outItr(outputImage, outputImage->GetLargestPossibleRegion());

    ImageType::IndexType index;
    ImageType::PointType point;
    VectorImageType::PixelType pixel(3);

    while (!inItr.IsAtEnd())
    {
        index = inItr.GetIndex();
        inputImage->TransformIndexToPhysicalPoint(index, point);
        for (int i = 0; i < 3; ++i)
            pixel[i] = point[i];
        outItr.Set(pixel);
        ++inItr;
        ++outItr;
    }

    anima::writeImage<VectorImageType>(outputArg.getValue(), outputImage);

    return EXIT_SUCCESS;
}
