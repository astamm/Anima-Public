#pragma once

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

#include <iostream>
#include <vector>

namespace anima
{

template <typename TInputPixelType>
class DiffusionExtractFromODFImageFilter :
public itk::ImageToImageFilter< itk::VectorImage<TInputPixelType, 3> , itk::Image <TInputPixelType, 3> >
{
public:
    /** Standard class typedefs. */
    using Self = DiffusionExtractFromODFImageFilter;
    using TInputImage = itk::VectorImage <TInputPixelType, 3>;
    using TOutputImage = itk::Image <TInputPixelType, 3>;
    using Superclass = itk::ImageToImageFilter< TInputImage, TOutputImage >;
    using Pointer = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer<const Self>;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(DiffusionExtractFromODFImageFilter, ImageToImageFilter);

    using InputImagePointer = typename TInputImage::Pointer;
    using InputImagePixel = typename TInputImage::PixelType;
    using OutputImagePointer = typename TOutputImage::Pointer;

    /** Superclass typedefs. */
    using OutputImageRegionType = typename Superclass::OutputImageRegionType;

protected:
    DiffusionExtractFromODFImageFilter()
    {
    }

    virtual ~DiffusionExtractFromODFImageFilter() {}

    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    bool isZero(InputImagePixel &testVal)
    {
        bool resVal = true;
        for (unsigned int i = 0;i < testVal.GetSize();++i)
        {
            if (testVal[i] != 0)
            {
                resVal = false;
                break;
            }
        }

        return resVal;
    }

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(DiffusionExtractFromODFImageFilter);
};
	
} // end of namespace anima

#include "animaDiffusionExtractFromODFImageFilter.hxx"
