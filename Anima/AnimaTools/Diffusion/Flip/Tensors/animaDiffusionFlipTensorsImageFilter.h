#pragma once

#include <animaMaskedImageToImageFilter.h>

#include <itkVectorImage.h>

namespace anima
{

template <class TPixelType, unsigned int TImageDimension = 3>
class DiffusionFlipTensorsImageFilter :
    public anima::MaskedImageToImageFilter< itk::VectorImage <TPixelType,TImageDimension>,
                                    itk::VectorImage <TPixelType,TImageDimension> >
{
public:
    /** Standard class typedefs. */
    using Self = DiffusionFlipTensorsImageFilter;
    
    using InputImageType = itk::VectorImage<TPixelType,TImageDimension>;
    using InputPixelType = typename InputImageType::PixelType;
    
    using OutputImageType = itk::VectorImage<TPixelType,TImageDimension>;
    using OutputPixelType = typename OutputImageType::PixelType;
    
    using Superclass = anima::MaskedImageToImageFilter<InputImageType, OutputImageType>;
    
    using Pointer = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer<const Self>;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    
    /** Run-time type information (and related methods) */
    itkTypeMacro(DiffusionFlipTensorsImageFilter, MaskedImageToImageFilter);
    
    /** Superclass typedefs. */
    using InputImageRegionType = typename Superclass::InputImageRegionType;
    using OutputImageRegionType = typename Superclass::OutputImageRegionType;
    using MaskImageType = typename Superclass::MaskImageType;
    
    using InputPointerType = typename InputImageType::Pointer;
    using OutputPointerType = typename OutputImageType::Pointer;
    
    itkSetMacro(FlippedAxis, std::string);

protected:
    DiffusionFlipTensorsImageFilter() : Superclass()
    {
        m_FlippedAxis = "";
    }
    
    virtual ~DiffusionFlipTensorsImageFilter() {}

    void GenerateOutputInformation() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(DiffusionFlipTensorsImageFilter);
    
    std::string m_FlippedAxis;
    
    static const unsigned int m_NumberOfComponents = 6;
};

} // end of namespace anima

#include "animaDiffusionFlipTensorsImageFilter.hxx"
