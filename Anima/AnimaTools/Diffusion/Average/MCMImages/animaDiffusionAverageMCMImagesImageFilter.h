#pragma once

#include <animaMCMImage.h>
#include <animaMCMWeightedAverager.h>
#include <animaMultiCompartmentModel.h>

#include <itkImageToImageFilter.h>

namespace anima
{

template <class TPixelType>
class DiffusionAverageMCMImagesImageFilter :
public itk::ImageToImageFilter< anima::MCMImage <TPixelType, 3>, anima::MCMImage<TPixelType, 3> >
{
public:
    /** Standard class type def */

    using Self = DiffusionAverageMCMImagesImageFilter;
    using InputImageType = anima::MCMImage <TPixelType, 3>;
    using OutputImageType = anima::MCMImage <TPixelType, 3>;
    using Superclass = itk::ImageToImageFilter <InputImageType, OutputImageType >;
    using Pointer = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer<const Self>;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(DiffusionAverageMCMImagesImageFilter, ImageToImageFilter);

    using InputImagePointer = typename InputImageType::ConstPointer;
    using OutputImagePointer = typename OutputImageType::Pointer;

    using InputRegionType = typename InputImageType::RegionType;
    using InputIndexType = typename InputImageType::IndexType;
    using PixelType = typename InputImageType::PixelType;

    using MaskImageType = itk::Image <unsigned char, 3>;
    using MaskImagePointer = typename MaskImageType::Pointer;

    // Multi-compartment models typedefs
    using MCModelType = anima::MultiCompartmentModel;
    using MCModelPointer = MCModelType::Pointer;

    using MCMAveragerType = anima::MCMWeightedAverager;
    using MCMAveragerPointer = typename MCMAveragerType::Pointer;

    void SetReferenceOutputModel(MCModelPointer &model);
    MCModelType *GetReferenceOutputModel() {return m_ReferenceOutputModel;}

    void AddMaskImage(MaskImageType *maskImage) {m_MaskImages.push_back(maskImage);}
    
    itkSetMacro(DDIAveragingMethod, int);

protected:
    DiffusionAverageMCMImagesImageFilter() {m_DDIAveragingMethod = 3;}
    virtual ~DiffusionAverageMCMImagesImageFilter() {}

    bool isZero(const itk::VariableLengthVector <double> &value) const
    {
        for (unsigned int i = 0;i < value.GetNumberOfElements();++i)
        {
            if (value[i] != 0)
                return false;
        }

        return true;
    }

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const InputRegionType &region) ITK_OVERRIDE;
    virtual MCMAveragerPointer CreateAverager();

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(DiffusionAverageMCMImagesImageFilter);

    std::vector <MCModelPointer> m_ReferenceInputModels;
    std::vector <MaskImagePointer> m_MaskImages;
    MCModelPointer m_ReferenceOutputModel;
    int m_DDIAveragingMethod;
};

} // end namespace anima

#include "animaDiffusionAverageMCMImagesImageFilter.hxx"
