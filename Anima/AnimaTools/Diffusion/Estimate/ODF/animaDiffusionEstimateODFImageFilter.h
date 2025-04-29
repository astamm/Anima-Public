#pragma once

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

#include <iostream>
#include <vector>

namespace anima
{

    template <typename TInputPixelType, typename TOutputPixelType>
    class DiffusionEstimateODFImageFilter : public itk::ImageToImageFilter<itk::Image<TInputPixelType, 3>, itk::VectorImage<TOutputPixelType, 3>>
    {
    public:
        /** Standard class typedefs. */
        using Self = DiffusionEstimateODFImageFilter;
        using TInputImage = itk::Image<TInputPixelType, 3>;
        using Image4DType = itk::Image<TInputPixelType, 4>;
        using OutputScalarImageType = itk::Image<TOutputPixelType, 3>;
        using TOutputImage = itk::VectorImage<TOutputPixelType, 3>;
        using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
        using Pointer = itk::SmartPointer<Self>;
        using ConstPointer = itk::SmartPointer<const Self>;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods) */
        itkTypeMacro(DiffusionEstimateODFImageFilter, ImageToImageFilter);

        using InputImagePointer = typename TInputImage::Pointer;
        using OutputImagePointer = typename TOutputImage::Pointer;
        using OutputScalarImagePointer = typename OutputScalarImageType::Pointer;

        /** Superclass typedefs. */
        using OutputImageRegionType = typename Superclass::OutputImageRegionType;

        void AddGradientDirection(unsigned int i, std::vector<double> &grad);
        void SetBValuesList(std::vector<double> bValuesList) { m_BValuesList = bValuesList; }
        itkSetMacro(BValueShellSelected, int);

        itkSetMacro(Lambda, double);
        itkSetMacro(LOrder, unsigned int);

        itkSetMacro(SharpnessRatio, double);
        itkSetMacro(Sharpen, bool);

        itkSetMacro(Normalize, bool);
        itkSetMacro(FileNameSphereTesselation, std::string);

        itkSetMacro(UseAganjEstimation, bool);
        itkSetMacro(DeltaAganjRegularization, double);

        itkGetMacro(EstimatedB0Image, OutputScalarImageType *);
        itkGetMacro(EstimatedVarianceImage, OutputScalarImageType *);

        void SetReferenceB0Image(TInputImage *refB0)
        {
            m_ReferenceB0Image = refB0;
        }

    protected:
        DiffusionEstimateODFImageFilter()
        {
            m_GradientDirections.clear();
            m_PVector.clear();
            m_ReferenceB0Image = nullptr;

            m_BValueShellSelected = -1;
            m_BValueShellTolerance = 20;

            m_Lambda = 0.006;
            m_DeltaAganjRegularization = 0.001;

            m_LOrder = 4;

            m_Sharpen = false;
            m_SharpnessRatio = 0.255;

            m_Normalize = false;
            m_SphereSHSampling.clear();

            m_UseAganjEstimation = false;
        }

        virtual ~DiffusionEstimateODFImageFilter() {}

        void GenerateOutputInformation() ITK_OVERRIDE;
        void BeforeThreadedGenerateData() ITK_OVERRIDE;
        void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    private:
        ITK_DISALLOW_COPY_AND_ASSIGN(DiffusionEstimateODFImageFilter);

        bool isZero(std::vector<double> &testVal)
        {
            bool resVal = true;
            for (unsigned int i = 0; i < testVal.size(); ++i)
            {
                if (testVal[i] != 0)
                {
                    resVal = false;
                    break;
                }
            }

            return resVal;
        }

        std::vector<std::vector<double>> m_GradientDirections;
        std::vector<double> m_BValuesList;
        InputImagePointer m_ReferenceB0Image;

        OutputScalarImagePointer m_EstimatedVarianceImage;
        OutputScalarImagePointer m_EstimatedB0Image;

        int m_BValueShellSelected;
        double m_BValueShellTolerance;
        std::vector<unsigned int> m_SelectedDWIIndexes;

        vnl_matrix<double> m_TMatrix; // evaluation matrix computed once and for all before threaded generate data
        vnl_matrix<double> m_BMatrix;
        std::vector<double> m_DeconvolutionVector;
        std::vector<double> m_PVector;

        std::vector<unsigned int> m_B0Indexes, m_GradientIndexes;

        bool m_Normalize;
        std::string m_FileNameSphereTesselation;
        std::vector<std::vector<double>> m_SphereSHSampling;

        double m_Lambda;
        double m_SharpnessRatio; // See Descoteaux et al. TMI 2009, article plus appendix
        bool m_Sharpen;
        bool m_UseAganjEstimation;
        double m_DeltaAganjRegularization;
        unsigned int m_LOrder;
    };

} // end of namespace anima

#include "animaDiffusionEstimateODFImageFilter.hxx"
