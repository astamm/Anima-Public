#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMExport.h>
#include <animaVanGelderenCompartment.h>

namespace anima
{

    class ANIMAMCM_EXPORT CylinderGPDPulsedGradientCompartment : public BaseCompartment
    {
    public:
        // Useful typedefs
        using Self = CylinderGPDPulsedGradientCompartment;
        using Superclass = BaseCompartment;
        using BasePointer = Superclass::Pointer;
        using Pointer = itk::SmartPointer<Self>;
        using ConstPointer = itk::SmartPointer<const Self>;
        using ModelOutputVectorType = Superclass::ModelOutputVectorType;
        using Vector3DType = Superclass::Vector3DType;
        using Matrix3DType = Superclass::Matrix3DType;

        // New macro
        itkNewMacro(Self);

        /** Run-time type information (and related methods) */
        itkTypeMacro(CylinderGPDPulsedGradientCompartment, BaseCompartment);

        DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE
        {
            return CylinderGPDPulsedGradient;
        }

        virtual double GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
        virtual ListType &GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
        virtual double GetLogDiffusionProfile(const Vector3DType &sample) ITK_OVERRIDE;

        virtual void SetParametersFromVector(const ListType &params) ITK_OVERRIDE;
        virtual ListType &GetParametersAsVector() ITK_OVERRIDE;

        virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
        virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

        // Set constraints
        void SetEstimateAxialDiffusivity(bool arg);
        void SetEstimateTissueRadius(bool arg);
        void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

        unsigned int GetCompartmentSize() ITK_OVERRIDE;
        unsigned int GetNumberOfParameters() ITK_OVERRIDE;
        ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

        void SetTissueRadius(double num) ITK_OVERRIDE;
        void SetAxialDiffusivity(double num) ITK_OVERRIDE;

        bool GetTensorCompatible() ITK_OVERRIDE { return false; }
        double GetApparentFractionalAnisotropy() ITK_OVERRIDE;

    protected:
        CylinderGPDPulsedGradientCompartment() : Superclass()
        {
            m_EstimateAxialDiffusivity = true;
            m_EstimateTissueRadius = true;
            m_ChangedConstraints = true;
            m_VanGelderenCompartment = anima::VanGelderenCompartment::New();
        }

        virtual ~CylinderGPDPulsedGradientCompartment() {}

    private:
        bool m_EstimateAxialDiffusivity, m_EstimateTissueRadius;
        bool m_ChangedConstraints;
        unsigned int m_NumberOfParameters;

        anima::VanGelderenCompartment::Pointer m_VanGelderenCompartment;
    };

} // end namespace anima
