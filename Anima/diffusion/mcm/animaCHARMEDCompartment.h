#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMExport.h>
#include <animaStaniszCompartment.h>

#include <tuple>
#include <map>

namespace anima
{

    class ANIMAMCM_EXPORT CHARMEDCompartment : public BaseCompartment
    {
    public:
        // Useful typedefs
        using Self = CHARMEDCompartment;
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
        itkTypeMacro(CHARMEDCompartment, BaseCompartment);

        DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE
        {
            return CHARMED;
        }

        virtual double GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
        virtual ListType &GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
        virtual double GetLogDiffusionProfile(const Vector3DType &sample) ITK_OVERRIDE;

        virtual void SetParametersFromVector(const ListType &params) ITK_OVERRIDE;
        virtual ListType &GetParametersAsVector() ITK_OVERRIDE;

        virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
        virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

        // Set constraints
        void SetEstimateDiffusivities(bool arg);
        void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

        unsigned int GetCompartmentSize() ITK_OVERRIDE;
        unsigned int GetNumberOfParameters() ITK_OVERRIDE;
        ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

        void SetTissueRadius(double num) ITK_OVERRIDE;
        void SetAxialDiffusivity(double num) ITK_OVERRIDE;

        bool GetTensorCompatible() ITK_OVERRIDE { return false; }
        double GetApparentFractionalAnisotropy() ITK_OVERRIDE;

    protected:
        CHARMEDCompartment() : Superclass()
        {
            m_EstimateDiffusivities = true;
            m_ChangedConstraints = true;
            m_StaniszCompartment = anima::StaniszCompartment::New();
            m_StaniszCompartment->SetEstimateTissueRadius(true);
        }

        virtual ~CHARMEDCompartment() {}

        void PrepareData(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient);

    private:
        bool m_EstimateDiffusivities;
        bool m_ChangedConstraints;
        unsigned int m_NumberOfParameters;
        
        // Parameters for data perparation
        anima::StaniszCompartment::Pointer m_StaniszCompartment;
        double m_AxialSignal, m_RadialHinderedSignal, m_RadialRestrictedSignal;
        double m_OrthogonalGradientStrength;
        Vector3DType m_OrthogonalGradient;
        double m_BValue, m_InnerProd;
    };

} // end namespace anima
