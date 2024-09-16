#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMExport.h>
#include <animaNeumanCylinderCompartment.h>
#include <animaStaniszCylinderCompartment.h>
#include <animaVanGelderenCylinderCompartment.h>
#include <animaZeppelinCompartment.h>

#include <map>
#include <tuple>


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
        void SetEstimateTissueRadius(bool arg);
        void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

        unsigned int GetCompartmentSize() ITK_OVERRIDE;
        unsigned int GetNumberOfParameters() ITK_OVERRIDE;
        ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

        void SetOrientationTheta(double num) ITK_OVERRIDE;
        void SetOrientationPhi(double num) ITK_OVERRIDE;
        void SetExtraAxonalFraction(double num) ITK_OVERRIDE;
        void SetAxialDiffusivity(double num) ITK_OVERRIDE;
        // void SetRadialDiffusivity1(double num) ITK_OVERRIDE;
        void SetTissueRadius(double num) ITK_OVERRIDE;

        bool GetTensorCompatible() ITK_OVERRIDE { return false; }
        double GetApparentFractionalAnisotropy() ITK_OVERRIDE;

    protected:
        CHARMEDCompartment() : Superclass()
        {
            m_EstimateDiffusivities = false;
            m_EstimateTissueRadius = false;
            m_ChangedConstraints = true;

            m_VanGelderenCylinderCompartment = anima::VanGelderenCylinderCompartment::New();
            m_VanGelderenCylinderCompartment->SetEstimateTissueRadius(false);
            m_VanGelderenCylinderCompartment->SetEstimateAxialDiffusivity(false);

            m_StaniszCylinderCompartment = anima::StaniszCylinderCompartment::New();
            m_StaniszCylinderCompartment->SetEstimateTissueRadius(false);
            m_StaniszCylinderCompartment->SetEstimateAxialDiffusivity(false);

            m_NeumanCylinderCompartment = anima::NeumanCylinderCompartment::New();
            m_NeumanCylinderCompartment->SetEstimateTissueRadius(false);
            m_NeumanCylinderCompartment->SetEstimateAxialDiffusivity(false);
            m_NeumanCylinderCompartment->SetEchoTime(57e-3);

            m_ZeppelinCompartment = anima::ZeppelinCompartment::New();
            m_ZeppelinCompartment->SetEstimateDiffusivities(false);

            // Genu of corpus callosum (from Aboitiz et al. 1992)
            m_RadiusValues = {2.73e-4, 3.44e-4, 4.02e-4, 4.57e-4, 5.12e-4, 5.72e-4, 6.41e-4, 7.29e-4, 8.63e-4};
            m_RadiusWeights = {0.0943, 0.1281, 0.1428, 0.1455, 0.1390, 0.1251, 0.1044, 0.0773, 0.0435};
        }

        virtual ~CHARMEDCompartment() {}

    private:
        bool m_EstimateDiffusivities, m_EstimateTissueRadius;
        bool m_ChangedConstraints;
        unsigned int m_NumberOfParameters;

        // Parameters for data perparation
        anima::StaniszCylinderCompartment::Pointer m_StaniszCylinderCompartment;
        anima::NeumanCylinderCompartment::Pointer m_NeumanCylinderCompartment;
        anima::VanGelderenCylinderCompartment::Pointer m_VanGelderenCylinderCompartment;
        anima::ZeppelinCompartment::Pointer m_ZeppelinCompartment;

        std::vector<double> m_RadiusValues;
        std::vector<double> m_RadiusWeights;
    };

} // end namespace anima
