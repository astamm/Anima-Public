#include "animaMultiCompartmentModelCreator.h"

#include <animaCHARMEDCompartment.h>
#include <animaDDICompartment.h>
#include <animaFreeWaterCompartment.h>
#include <animaNODDICompartment.h>
#include <animaSphereGPDPulsedGradientCompartment.h>
#include <animaPlaneSGPPulsedGradientCompartment.h>
#include <animaStickCompartment.h>
#include <animaTensorCompartment.h>
#include <animaZeppelinCompartment.h>

namespace anima
{
    const double MultiCompartmentModelCreator::m_FreeWaterDiffusivity = 3.0e-3;
    const double MultiCompartmentModelCreator::m_SmallGlialCellRadius = 5e-3;
    const double MultiCompartmentModelCreator::m_MediumGlialCellRadius = 25e-3;
    const double MultiCompartmentModelCreator::m_LargeGlialCellRadius = 50e-3;

    MultiCompartmentModelCreator::MultiCompartmentModelCreator()
    {
        m_SphereCompartmentType = SphereGPDPulsedGradient;
        m_CylinderCompartmentType = Tensor;
        m_ModelWithFreeWaterComponent = false;
        m_ModelWithSphereComponent = false;

        m_NumberOfCompartments = 1;
        m_VariableProjectionEstimationMode = true;

        m_UseConstrainedDiffusivity = false;
        m_UseConstrainedSphereDiffusivity = true;
        m_UseConstrainedSphereRadius = true;
        m_UseConstrainedCylinderRadius = true;
        m_UseConstrainedOrientationConcentration = false;
        m_UseConstrainedExtraAxonalFraction = false;

        m_UseCommonDiffusivities = false;
        m_UseCommonConcentrations = false;
        m_UseCommonExtraAxonalFractions = false;

        m_AxialDiffusivity = 1.71e-3;
        m_SphereDiffusivity = 1.71e-3;
        m_RadialDiffusivity1 = 1.9e-4;
        m_RadialDiffusivity2 = 1.5e-4;
        m_ExtraAxonalFraction = 0.1;
        m_OrientationConcentration = 10.0;
        m_SphereRadius = 0.015;
        m_CylinderRadius = 0.0005;
    }

    MultiCompartmentModelCreator::MCMPointer MultiCompartmentModelCreator::GetNewMultiCompartmentModel()
    {
        MCMPointer outputMCM = MCMType::New();
        outputMCM->SetOptimizeWeights(!m_VariableProjectionEstimationMode);
        outputMCM->SetCommonDiffusivityParameters(m_UseCommonDiffusivities);
        outputMCM->SetCommonConcentrationParameters(m_UseCommonConcentrations);
        outputMCM->SetCommonExtraAxonalFractionParameters(m_UseCommonExtraAxonalFractions);

        double numCompartments = m_ModelWithFreeWaterComponent + m_ModelWithSphereComponent + m_NumberOfCompartments;
        double defaultWeight = 1.0 / numCompartments;

        if (m_ModelWithFreeWaterComponent)
        {
            using FreeWaterType = anima::FreeWaterCompartment;
            FreeWaterType::Pointer fwComp = FreeWaterType::New();
            fwComp->SetEstimateAxialDiffusivity(false);
            fwComp->SetAxialDiffusivity(m_FreeWaterDiffusivity);

            outputMCM->AddCompartment(defaultWeight, fwComp);
        }

        if (m_ModelWithSmallGlialCellComponent)
        {
            using SphereType = anima::SphereGPDPulsedGradientCompartment;
            SphereType::Pointer sphereComp = SphereType::New();
            sphereComp->SetEstimateAxialDiffusivity(false);
            sphereComp->SetEstimateTissueRadius(false);
            sphereComp->SetAxialDiffusivity(m_SphereDiffusivity);
            sphereComp->SetTissueRadius(m_SmallGlialCellRadius);

            outputMCM->AddCompartment(defaultWeight, sphereComp);
        }

        if (m_ModelWithMediumGlialCellComponent)
        {
            using SphereType = anima::SphereGPDPulsedGradientCompartment;
            SphereType::Pointer sphereComp = SphereType::New();
            sphereComp->SetEstimateAxialDiffusivity(false);
            sphereComp->SetEstimateTissueRadius(false);
            sphereComp->SetAxialDiffusivity(m_SphereDiffusivity);
            sphereComp->SetTissueRadius(m_MediumGlialCellRadius);

            outputMCM->AddCompartment(defaultWeight, sphereComp);
        }

        if (m_ModelWithLargeGlialCellComponent)
        {
            using SphereType = anima::SphereGPDPulsedGradientCompartment;
            SphereType::Pointer sphereComp = SphereType::New();
            sphereComp->SetEstimateAxialDiffusivity(false);
            sphereComp->SetEstimateTissueRadius(false);
            sphereComp->SetAxialDiffusivity(m_SphereDiffusivity);
            sphereComp->SetTissueRadius(m_LargeGlialCellRadius);

            outputMCM->AddCompartment(defaultWeight, sphereComp);
        }

        if (m_ModelWithSphereComponent)
        {
            anima::BaseCompartment::Pointer sphereComp;

            switch (m_SphereCompartmentType)
            {
                case SphereGPDPulsedGradient:
                    this->CreateSphereGPDPulsedGradientCompartment(sphereComp);
                    break;
                
                case PlaneSGPPulsedGradient:
                    this->CreatePlaneSGPPulsedGradientCompartment(sphereComp);
                    break;
                
                default:
                    throw itk::ExceptionObject(__FILE__, __LINE__, "The sphere compartment must be one of SphereGPDPulsedGradient or PlaneSGPPulsedGradient.", ITK_LOCATION);
                    break;
            }

            outputMCM->AddCompartment(defaultWeight, sphereComp);
        }

        for (unsigned int i = 0; i < m_NumberOfCompartments; ++i)
        {
            anima::BaseCompartment::Pointer tmpPointer;
            bool applyCommonConstraints = (i > 0);

            switch (m_CylinderCompartmentType)
            {
            case Stick:
                this->CreateStickCompartment(tmpPointer, applyCommonConstraints);
                break;

            case Zeppelin:
                this->CreateZeppelinCompartment(tmpPointer, applyCommonConstraints);
                break;

            case Tensor:
                this->CreateTensorCompartment(tmpPointer, applyCommonConstraints);
                break;

            case NODDI:
                this->CreateNODDICompartment(tmpPointer, applyCommonConstraints);
                break;

            case DDI:
                this->CreateDDICompartment(tmpPointer, applyCommonConstraints);
                break;

            case CHARMED:
                this->CreateCHARMEDCompartment(tmpPointer, applyCommonConstraints);
                break;

            default:
                throw itk::ExceptionObject(__FILE__, __LINE__, "Cylinder compartments must be one of Stick, Zeppelin, Tensor, NODDI, DDI or CHARMED.", ITK_LOCATION);
                break;
            }

            // Kind of ugly but required for optimization, otherwise initialization from simplified models may fail
            tmpPointer->SetOrientationConcentration(m_OrientationConcentration);

            outputMCM->AddCompartment(defaultWeight, tmpPointer);
        }

        return outputMCM;
    }

    void MultiCompartmentModelCreator::CreateSphereGPDPulsedGradientCompartment(BaseCompartmentPointer &compartmentPointer)
    {
        using SphereType = anima::SphereGPDPulsedGradientCompartment;
        SphereType::Pointer sphereComp = SphereType::New();
        sphereComp->SetEstimateAxialDiffusivity(!m_UseConstrainedSphereDiffusivity);
        sphereComp->SetEstimateTissueRadius(!m_UseConstrainedSphereRadius);
        sphereComp->SetAxialDiffusivity(m_SphereDiffusivity);
        sphereComp->SetTissueRadius(m_SphereRadius);
        compartmentPointer = sphereComp;
    }

    void MultiCompartmentModelCreator::CreatePlaneSGPPulsedGradientCompartment(BaseCompartmentPointer &compartmentPointer)
    {
        using SphereType = anima::PlaneSGPPulsedGradientCompartment;
        SphereType::Pointer sphereComp = SphereType::New();
        sphereComp->SetEstimateAxialDiffusivity(!m_UseConstrainedSphereDiffusivity);
        sphereComp->SetEstimateTissueRadius(!m_UseConstrainedSphereRadius);
        sphereComp->SetAxialDiffusivity(m_SphereDiffusivity);
        sphereComp->SetTissueRadius(m_SphereRadius);
        compartmentPointer = sphereComp;
    }

    void MultiCompartmentModelCreator::CreateStickCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints)
    {
        using StickType = anima::StickCompartment;

        StickType::Pointer stickComp = StickType::New();
        stickComp->SetEstimateAxialDiffusivity(!m_UseConstrainedDiffusivity);

        stickComp->SetAxialDiffusivity(m_AxialDiffusivity);
        stickComp->SetRadialDiffusivity1((m_RadialDiffusivity1 + m_RadialDiffusivity2) / 2.0);

        if (applyConstraints)
        {
            if (m_UseCommonDiffusivities)
                stickComp->SetEstimateAxialDiffusivity(false);
        }

        compartmentPointer = stickComp;
    }

    void MultiCompartmentModelCreator::CreateZeppelinCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints)
    {
        using ZeppelinType = anima::ZeppelinCompartment;

        ZeppelinType::Pointer zepComp = ZeppelinType::New();
        zepComp->SetEstimateDiffusivities(!m_UseConstrainedDiffusivity);

        zepComp->SetAxialDiffusivity(m_AxialDiffusivity);
        zepComp->SetRadialDiffusivity1((m_RadialDiffusivity1 + m_RadialDiffusivity2) / 2.0);

        if (applyConstraints)
        {
            if (m_UseCommonDiffusivities)
                zepComp->SetEstimateDiffusivities(false);
        }

        compartmentPointer = zepComp;
    }

    void MultiCompartmentModelCreator::CreateTensorCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints)
    {
        using TensorType = anima::TensorCompartment;

        TensorType::Pointer tensComp = TensorType::New();
        tensComp->SetEstimateDiffusivities(!m_UseConstrainedDiffusivity);

        tensComp->SetAxialDiffusivity(m_AxialDiffusivity);
        tensComp->SetRadialDiffusivity1(m_RadialDiffusivity1);
        tensComp->SetRadialDiffusivity2(m_RadialDiffusivity2);

        if (applyConstraints)
        {
            if (m_UseCommonDiffusivities)
                tensComp->SetEstimateDiffusivities(false);
        }

        compartmentPointer = tensComp;
    }

    void MultiCompartmentModelCreator::CreateNODDICompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints)
    {
        using NODDIType = anima::NODDICompartment;

        NODDIType::Pointer noddiComp = NODDIType::New();
        noddiComp->SetEstimateAxialDiffusivity(!m_UseConstrainedDiffusivity);

        noddiComp->SetOrientationConcentration(m_OrientationConcentration);
        noddiComp->SetExtraAxonalFraction(m_ExtraAxonalFraction);
        noddiComp->SetAxialDiffusivity(m_AxialDiffusivity);

        if (applyConstraints)
        {
            if (m_UseCommonDiffusivities)
                noddiComp->SetEstimateAxialDiffusivity(false);

            if (this->GetUseCommonConcentrations())
                noddiComp->SetEstimateOrientationConcentration(false);

            if (this->GetUseCommonExtraAxonalFractions())
                noddiComp->SetEstimateExtraAxonalFraction(false);
        }

        compartmentPointer = noddiComp;
    }

    void MultiCompartmentModelCreator::CreateDDICompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints)
    {
        using DDIType = anima::DDICompartment;

        DDIType::Pointer ddiComp = DDIType::New();
        ddiComp->SetEstimateOrientationConcentration(!this->GetUseConstrainedOrientationConcentration());
        ddiComp->SetEstimateAxialDiffusivity(!this->GetUseConstrainedDiffusivity());
        ddiComp->SetEstimateExtraAxonalFraction(!this->GetUseConstrainedExtraAxonalFraction());

        ddiComp->SetOrientationConcentration(this->GetOrientationConcentration());
        ddiComp->SetAxialDiffusivity(this->GetAxialDiffusivity());
        ddiComp->SetRadialDiffusivity1((this->GetRadialDiffusivity1() + this->GetRadialDiffusivity2()) / 2.0);
        ddiComp->SetExtraAxonalFraction(this->GetExtraAxonalFraction());

        if (applyConstraints)
        {
            if (this->GetUseCommonDiffusivities())
                ddiComp->SetEstimateAxialDiffusivity(false);

            if (this->GetUseCommonConcentrations())
                ddiComp->SetEstimateOrientationConcentration(false);

            if (this->GetUseCommonExtraAxonalFractions())
                ddiComp->SetEstimateExtraAxonalFraction(false);
        }

        compartmentPointer = ddiComp;
    }

    void MultiCompartmentModelCreator::CreateCHARMEDCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints)
    {
        using CHARMEDType = anima::CHARMEDCompartment;

        CHARMEDType::Pointer charmedComp = CHARMEDType::New();
        charmedComp->SetEstimateDiffusivities(!m_UseConstrainedDiffusivity);
        charmedComp->SetEstimateTissueRadius(!m_UseConstrainedCylinderRadius);
        charmedComp->SetAxialDiffusivity(m_AxialDiffusivity);
        charmedComp->SetRadialDiffusivity1((m_RadialDiffusivity1 + m_RadialDiffusivity2) / 2.0);
        charmedComp->SetTissueRadius(m_CylinderRadius);

        if (applyConstraints)
        {
            if (m_UseCommonDiffusivities)
                charmedComp->SetEstimateDiffusivities(false);
        }

        compartmentPointer = charmedComp;
    }

} // end namespace anima
