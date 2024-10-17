#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMExport.h>
#include <animaMultiCompartmentModel.h>

namespace anima
{

    //! Really this class is some simplified factory that creates the MCM that it knows
    // One could think of making it a singleton, but parameter are somewhat in the way
    class ANIMAMCM_EXPORT MultiCompartmentModelCreator
    {
    public:
        MultiCompartmentModelCreator();
        virtual ~MultiCompartmentModelCreator() {}

        using MCMType = anima::MultiCompartmentModel;
        using MCMPointer = MCMType::Pointer;

        using BaseCompartmentType = anima::BaseCompartment;
        using BaseCompartmentPointer = BaseCompartmentType::Pointer;
        using CompartmentType = anima::DiffusionModelCompartmentType;

        // Isotropic compartments (one free water, one sphere)
        void SetModelWithFreeWaterComponent(bool arg) { m_ModelWithFreeWaterComponent = arg; }
        void SetModelWithSphereComponent(bool arg) { m_ModelWithSphereComponent = arg; }

        void SetCompartmentType(CompartmentType arg) { m_CompartmentType = arg; }
        void SetNumberOfCompartments(unsigned int num) { m_NumberOfCompartments = num; }
        void SetVariableProjectionEstimationMode(bool arg) { m_VariableProjectionEstimationMode = arg; }

        void SetUseConstrainedDiffusivity(bool arg) { m_UseConstrainedDiffusivity = arg; }
        void SetUseConstrainedOrientationConcentration(bool arg) { m_UseConstrainedOrientationConcentration = arg; }
        void SetUseConstrainedExtraAxonalFraction(bool arg) { m_UseConstrainedExtraAxonalFraction = arg; }
        void SetUseConstrainedFreeWaterDiffusivity(bool arg) { m_UseConstrainedFreeWaterDiffusivity = arg; }
        void SetUseConstrainedStaniszDiffusivity(bool arg) { m_UseConstrainedStaniszDiffusivity = arg; }
        void SetUseConstrainedStaniszRadius(bool arg) { m_UseConstrainedStaniszRadius = arg; }
        void SetUseConstrainedCylinderRadius(bool arg) { m_UseConstrainedCylinderRadius = arg; }

        bool GetUseConstrainedDiffusivity() { return m_UseConstrainedDiffusivity; }
        bool GetUseConstrainedOrientationConcentration() { return m_UseConstrainedOrientationConcentration; }
        bool GetUseConstrainedExtraAxonalFraction() { return m_UseConstrainedExtraAxonalFraction; }

        void SetUseCommonDiffusivities(bool arg) { m_UseCommonDiffusivities = arg; }
        void SetUseCommonConcentrations(bool arg) { m_UseCommonConcentrations = arg; }
        void SetUseCommonExtraAxonalFractions(bool arg) { m_UseCommonExtraAxonalFractions = arg; }

        bool GetUseCommonDiffusivities() { return m_UseCommonDiffusivities; }
        bool GetUseCommonConcentrations() { return m_UseCommonConcentrations; }
        bool GetUseCommonExtraAxonalFractions() { return m_UseCommonExtraAxonalFractions; }

        void SetFreeWaterDiffusivityValue(double arg) { m_FreeWaterDiffusivity = arg; }
        void SetSphereDiffusivityValue(double arg) { m_SphereDiffusivity = arg; }
        void SetAxialDiffusivityValue(double arg) { m_AxialDiffusivity = arg; }
        void SetRadialDiffusivity1Value(double arg) { m_RadialDiffusivity1 = arg; }
        void SetRadialDiffusivity2Value(double arg) { m_RadialDiffusivity2 = arg; }
        void SetSphereRadiusValue(double arg) { m_SphereRadius = arg; }
        void SetCylinderRadiusValue(double arg) { m_CylinderRadius = arg; }

        double GetOrientationConcentration() { return m_OrientationConcentration; }
        double GetAxialDiffusivity() { return m_AxialDiffusivity; }
        double GetRadialDiffusivity1() { return m_RadialDiffusivity1; }
        double GetRadialDiffusivity2() { return m_RadialDiffusivity2; }
        double GetExtraAxonalFraction() { return m_ExtraAxonalFraction; }
        double GetSphereRadius() { return m_SphereRadius; }
        double GetCylinderRadius() { return m_CylinderRadius; }

        MCMPointer GetNewMultiCompartmentModel();

    private:
        void CreateStickCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
        void CreateZeppelinCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
        void CreateTensorCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
        void CreateNODDICompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
        void CreateDDICompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
        void CreateCHARMEDCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);

        CompartmentType m_CompartmentType;
        bool m_ModelWithFreeWaterComponent;
        bool m_ModelWithSphereComponent;
        unsigned int m_NumberOfCompartments;

        bool m_VariableProjectionEstimationMode;
        bool m_UseConstrainedDiffusivity;
        bool m_UseConstrainedOrientationConcentration;
        bool m_UseConstrainedExtraAxonalFraction;
        bool m_UseConstrainedFreeWaterDiffusivity;
        bool m_UseConstrainedStaniszDiffusivity;
        bool m_UseConstrainedStaniszRadius;
        bool m_UseConstrainedCylinderRadius;

        bool m_UseCommonDiffusivities;
        bool m_UseCommonConcentrations;
        bool m_UseCommonExtraAxonalFractions;

        double m_FreeWaterDiffusivity, m_SphereDiffusivity;
        double m_OrientationConcentration, m_ExtraAxonalFraction;
        double m_AxialDiffusivity;
        double m_RadialDiffusivity1, m_RadialDiffusivity2;
        double m_SphereRadius, m_CylinderRadius;
    };

} // end namespace anima
