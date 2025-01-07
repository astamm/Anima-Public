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
        void SetModelWithSmallGlialCellComponent(bool arg) { m_ModelWithSmallGlialCellComponent = arg; }
        void SetModelWithMediumGlialCellComponent(bool arg) { m_ModelWithMediumGlialCellComponent = arg; }
        void SetModelWithLargeGlialCellComponent(bool arg) { m_ModelWithLargeGlialCellComponent = arg; }
        void SetModelWithSphereComponent(bool arg) { m_ModelWithSphereComponent = arg; }

        void SetSphereCompartmentType(CompartmentType arg) { m_SphereCompartmentType = arg; }
        void SetCylinderCompartmentType(CompartmentType arg) { m_CylinderCompartmentType = arg; }
        void SetNumberOfCompartments(unsigned int num) { m_NumberOfCompartments = num; }
        void SetVariableProjectionEstimationMode(bool arg) { m_VariableProjectionEstimationMode = arg; }

        void SetUseConstrainedDiffusivity(bool arg) { m_UseConstrainedDiffusivity = arg; }
        void SetUseConstrainedOrientationConcentration(bool arg) { m_UseConstrainedOrientationConcentration = arg; }
        void SetUseConstrainedExtraAxonalFraction(bool arg) { m_UseConstrainedExtraAxonalFraction = arg; }
        void SetUseConstrainedSphereDiffusivity(bool arg) { m_UseConstrainedSphereDiffusivity = arg; }
        void SetUseConstrainedSphereRadius(bool arg) { m_UseConstrainedSphereRadius = arg; }
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
        void CreateSphereGPDPulsedGradientCompartment(BaseCompartmentPointer &compartmentPointer);
        void CreatePlaneSGPPulsedGradientCompartment(BaseCompartmentPointer &compartmentPointer);
        void CreateStickCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
        void CreateZeppelinCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
        void CreateTensorCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
        void CreateNODDICompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
        void CreateDDICompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);
        void CreateCHARMEDCompartment(BaseCompartmentPointer &compartmentPointer, bool applyConstraints);

        CompartmentType m_SphereCompartmentType, m_CylinderCompartmentType;
        bool m_ModelWithFreeWaterComponent;
        bool m_ModelWithSphereComponent;
        bool m_ModelWithSmallGlialCellComponent;
        bool m_ModelWithMediumGlialCellComponent;
        bool m_ModelWithLargeGlialCellComponent;
        unsigned int m_NumberOfCompartments;

        bool m_VariableProjectionEstimationMode;
        bool m_UseConstrainedDiffusivity;
        bool m_UseConstrainedOrientationConcentration;
        bool m_UseConstrainedExtraAxonalFraction;
        bool m_UseConstrainedSphereDiffusivity;
        bool m_UseConstrainedSphereRadius;
        bool m_UseConstrainedCylinderRadius;

        bool m_UseCommonDiffusivities;
        bool m_UseCommonConcentrations;
        bool m_UseCommonExtraAxonalFractions;

        double m_SphereDiffusivity;
        double m_OrientationConcentration, m_ExtraAxonalFraction;
        double m_AxialDiffusivity;
        double m_RadialDiffusivity1, m_RadialDiffusivity2;
        double m_SphereRadius, m_CylinderRadius;

        static const double m_FreeWaterDiffusivity;
        static const double m_SmallGlialCellRadius;
        static const double m_MediumGlialCellRadius;
        static const double m_LargeGlialCellRadius;
    };

} // end namespace anima
