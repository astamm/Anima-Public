#include "animaCHARMEDCompartment.h"
#include <animaMCMConstants.h>

#include <limits>

namespace anima
{

    double CHARMEDCompartment::GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
    {
        double hinderedSignal = m_ZeppelinCompartment->GetFourierTransformedDiffusionProfile(smallDelta, bigDelta, gradientStrength, gradient);
        
        // double restrictedSignal = 0.0;
        // for (unsigned int i = 0;i < m_RadiusValues.size();++i)
        // {
        //     m_VanGelderenCylinderCompartment->SetTissueRadius(m_RadiusValues[i]);
        //     restrictedSignal += m_RadiusWeights[i] * m_VanGelderenCylinderCompartment->GetFourierTransformedDiffusionProfile(smallDelta, bigDelta, gradientStrength, gradient);
        // }
        double restrictedSignal = m_VanGelderenCylinderCompartment->GetFourierTransformedDiffusionProfile(smallDelta, bigDelta, gradientStrength, gradient);
        double fhValue = this->GetExtraAxonalFraction();
        return fhValue * hinderedSignal + (1.0 - fhValue) * restrictedSignal;
    }

    CHARMEDCompartment::ListType &CHARMEDCompartment::GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
    {
        m_JacobianVector.resize(this->GetNumberOfParameters());
        return m_JacobianVector;
    }

    double CHARMEDCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
    {
        return 0.0;
    }

    void CHARMEDCompartment::SetOrientationTheta(double num)
    {
        if (num != this->GetOrientationTheta())
        {
            m_VanGelderenCylinderCompartment->SetOrientationTheta(num);
            m_ZeppelinCompartment->SetOrientationTheta(num);
            this->Superclass::SetOrientationTheta(num);
        }
    }

    void CHARMEDCompartment::SetOrientationPhi(double num)
    {
        if (num != this->GetOrientationPhi())
        {
            m_VanGelderenCylinderCompartment->SetOrientationPhi(num);
            m_ZeppelinCompartment->SetOrientationPhi(num);
            this->Superclass::SetOrientationPhi(num);
        }
    }

    void CHARMEDCompartment::SetExtraAxonalFraction(double num)
    {
        if (num != this->GetExtraAxonalFraction())
        {
            this->Superclass::SetExtraAxonalFraction(num);
            this->SetRadialDiffusivity1(num);
        }
    }

    void CHARMEDCompartment::SetAxialDiffusivity(double num)
    {
        if (num != this->GetAxialDiffusivity())
        {
            m_VanGelderenCylinderCompartment->SetAxialDiffusivity(num);
            m_ZeppelinCompartment->SetAxialDiffusivity(num);
            this->Superclass::SetAxialDiffusivity(num);
        }
    }

    void CHARMEDCompartment::SetRadialDiffusivity1(double num)
    {
        if (!m_UseTortuosityModel)
        {
            m_ZeppelinCompartment->SetRadialDiffusivity1(num);
            this->Superclass::SetRadialDiffusivity1(num);
            return;
        }

        double axialDiff = this->GetAxialDiffusivity();
        double fr = 1.0 - this->GetExtraAxonalFraction();
        double radialDiff = axialDiff * (1.0 - 0.8 * fr);
        m_ZeppelinCompartment->SetRadialDiffusivity1(radialDiff);
        this->Superclass::SetRadialDiffusivity1(radialDiff);
    }

    void CHARMEDCompartment::SetTissueRadius(double num)
    {
        if (num != this->GetTissueRadius())
        {
            double radius = 0.0005;
            m_VanGelderenCylinderCompartment->SetTissueRadius(radius);
            this->Superclass::SetTissueRadius(radius);
        }
    }

    void CHARMEDCompartment::SetParametersFromVector(const ListType &params)
    {
        if (params.size() != this->GetNumberOfParameters())
            return;

        unsigned int pos = 0;

        this->SetOrientationTheta(params[pos]);
        ++pos;

        this->SetOrientationPhi(params[pos]);
        ++pos;

        this->SetExtraAxonalFraction(params[pos]);
        ++pos;

        if (m_EstimateTissueRadius)
        {
            this->SetTissueRadius(params[pos]);
            ++pos;
        }

        if (m_EstimateDiffusivities)
        {
            this->SetAxialDiffusivity(params[pos] + params[pos + 1]);
            ++pos;
            this->SetRadialDiffusivity1(params[pos]);
        }
    }

    CHARMEDCompartment::ListType &CHARMEDCompartment::GetParametersAsVector()
    {
        m_ParametersVector.resize(this->GetNumberOfParameters());

        unsigned int pos = 0;

        m_ParametersVector[pos] = this->GetOrientationTheta();
        ++pos;

        m_ParametersVector[pos] = this->GetOrientationPhi();
        ++pos;

        m_ParametersVector[pos] = this->GetExtraAxonalFraction();
        ++pos;

        if (m_EstimateTissueRadius)
        {
            m_ParametersVector[pos] = this->GetTissueRadius();
            ++pos;
        }

        if (m_EstimateDiffusivities)
        {
            m_ParametersVector[pos] = this->GetAxialDiffusivity() - this->GetRadialDiffusivity1();
            ++pos;
            m_ParametersVector[pos] = this->GetRadialDiffusivity1();
        }

        return m_ParametersVector;
    }

    CHARMEDCompartment::ListType &CHARMEDCompartment::GetParameterLowerBounds()
    {
        m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());

        unsigned int pos = 0;

        m_ParametersLowerBoundsVector[pos] = anima::MCMZeroLowerBound;
        ++pos;

        m_ParametersLowerBoundsVector[pos] = anima::MCMZeroLowerBound;
        ++pos;

        m_ParametersLowerBoundsVector[pos] = anima::MCMZeroLowerBound;
        ++pos;

        if (m_EstimateTissueRadius)
        {
            m_ParametersLowerBoundsVector[pos] = anima::MCMTissueRadiusLowerBound;
            ++pos;
        }

        if (m_EstimateDiffusivities)
        {
            m_ParametersLowerBoundsVector[pos] = anima::MCMZeroLowerBound;
            ++pos;
            m_ParametersLowerBoundsVector[pos] = anima::MCMDiffusivityLowerBound;
        }

        return m_ParametersLowerBoundsVector;
    }

    CHARMEDCompartment::ListType &CHARMEDCompartment::GetParameterUpperBounds()
    {
        m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

        unsigned int pos = 0;

        m_ParametersUpperBoundsVector[pos] = anima::MCMPolarAngleUpperBound;
        ++pos;

        m_ParametersUpperBoundsVector[pos] = anima::MCMAzimuthAngleUpperBound;
        ++pos;

        m_ParametersUpperBoundsVector[pos] = anima::MCMFractionUpperBound;
        ++pos;

        if (m_EstimateTissueRadius)
        {
            m_ParametersUpperBoundsVector[pos] = anima::MCMTissueRadiusUpperBound;
            ++pos;
        }

        if (m_EstimateDiffusivities)
        {
            m_ParametersUpperBoundsVector[pos] = anima::MCMDiffusivityUpperBound;
            ++pos;
            m_ParametersUpperBoundsVector[pos] = anima::MCMRadialDiffusivityUpperBound;
        }

        return m_ParametersUpperBoundsVector;
    }

    void CHARMEDCompartment::SetEstimateDiffusivities(bool arg)
    {
        if (m_EstimateDiffusivities == arg)
            return;

        m_VanGelderenCylinderCompartment->SetEstimateAxialDiffusivity(arg);
        m_ZeppelinCompartment->SetEstimateDiffusivities(arg);
        m_EstimateDiffusivities = arg;
        m_ChangedConstraints = true;
    }

    void CHARMEDCompartment::SetEstimateTissueRadius(bool arg)
    {
        if (m_EstimateTissueRadius == arg)
            return;

        m_VanGelderenCylinderCompartment->SetEstimateTissueRadius(arg);
        m_EstimateTissueRadius = arg;
        m_ChangedConstraints = true;
    }

    void CHARMEDCompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
    {
        if (compartmentVector.GetSize() != this->GetCompartmentSize())
            itkExceptionMacro("The input vector size does not match the size of the compartment");

        this->SetOrientationTheta(compartmentVector[0]);
        this->SetOrientationPhi(compartmentVector[1]);
        this->SetExtraAxonalFraction(compartmentVector[2]);
        this->SetTissueRadius(compartmentVector[3]);
        this->SetAxialDiffusivity(compartmentVector[4]);
        this->SetRadialDiffusivity1(compartmentVector[5]);
    }

    unsigned int CHARMEDCompartment::GetCompartmentSize()
    {
        return 6;
    }

    unsigned int CHARMEDCompartment::GetNumberOfParameters()
    {
        if (!m_ChangedConstraints)
            return m_NumberOfParameters;

        m_NumberOfParameters = 6;

        if (!m_EstimateTissueRadius)
            m_NumberOfParameters -= 1;

        if (!m_EstimateDiffusivities)
            m_NumberOfParameters -= 2;

        m_ChangedConstraints = false;
        return m_NumberOfParameters;
    }

    CHARMEDCompartment::ModelOutputVectorType &CHARMEDCompartment::GetCompartmentVector()
    {
        if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
            m_CompartmentVector.SetSize(this->GetCompartmentSize());

        // params are
        // 1) m_OrientationTheta
        // 2) m_OrientationPhi
        // 3) m_ExtraAxonalFraction
        // 4) m_TissueRadius
        // 5) m_AxialDiffusivity - m_RadialDiffusivity1 --> choose to fix it
        // 6) m_RadialDiffusivity1 --> choose to fix it

        m_CompartmentVector[0] = this->GetOrientationTheta();
        m_CompartmentVector[1] = this->GetOrientationPhi();
        m_CompartmentVector[2] = this->GetExtraAxonalFraction();
        m_CompartmentVector[3] = this->GetTissueRadius();
        m_CompartmentVector[4] = this->GetAxialDiffusivity();
        m_CompartmentVector[5] = this->GetRadialDiffusivity1();

        return m_CompartmentVector;
    }

    double CHARMEDCompartment::GetApparentFractionalAnisotropy()
    {
        return 0.0;
    }

} // end namespace anima
