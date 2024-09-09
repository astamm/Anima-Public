#include "animaNeumanCompartment.h"
#include <animaMCMConstants.h>

#include <limits>

namespace anima
{

    double NeumanCompartment::GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
    {
        double intrinsicDiffusivity = this->GetAxialDiffusivity();
        double radiusSq = this->GetTissueRadius() * this->GetTissueRadius();
        double surfValue = intrinsicDiffusivity * m_EchoTime;
        double inExpValue = 7.0 * radiusSq * radiusSq;
        inExpValue *= smallDelta * smallDelta;
        inExpValue *= anima::DiffusionGyromagneticRatio * anima::DiffusionGyromagneticRatio;
        inExpValue *= gradientStrength * gradientStrength;
        inExpValue /= (48.0 * surfValue);
        inExpValue *= std::abs(2.0 - 99.0 * radiusSq / (56.0 * surfValue));
        return std::exp(-inExpValue);
    }

    NeumanCompartment::ListType &NeumanCompartment::GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
    {
        m_JacobianVector.resize(this->GetNumberOfParameters());
        return m_JacobianVector;
    }

    double NeumanCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
    {
        return 0.0;
    }

    void NeumanCompartment::SetParametersFromVector(const ListType &params)
    {
        if (params.size() != this->GetNumberOfParameters())
            return;

        unsigned int pos = 0;
        if (m_EstimateTissueRadius)
        {
            this->SetTissueRadius(params[pos]);
            ++pos;
        }

        if (m_EstimateAxialDiffusivity)
            this->SetAxialDiffusivity(params[pos]);
    }

    NeumanCompartment::ListType &NeumanCompartment::GetParametersAsVector()
    {
        m_ParametersVector.resize(this->GetNumberOfParameters());

        unsigned int pos = 0;
        if (m_EstimateTissueRadius)
        {
            m_ParametersVector[pos] = this->GetTissueRadius();
            ++pos;
        }

        if (m_EstimateAxialDiffusivity)
            m_ParametersVector[pos] = this->GetAxialDiffusivity();

        return m_ParametersVector;
    }

    NeumanCompartment::ListType &NeumanCompartment::GetParameterLowerBounds()
    {
        m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());

        unsigned int pos = 0;
        if (m_EstimateTissueRadius)
        {
            m_ParametersLowerBoundsVector[pos] = anima::MCMTissueRadiusLowerBound;
            ++pos;
        }

        if (m_EstimateAxialDiffusivity)
            m_ParametersLowerBoundsVector[pos] = anima::MCMDiffusivityLowerBound;

        return m_ParametersLowerBoundsVector;
    }

    NeumanCompartment::ListType &NeumanCompartment::GetParameterUpperBounds()
    {
        m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

        unsigned int pos = 0;
        if (m_EstimateTissueRadius)
        {
            m_ParametersUpperBoundsVector[pos] = anima::MCMTissueRadiusUpperBound;
            ++pos;
        }

        if (m_EstimateAxialDiffusivity)
            m_ParametersUpperBoundsVector[pos] = anima::MCMDiffusivityUpperBound;

        return m_ParametersUpperBoundsVector;
    }

    void NeumanCompartment::SetEstimateAxialDiffusivity(bool arg)
    {
        if (m_EstimateAxialDiffusivity == arg)
            return;

        m_EstimateAxialDiffusivity = arg;
        m_ChangedConstraints = true;
    }

    void NeumanCompartment::SetEstimateTissueRadius(bool arg)
    {
        if (m_EstimateTissueRadius == arg)
            return;

        m_EstimateTissueRadius = arg;
        m_ChangedConstraints = true;
    }

    void NeumanCompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
    {
        if (compartmentVector.GetSize() != this->GetCompartmentSize())
            itkExceptionMacro("The input vector size does not match the size of the compartment");

        this->SetTissueRadius(compartmentVector[0]);
        this->SetAxialDiffusivity(compartmentVector[1]);
    }

    unsigned int NeumanCompartment::GetCompartmentSize()
    {
        return 2;
    }

    unsigned int NeumanCompartment::GetNumberOfParameters()
    {
        if (!m_ChangedConstraints)
            return m_NumberOfParameters;

        m_NumberOfParameters = 2;

        if (!m_EstimateTissueRadius)
            --m_NumberOfParameters;

        if (!m_EstimateAxialDiffusivity)
            --m_NumberOfParameters;

        m_ChangedConstraints = false;
        return m_NumberOfParameters;
    }

    NeumanCompartment::ModelOutputVectorType &NeumanCompartment::GetCompartmentVector()
    {
        if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
            m_CompartmentVector.SetSize(this->GetCompartmentSize());

        m_CompartmentVector[0] = this->GetTissueRadius();
        m_CompartmentVector[1] = this->GetAxialDiffusivity();

        return m_CompartmentVector;
    }

    double NeumanCompartment::GetApparentFractionalAnisotropy()
    {
        return 0.0;
    }

} // end namespace anima
