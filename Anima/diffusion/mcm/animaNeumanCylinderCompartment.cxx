#include "animaNeumanCylinderCompartment.h"
#include <animaMCMConstants.h>
#include <animaVectorOperations.h>

#include <limits>

namespace anima
{

    double NeumanCylinderCompartment::GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
    {
        double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, bigDelta, gradientStrength);
        Vector3DType sphCoords, carCoords;
        sphCoords[0] = this->GetOrientationTheta();
        sphCoords[1] = this->GetOrientationPhi();
        sphCoords[2] = 1.0;
        anima::TransformSphericalToCartesianCoordinates(sphCoords, carCoords);
        double cosAngle = anima::ComputeScalarProduct(gradient, carCoords);
        double axialSignal = std::exp(-bValue * this->GetAxialDiffusivity() * cosAngle * cosAngle);
        double radialSignal = m_NeumanCompartment->GetFourierTransformedDiffusionProfile(smallDelta, bigDelta, gradientStrength * std::sqrt(1.0 - cosAngle * cosAngle), gradient);
        return axialSignal * radialSignal;
    }

    NeumanCylinderCompartment::ListType &NeumanCylinderCompartment::GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
    {
        m_JacobianVector.resize(this->GetNumberOfParameters());

        return m_JacobianVector;
    }

    double NeumanCylinderCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
    {
        return 0.0;
    }

    void NeumanCylinderCompartment::SetTissueRadius(double num)
    {
        if (num != this->GetTissueRadius())
        {
            m_NeumanCompartment->SetTissueRadius(num);
            this->Superclass::SetTissueRadius(num);
        }
    }

    void NeumanCylinderCompartment::SetAxialDiffusivity(double num)
    {
        if (num != this->GetAxialDiffusivity())
        {
            m_NeumanCompartment->SetAxialDiffusivity(num);
            this->Superclass::SetAxialDiffusivity(num);
        }
    }

    void NeumanCylinderCompartment::SetParametersFromVector(const ListType &params)
    {
        if (params.size() != this->GetNumberOfParameters())
            return;

        unsigned int pos = 0;

        this->SetOrientationTheta(params[pos]);
        ++pos;

        this->SetOrientationPhi(params[pos]);
        ++pos;

        if (m_EstimateTissueRadius)
        {
            this->SetTissueRadius(params[pos]);
            ++pos;
        }

        if (m_EstimateAxialDiffusivity)
            this->SetAxialDiffusivity(params[pos]);
    }

    NeumanCylinderCompartment::ListType &NeumanCylinderCompartment::GetParametersAsVector()
    {
        m_ParametersVector.resize(this->GetNumberOfParameters());

        unsigned int pos = 0;

        m_ParametersVector[pos] = this->GetOrientationTheta();
        ++pos;

        m_ParametersVector[pos] = this->GetOrientationPhi();
        ++pos;

        if (m_EstimateTissueRadius)
        {
            m_ParametersVector[pos] = this->GetTissueRadius();
            ++pos;
        }

        if (m_EstimateAxialDiffusivity)
            m_ParametersVector[pos] = this->GetAxialDiffusivity();

        return m_ParametersVector;
    }

    NeumanCylinderCompartment::ListType &NeumanCylinderCompartment::GetParameterLowerBounds()
    {
        m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());

        unsigned int pos = 0;

        m_ParametersLowerBoundsVector[pos] = anima::MCMZeroLowerBound;
        ++pos;

        m_ParametersLowerBoundsVector[pos] = anima::MCMZeroLowerBound;
        ++pos;

        if (m_EstimateTissueRadius)
        {
            m_ParametersLowerBoundsVector[pos] = anima::MCMTissueRadiusLowerBound;
            ++pos;
        }

        if (m_EstimateAxialDiffusivity)
            m_ParametersLowerBoundsVector[pos] = anima::MCMDiffusivityLowerBound;

        return m_ParametersLowerBoundsVector;
    }

    NeumanCylinderCompartment::ListType &NeumanCylinderCompartment::GetParameterUpperBounds()
    {
        m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

        unsigned int pos = 0;

        m_ParametersUpperBoundsVector[pos] = anima::MCMPolarAngleUpperBound;
        ++pos;

        m_ParametersUpperBoundsVector[pos] = anima::MCMAzimuthAngleUpperBound;
        ++pos;

        if (m_EstimateTissueRadius)
        {
            m_ParametersUpperBoundsVector[pos] = anima::MCMTissueRadiusUpperBound;
            ++pos;
        }

        if (m_EstimateAxialDiffusivity)
            m_ParametersUpperBoundsVector[pos] = anima::MCMDiffusivityUpperBound;

        return m_ParametersUpperBoundsVector;
    }

    void NeumanCylinderCompartment::SetEstimateAxialDiffusivity(bool arg)
    {
        if (m_EstimateAxialDiffusivity == arg)
            return;

        m_NeumanCompartment->SetEstimateAxialDiffusivity(arg);
        m_EstimateAxialDiffusivity = arg;
        m_ChangedConstraints = true;
    }

    void NeumanCylinderCompartment::SetEstimateTissueRadius(bool arg)
    {
        if (m_EstimateTissueRadius == arg)
            return;

        m_NeumanCompartment->SetEstimateTissueRadius(arg);
        m_EstimateTissueRadius = arg;
        m_ChangedConstraints = true;
    }

    void NeumanCylinderCompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
    {
        if (compartmentVector.GetSize() != this->GetCompartmentSize())
            itkExceptionMacro("The input vector size does not match the size of the compartment");

        this->SetOrientationTheta(compartmentVector[0]);
        this->SetOrientationPhi(compartmentVector[1]);
        this->SetTissueRadius(compartmentVector[2]);
        this->SetAxialDiffusivity(compartmentVector[3]);
    }

    unsigned int NeumanCylinderCompartment::GetCompartmentSize()
    {
        return 4;
    }

    unsigned int NeumanCylinderCompartment::GetNumberOfParameters()
    {
        if (!m_ChangedConstraints)
            return m_NumberOfParameters;

        m_NumberOfParameters = 4;

        if (!m_EstimateTissueRadius)
            m_NumberOfParameters--;

        if (!m_EstimateAxialDiffusivity)
            m_NumberOfParameters--;

        m_ChangedConstraints = false;
        return m_NumberOfParameters;
    }

    NeumanCylinderCompartment::ModelOutputVectorType &NeumanCylinderCompartment::GetCompartmentVector()
    {
        if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
            m_CompartmentVector.SetSize(this->GetCompartmentSize());

        // params are
        // 1) m_OrientationTheta
        // 2) m_OrientationPhi
        // 3) m_TissueRadius --> choose to fix it
        // 4) m_AxialDiffusivity --> choose to fix it

        m_CompartmentVector[0] = this->GetOrientationTheta();
        m_CompartmentVector[1] = this->GetOrientationPhi();
        m_CompartmentVector[2] = this->GetTissueRadius();
        m_CompartmentVector[3] = this->GetAxialDiffusivity();

        return m_CompartmentVector;
    }

    double NeumanCylinderCompartment::GetApparentFractionalAnisotropy()
    {
        return 0.0;
    }

} // end namespace anima
