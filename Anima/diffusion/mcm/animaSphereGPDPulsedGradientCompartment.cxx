#include "animaSphereGPDPulsedGradientCompartment.h"
#include <animaMCMConstants.h>

#include <limits>

namespace anima
{

    SphereGPDPulsedGradientCompartment::KeyType SphereGPDPulsedGradientCompartment::GenerateKey(double smallDelta, double bigDelta, double gradientStrength)
    {
        unsigned int smallDeltaInt = std::floor(smallDelta * 1.0e6);
        unsigned int bigDeltaInt = std::floor(bigDelta * 1.0e6);
        unsigned int gradientStrengthInt = std::floor(gradientStrength * 1.0e6);

        KeyType outValue(smallDeltaInt, bigDeltaInt, gradientStrengthInt);

        return outValue;
    }

    void SphereGPDPulsedGradientCompartment::UpdateSignals(double smallDelta, double bigDelta, double gradientStrength)
    {
        KeyType gradientKey = this->GenerateKey(smallDelta, bigDelta, gradientStrength);
        if (m_FirstSummations.find(gradientKey) != m_FirstSummations.end())
            return;

        double firstSummation = 0.0, secondSummation = 0.0, thirdSummation = 0.0, fourthSummation = 0.0;
        for (unsigned int m = 0; m < m_MaximumNumberOfSumElements; ++m)
        {
            double sqRoot = m_FunctorZeros[m] * m_FunctorZeros[m];
            double alphaSq = sqRoot / (this->GetTissueRadius() * this->GetTissueRadius());
            double workValue = this->GetAxialDiffusivity() * alphaSq;
            double alpha4Inv = 1.0 / (alphaSq * alphaSq);
            double numElement = 2.0;
            numElement -= 2.0 * std::exp(-workValue * smallDelta);
            numElement -= 2.0 * std::exp(-workValue * bigDelta);
            numElement += std::exp(-workValue * (bigDelta - smallDelta));
            numElement += std::exp(-workValue * (bigDelta + smallDelta));
            double sumElement = 2.0 * smallDelta - numElement / workValue;
            sumElement *= alpha4Inv / (sqRoot - 2.0);
            firstSummation += sumElement;

            if (std::abs(sumElement) < m_SignalSummationTolerance * std::abs(firstSummation))
                break;
        }

        if (std::abs(firstSummation) < std::numeric_limits<double>::epsilon())
            firstSummation = 0;
        if (std::abs(secondSummation) < std::numeric_limits<double>::epsilon())
            secondSummation = 0;
        if (std::abs(thirdSummation) < std::numeric_limits<double>::epsilon())
            thirdSummation = 0;
        if (std::abs(fourthSummation) < std::numeric_limits<double>::epsilon())
            fourthSummation = 0;

        m_FirstSummations[gradientKey] = firstSummation;
        m_SecondSummations[gradientKey] = secondSummation;
        m_ThirdSummations[gradientKey] = thirdSummation;
        m_FourthSummations[gradientKey] = fourthSummation;
    }

    double SphereGPDPulsedGradientCompartment::GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
    {
        this->UpdateSignals(smallDelta, bigDelta, gradientStrength);
        KeyType gradientKey = this->GenerateKey(smallDelta, bigDelta, gradientStrength);

        double logSignalValue = -2.0 * anima::DiffusionGyromagneticRatio * anima::DiffusionGyromagneticRatio;
        logSignalValue /= this->GetAxialDiffusivity();
        logSignalValue *= (gradientStrength * gradientStrength);
        logSignalValue *= m_FirstSummations[gradientKey];

        return std::exp(logSignalValue);
    }

    SphereGPDPulsedGradientCompartment::ListType &SphereGPDPulsedGradientCompartment::GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
    {
        this->UpdateSignals(smallDelta, bigDelta, gradientStrength);
        KeyType gradientKey = this->GenerateKey(smallDelta, bigDelta, gradientStrength);

        m_JacobianVector.resize(this->GetNumberOfParameters());

        double tissueRadius = this->GetTissueRadius();
        double axialDiff = this->GetAxialDiffusivity();
        double alpha = anima::DiffusionGyromagneticRatio * smallDelta * gradientStrength;
        double alphaRs = alpha * tissueRadius;
        double alphaRsSquare = alphaRs * alphaRs;
        double piSquare = M_PI * M_PI;
        double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, bigDelta, gradientStrength);

        // Derivative w.r.t. R_S
        unsigned int pos = 0;

        if (m_EstimateTissueRadius)
        {
            if (alphaRs < 1.0e-8)
                m_JacobianVector[pos] = -alphaRsSquare / (6.0 * tissueRadius);
            else
                m_JacobianVector[pos] = 2.0 * (alphaRs * std::sin(alphaRs) - 2.0 * (1.0 - std::cos(alphaRs))) / (alphaRsSquare * tissueRadius);

            m_JacobianVector[pos] += 8.0 * alpha * alphaRs * m_FirstSummations[gradientKey];
            m_JacobianVector[pos] += 8.0 * piSquare * bValue * axialDiff * m_SecondSummations[gradientKey] / tissueRadius;
            m_JacobianVector[pos] += 4.0 * alpha * alphaRsSquare * m_ThirdSummations[gradientKey];
            m_JacobianVector[pos] -= 16.0 * alpha * alphaRs * alphaRsSquare * m_FourthSummations[gradientKey];

            ++pos;
        }

        if (m_EstimateAxialDiffusivity)
        {
            // Derivative w.r.t. D_I
            m_JacobianVector[pos] = -4.0 * bValue * piSquare * m_SecondSummations[gradientKey];
        }

        return m_JacobianVector;
    }

    double SphereGPDPulsedGradientCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
    {
        // Compute equivalent isotropic term for default delta values and gradient strength
        double bValue = 1000.0;
        double gradientStrength = anima::GetGradientStrengthFromBValue(bValue, anima::DiffusionSmallDelta, anima::DiffusionBigDelta);
        this->UpdateSignals(anima::DiffusionSmallDelta, anima::DiffusionBigDelta, gradientStrength);
        KeyType gradientKey = this->GenerateKey(anima::DiffusionSmallDelta, anima::DiffusionBigDelta, gradientStrength);

        double alpha = anima::DiffusionGyromagneticRatio * anima::DiffusionSmallDelta * gradientStrength;
        double alphaRs = alpha * this->GetTissueRadius();
        double alphaRsSquare = alphaRs * alphaRs;
        double signalValue = 4.0 * alphaRsSquare * m_FirstSummations[gradientKey];

        if (alphaRs < 1.0e-8)
            signalValue += 1.0 - alphaRsSquare / 12.0;
        else
            signalValue += 2.0 * (1.0 - std::cos(alphaRs)) / alphaRsSquare;

        double equivalentGaussianDiffusivity = -std::log(signalValue) / bValue;

        double resVal = -1.5 * std::log(2.0 * M_PI * equivalentGaussianDiffusivity);

        resVal -= sample.squared_magnitude() / (2.0 * equivalentGaussianDiffusivity);

        return resVal;
    }

    void SphereGPDPulsedGradientCompartment::SetTissueRadius(double num)
    {
        if (num != this->GetTissueRadius())
        {
            m_FirstSummations.clear();
            m_SecondSummations.clear();
            m_ThirdSummations.clear();
            m_FourthSummations.clear();
            this->Superclass::SetTissueRadius(num);
        }
    }

    void SphereGPDPulsedGradientCompartment::SetAxialDiffusivity(double num)
    {
        if (num != this->GetAxialDiffusivity())
        {
            m_FirstSummations.clear();
            m_SecondSummations.clear();
            m_ThirdSummations.clear();
            m_FourthSummations.clear();
            this->Superclass::SetAxialDiffusivity(num);
        }
    }

    void SphereGPDPulsedGradientCompartment::SetParametersFromVector(const ListType &params)
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

    SphereGPDPulsedGradientCompartment::ListType &SphereGPDPulsedGradientCompartment::GetParametersAsVector()
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

    SphereGPDPulsedGradientCompartment::ListType &SphereGPDPulsedGradientCompartment::GetParameterLowerBounds()
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

    SphereGPDPulsedGradientCompartment::ListType &SphereGPDPulsedGradientCompartment::GetParameterUpperBounds()
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

    void SphereGPDPulsedGradientCompartment::SetEstimateAxialDiffusivity(bool arg)
    {
        if (m_EstimateAxialDiffusivity == arg)
            return;

        m_EstimateAxialDiffusivity = arg;
        m_ChangedConstraints = true;
    }

    void SphereGPDPulsedGradientCompartment::SetEstimateTissueRadius(bool arg)
    {
        if (m_EstimateTissueRadius == arg)
            return;

        m_EstimateTissueRadius = arg;
        m_ChangedConstraints = true;
    }

    void SphereGPDPulsedGradientCompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
    {
        if (compartmentVector.GetSize() != this->GetCompartmentSize())
            itkExceptionMacro("The input vector size does not match the size of the compartment");

        this->SetTissueRadius(compartmentVector[0]);
        this->SetAxialDiffusivity(compartmentVector[1]);
    }

    unsigned int SphereGPDPulsedGradientCompartment::GetCompartmentSize()
    {
        return 2;
    }

    unsigned int SphereGPDPulsedGradientCompartment::GetNumberOfParameters()
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

    SphereGPDPulsedGradientCompartment::ModelOutputVectorType &SphereGPDPulsedGradientCompartment::GetCompartmentVector()
    {
        if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
            m_CompartmentVector.SetSize(this->GetCompartmentSize());

        m_CompartmentVector[0] = this->GetTissueRadius();
        m_CompartmentVector[1] = this->GetAxialDiffusivity();

        return m_CompartmentVector;
    }

    double SphereGPDPulsedGradientCompartment::GetApparentFractionalAnisotropy()
    {
        return 0.0;
    }

} // end namespace anima
