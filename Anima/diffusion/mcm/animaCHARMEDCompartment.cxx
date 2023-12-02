#include "animaCHARMEDCompartment.h"
#include <animaMCMConstants.h>
#include <animaVectorOperations.h>

#include <limits>

namespace anima
{

    void CHARMEDCompartment::PrepareData(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
    {
        m_BValue = anima::GetBValueFromAcquisitionParameters(smallDelta, bigDelta, gradientStrength);
        Vector3DType sphCoords, carCoords;
        sphCoords[0] = this->GetOrientationTheta();
        sphCoords[1] = this->GetOrientationPhi();
        sphCoords[2] = 1.0;
        anima::TransformSphericalToCartesianCoordinates(sphCoords, carCoords);
        m_InnerProd = anima::ComputeScalarProduct(gradient, carCoords);
        double axialDiff = this->GetAxialDiffusivity();
        double radialDiff = this->GetRadialDiffusivity1();
        m_OrthogonalGradient = gradient - m_InnerProd * carCoords;
        double orthogonalGradientNorm = m_OrthogonalGradient.magnitude();
        m_OrthogonalGradientStrength = gradientStrength * orthogonalGradientNorm;
        m_OrthogonalGradient.normalize();

        // Axial contribution
        m_AxialSignal = std::exp(-m_BValue * axialDiff * m_InnerProd * m_InnerProd);

        // Hindered compartment, radial contribution
        m_RadialHinderedSignal = std::exp(-m_BValue * radialDiff * orthogonalGradientNorm * orthogonalGradientNorm);

        // Restricted compartment, radial contribution
        m_RadialRestrictedSignal = m_StaniszCompartment->GetFourierTransformedDiffusionProfile(smallDelta, bigDelta, m_OrthogonalGradientStrength, m_OrthogonalGradient);
    }

    double CHARMEDCompartment::GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
    {
        this->PrepareData(smallDelta, bigDelta, gradientStrength, gradient);

        // CHARMED signal
        double extraAxonalFraction = this->GetExtraAxonalFraction();
        double signalValue = extraAxonalFraction * m_RadialHinderedSignal;
        signalValue += (1.0 - extraAxonalFraction) * m_RadialRestrictedSignal;
        signalValue *= m_AxialSignal;

        return signalValue;
    }

    CHARMEDCompartment::ListType &CHARMEDCompartment::GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
    {
        m_JacobianVector.resize(this->GetNumberOfParameters());

        this->PrepareData(smallDelta, bigDelta, gradientStrength, gradient);

        // params are
        // 1) m_OrientationTheta
        // 2) m_OrientationPhi
        // 3) m_ExtraAxonalFraction
        // 4) m_TissueRadius --> choose to fix it
        // 5) m_AxialDiffusivity - m_RadialDiffusivity1 --> choose to fix it
        // 6) m_RadialDiffusivity1 --> choose to fix it

        double signalValue = this->GetFourierTransformedDiffusionProfile(smallDelta, bigDelta, gradientStrength, gradient);
        ListType staniszJacobian = m_StaniszCompartment->GetSignalAttenuationJacobian(smallDelta, bigDelta, m_OrthogonalGradientStrength, m_OrthogonalGradient);

        double angleCommonValue = 2.0 * this->GetExtraAxonalFraction() * m_BValue * this->GetRadialDiffusivity1() * m_RadialHinderedSignal;
        if (std::abs(1.0 - m_InnerProd * m_InnerProd) > std::sqrt(std::numeric_limits<double>::epsilon()))
            angleCommonValue -= (1.0 - this->GetExtraAxonalFraction()) * this->GetTissueRadius() * staniszJacobian[0] / (1.0 - m_InnerProd * m_InnerProd);
        angleCommonValue *= m_AxialSignal;
        angleCommonValue -= 2.0 * m_BValue * this->GetAxialDiffusivity() * signalValue;
        angleCommonValue *= m_InnerProd;

        double extraAxonalFraction = this->GetExtraAxonalFraction();

        unsigned int pos = 0;

        // Derivative w.r.t. m_OrientationTheta
        double innerProdThetaDerivative = gradient[0] * std::cos(this->GetOrientationPhi());
        innerProdThetaDerivative += gradient[1] * std::sin(this->GetOrientationPhi());
        innerProdThetaDerivative *= std::cos(this->GetOrientationTheta());
        innerProdThetaDerivative -= gradient[2] * std::sin(this->GetOrientationTheta());
        m_JacobianVector[pos] = angleCommonValue * innerProdThetaDerivative;
        ++pos;

        // Derivative w.r.t. m_OrientationPhi
        double innerProdPhiDerivative = gradient[1] * std::cos(this->GetOrientationPhi());
        innerProdPhiDerivative -= gradient[0] * std::sin(this->GetOrientationPhi());
        innerProdPhiDerivative *= std::sin(this->GetOrientationTheta());
        m_JacobianVector[pos] += angleCommonValue * innerProdPhiDerivative;
        ++pos;

        // Derivative w.r.t. m_ExtraAxonalFraction
        m_JacobianVector[pos] = m_AxialSignal * (m_RadialHinderedSignal - m_RadialRestrictedSignal);
        ++pos;

        // Derivative w.r.t. m_TissueRadius
        m_JacobianVector[pos] = (1.0 - extraAxonalFraction) * m_AxialSignal * staniszJacobian[0];
        ++pos;

        if (m_EstimateDiffusivities)
        {
            double perpDeriv = -m_BValue * (1.0 - m_InnerProd * m_InnerProd) * m_AxialSignal * extraAxonalFraction * m_RadialHinderedSignal;
            double paraDeriv = -m_BValue * m_InnerProd * m_InnerProd * signalValue;
            paraDeriv += m_AxialSignal * (1.0 - extraAxonalFraction) * staniszJacobian[1];

            // Derivative w.r.t. m_AxialDiffusivity - m_RadialDiffusivity1
            m_JacobianVector[pos] = paraDeriv - perpDeriv;
            ++pos;

            // Derivative w.r.t. m_RadialDiffusivity1
            m_JacobianVector[pos] += perpDeriv;
        }

        return m_JacobianVector;
    }

    double CHARMEDCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
    {
        // Equivalent tensor for the hindered compartment
        Matrix3DType hinderedTensor;
        hinderedTensor.Fill(0.0);
        Vector3DType sphCoords, carCoords;
        sphCoords[0] = this->GetOrientationTheta();
        sphCoords[1] = this->GetOrientationPhi();
        sphCoords[2] = 1.0;
        anima::TransformSphericalToCartesianCoordinates(sphCoords, carCoords);

        for (unsigned int i = 0; i < 3; ++i)
        {
            hinderedTensor(i, i) = this->GetRadialDiffusivity1();
            for (unsigned int j = i; j < 3; ++j)
            {
                hinderedTensor(i, j) += (this->GetAxialDiffusivity() - this->GetRadialDiffusivity1()) * carCoords[i] * carCoords[j];
                if (i != j)
                    hinderedTensor(j, i) = hinderedTensor(i, j);
            }
        }

        hinderedTensor = hinderedTensor.GetInverse();
        double hinderedDet = this->GetAxialDiffusivity() * this->GetRadialDiffusivity1() * this->GetRadialDiffusivity1();

        double hinderedQuad = 0.0;
        for (unsigned int i = 0; i < 3; ++i)
        {
            for (unsigned int j = i; j < 3; ++j)
            {
                hinderedQuad += sample[i] * sample[j] * hinderedTensor(i, j);
                if (i != j)
                    hinderedQuad += sample[j] * sample[i] * hinderedTensor(j, i);
            }
        }

        double hinderedDensity = std::exp(-0.5 * hinderedQuad) / std::sqrt(hinderedDet * std::pow(2.0 * M_PI, 3.0));

        // Equivalent tensor for the restricted compartment
        Matrix3DType restrictedTensor;
        restrictedTensor.Fill(0.0);
        double axialDiff = this->GetAxialDiffusivity();
        Vector3DType orthogonalVector;
        orthogonalVector[0] = -carCoords[1];
        orthogonalVector[1] = carCoords[0];
        orthogonalVector[2] = 0.0;
        double bValue = 1000.0;
        double gradientStrength = anima::GetGradientStrengthFromBValue(bValue, anima::DiffusionSmallDelta, anima::DiffusionBigDelta);
        double signalValue = m_StaniszCompartment->GetFourierTransformedDiffusionProfile(anima::DiffusionSmallDelta, anima::DiffusionBigDelta, gradientStrength, orthogonalVector);
        double radialDiff = -std::log(signalValue) / bValue;

        for (unsigned int i = 0; i < 3; ++i)
        {
            restrictedTensor(i, i) = radialDiff;
            for (unsigned int j = i; j < 3; ++j)
            {
                restrictedTensor(i, j) += (axialDiff - radialDiff) * carCoords[i] * carCoords[j];
                if (i != j)
                    restrictedTensor(j, i) = restrictedTensor(i, j);
            }
        }

        restrictedTensor = restrictedTensor.GetInverse();
        double restrictedDet = axialDiff * radialDiff * radialDiff;

        double restrictedQuad = 0.0;
        for (unsigned int i = 0; i < 3; ++i)
        {
            for (unsigned int j = i; j < 3; ++j)
            {
                restrictedQuad += sample[i] * sample[j] * restrictedTensor(i, j);
                if (i != j)
                    restrictedQuad += sample[j] * sample[i] * restrictedTensor(j, i);
            }
        }

        double restrictedDensity = std::exp(-0.5 * restrictedQuad) / std::sqrt(restrictedDet * std::pow(2.0 * M_PI, 3.0));

        double densityValue = this->GetExtraAxonalFraction() * hinderedDensity;
        densityValue += (1.0 - this->GetExtraAxonalFraction()) * restrictedDensity;

        return densityValue;
    }

    void CHARMEDCompartment::SetTissueRadius(double num)
    {
        if (num != this->GetTissueRadius())
        {
            m_StaniszCompartment->SetTissueRadius(num);
            this->Superclass::SetTissueRadius(num);
        }
    }

    void CHARMEDCompartment::SetAxialDiffusivity(double num)
    {
        if (num != this->GetAxialDiffusivity())
        {
            m_StaniszCompartment->SetAxialDiffusivity(num);
            this->Superclass::SetAxialDiffusivity(num);
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

        this->SetTissueRadius(params[pos]);
        ++pos;

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

        m_ParametersVector[pos] = this->GetTissueRadius();
        ++pos;

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

        m_ParametersUpperBoundsVector[pos] = anima::MCMZeroLowerBound;
        ++pos;

        m_ParametersUpperBoundsVector[pos] = anima::MCMZeroLowerBound;
        ++pos;

        m_ParametersUpperBoundsVector[pos] = anima::MCMZeroLowerBound;
        ++pos;

        m_ParametersUpperBoundsVector[pos] = anima::MCMTissueRadiusLowerBound;
        ++pos;

        if (m_EstimateDiffusivities)
        {
            m_ParametersUpperBoundsVector[pos] = anima::MCMZeroLowerBound;
            ++pos;

            m_ParametersUpperBoundsVector[pos] = anima::MCMDiffusivityLowerBound;
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

        m_ParametersUpperBoundsVector[pos] = anima::MCMTissueRadiusUpperBound;
        ++pos;

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

        m_EstimateDiffusivities = arg;
        m_StaniszCompartment->SetEstimateAxialDiffusivity(arg);
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
        this->SetAxialDiffusivity(compartmentVector[4] + compartmentVector[5]);
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
        m_CompartmentVector[4] = this->GetAxialDiffusivity() - this->GetRadialDiffusivity1();
        m_CompartmentVector[5] = this->GetRadialDiffusivity1();

        return m_CompartmentVector;
    }

    double CHARMEDCompartment::GetApparentFractionalAnisotropy()
    {
        return 0.0;
    }

} // end namespace anima
