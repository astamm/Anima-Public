#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMExport.h>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/tools/roots.hpp>

#include <tuple>
#include <map>

// It implements the MR signal attenuation induced by restricted diffusion
// - in spherical geometry
// - under the Gaussian Phase Distribution (GPD) approximation 
// - for pulsed-gradient spin-echo (PGSE) sequences.
//
// [1] Balinov, B., Jonsson, B., Linse, P., & Soderman, O. (1993). The NMR 
// self-diffusion method applied to restricted diffusion. Simulation of 
// echo attenuation from molecules in spheres and between planes. Journal 
// of Magnetic Resonance, Series A, 104(1), 17-25.
//
// [2] Murday, J. S., & Cotts, R. M. (1968). Self‐diffusion coefficient of 
// liquid lithium. The Journal of Chemical Physics, 48(11), 4938-4945.

namespace anima
{

    class SphereGPDPulsedGradientFunctor
    {
    public:
        SphereGPDPulsedGradientFunctor()
        {
            m_BesselOrder = 1.0;
        }

        void SetBesselOrder(double val)
        {
            m_BesselOrder = val;
        }

        double GetBesselOrder()
        {
            return m_BesselOrder;
        }

        double operator()(double const &x)
        {
            return x * boost::math::cyl_bessel_j_prime(m_BesselOrder, x) - boost::math::cyl_bessel_j(m_BesselOrder, x) / 2.0;
        }

    private:
        double m_BesselOrder;
    };

    class ANIMAMCM_EXPORT SphereGPDPulsedGradientCompartment : public BaseCompartment
    {
    public:
        // Useful typedefs
        using Self = SphereGPDPulsedGradientCompartment;
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
        itkTypeMacro(SphereGPDPulsedGradientCompartment, BaseCompartment);

        DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE
        {
            return SphereGPDPulsedGradient;
        }

        virtual double GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
        virtual ListType &GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
        virtual double GetLogDiffusionProfile(const Vector3DType &sample) ITK_OVERRIDE;

        virtual void SetParametersFromVector(const ListType &params) ITK_OVERRIDE;
        virtual ListType &GetParametersAsVector() ITK_OVERRIDE;

        virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
        virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

        // Set constraints
        void SetEstimateAxialDiffusivity(bool arg);
        void SetEstimateTissueRadius(bool arg);
        void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

        unsigned int GetCompartmentSize() ITK_OVERRIDE;
        unsigned int GetNumberOfParameters() ITK_OVERRIDE;
        ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

        void SetTissueRadius(double num) ITK_OVERRIDE;
        void SetAxialDiffusivity(double num) ITK_OVERRIDE;

        bool GetTensorCompatible() ITK_OVERRIDE { return false; }
        double GetApparentFractionalAnisotropy() ITK_OVERRIDE;

    protected:
        SphereGPDPulsedGradientCompartment() : Superclass()
        {
            m_EstimateAxialDiffusivity = true;
            m_EstimateTissueRadius = true;
            m_ChangedConstraints = true;

            m_SignalSummationTolerance = 1.0e-4;

            double besselOrder = 1.5;
            m_FunctorZeros.resize(m_MaximumNumberOfSumElements);
            SphereGPDPulsedGradientFunctor functorObj;
            functorObj.SetBesselOrder(besselOrder);

            const std::uintmax_t maxit = 20;                  // Limit to maximum iterations.
            int digits = std::numeric_limits<double>::digits; // Maximum possible binary digits accuracy for type T.
            // Some fraction of digits is used to control how accurate to try to make the result.
            int get_digits = digits - 3;                               // We have to have a non-zero interval at each step, so
                                                                       // maximum accuracy is digits - 1.  But we also have to
                                                                       // allow for inaccuracy in f(x), otherwise the last few
                                                                       // iterations just thrash around.
            boost::math::tools::eps_tolerance<double> tol(get_digits); // Set the tolerance.

            std::pair<double, double> r;
            double lowerBound = 1.5;
            double upperBound = lowerBound + M_PI;
            for (unsigned int i = 0; i < m_MaximumNumberOfSumElements; ++i)
            {
                std::uintmax_t it = maxit; // Initially our chosen max iterations, but updated with actual.
                r = boost::math::tools::toms748_solve(functorObj, lowerBound, upperBound, tol, it);
                m_FunctorZeros[i] = r.first + (r.second - r.first) / 2;
                lowerBound += M_PI;
                upperBound += M_PI;
            }
        }

        virtual ~SphereGPDPulsedGradientCompartment() {}

        using KeyType = std::tuple<unsigned int, unsigned int, unsigned int>;
        using MapType = std::map<KeyType, double>;

        KeyType GenerateKey(double smallDelta, double bigDelta, double gradientStrength);

        void UpdateSignals(double smallDelta, double bigDelta, double gradientStrength);

    private:
        bool m_EstimateAxialDiffusivity, m_EstimateTissueRadius;
        bool m_ChangedConstraints;
        unsigned int m_NumberOfParameters;

        MapType m_FirstSummations, m_SecondSummations;
        MapType m_ThirdSummations, m_FourthSummations;

        double m_SignalSummationTolerance;

        std::vector<double> m_FunctorZeros;
        const unsigned int m_MaximumNumberOfSumElements = 100;
    };

} // end namespace anima
