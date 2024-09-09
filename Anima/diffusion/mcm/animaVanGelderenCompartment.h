#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMExport.h>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/tools/roots.hpp>

#include <tuple>
#include <map>

namespace anima
{

    class BesselPrimeFunctor
    {
    public:
        BesselPrimeFunctor()
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
            return boost::math::cyl_bessel_j_prime(m_BesselOrder, x);
        }

    private:
        double m_BesselOrder;
    };

    class ANIMAMCM_EXPORT VanGelderenCompartment : public BaseCompartment
    {
    public:
        // Useful typedefs
        using Self = VanGelderenCompartment;
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
        itkTypeMacro(VanGelderenCompartment, BaseCompartment);

        DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE
        {
            return VanGelderen;
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
        VanGelderenCompartment() : Superclass()
        {
            m_EstimateAxialDiffusivity = true;
            m_EstimateTissueRadius = true;
            m_ChangedConstraints = true;

            m_SignalSummationTolerance = 1.0e-4;

            // Compute and store the zeros of first-order Bessel prime function
            double besselOrder = 1;
            std::vector<double> roots;
            boost::math::cyl_bessel_j_zero(besselOrder, 0, m_MaximumNumberOfSumElements + 1, std::back_inserter(roots));
            m_BesselPrimeZeros.resize(m_MaximumNumberOfSumElements);

            BesselPrimeFunctor besselPrime;
            besselPrime.SetBesselOrder(besselOrder);

            const std::uintmax_t maxit = 20;                  // Limit to maximum iterations.
            int digits = std::numeric_limits<double>::digits; // Maximum possible binary digits accuracy for type T.
            // Some fraction of digits is used to control how accurate to try to make the result.
            int get_digits = digits - 3;                               // We have to have a non-zero interval at each step, so
                                                                       // maximum accuracy is digits - 1.  But we also have to
                                                                       // allow for inaccuracy in f(x), otherwise the last few
                                                                       // iterations just thrash around.
            boost::math::tools::eps_tolerance<double> tol(get_digits); // Set the tolerance.

            std::pair<double, double> r;
            for (unsigned int i = 0; i < m_MaximumNumberOfSumElements; ++i)
            {
                std::uintmax_t it = maxit; // Initially our chosen max iterations, but updated with actual.
                r = boost::math::tools::toms748_solve(besselPrime, roots[i], roots[i + 1], tol, it);
                m_BesselPrimeZeros[i] = r.first + (r.second - r.first) / 2;
            }
        }

        virtual ~VanGelderenCompartment() {}

        typedef std::tuple<unsigned int, unsigned int, unsigned int> KeyType;
        typedef std::map<KeyType, double> MapType;

        KeyType GenerateKey(double smallDelta, double bigDelta, double gradientStrength);

        void UpdateSignals(double smallDelta, double bigDelta, double gradientStrength);

    private:
        bool m_EstimateAxialDiffusivity, m_EstimateTissueRadius;
        bool m_ChangedConstraints;
        unsigned int m_NumberOfParameters;

        MapType m_FirstSummations, m_SecondSummations;
        MapType m_ThirdSummations, m_FourthSummations;

        double m_SignalSummationTolerance;

        std::vector<double> m_BesselPrimeZeros;
        const unsigned int m_MaximumNumberOfSumElements = 2000;
    };

} // end namespace anima
