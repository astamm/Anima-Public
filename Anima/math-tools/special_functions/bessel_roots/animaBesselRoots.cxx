#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/tools/roots.hpp>

#include <tclap/CmdLine.h>

struct BesselPrimeFunctor
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

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg<unsigned int> nrootsArg(
        "n", "nb-roots",
        "An integer value specifying the number of roots to compute (default: 1000).",
        false, 1000, "number of roots", cmd);

    TCLAP::ValueArg<double> orderArg(
        "o", "order",
        "A double value specifying the order of the Bessel function (default: 1.0).",
        false, 1.0, "order of Bessel function", cmd);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    unsigned int numberOfRoots = nrootsArg.getValue();
    double order = orderArg.getValue();

    std::vector<double> roots;
    boost::math::cyl_bessel_j_zero(order, 0, numberOfRoots, std::back_inserter(roots));

    std::cout << "--- ROOTS ---" << std::endl;
    std::copy(roots.begin(), roots.end(), std::ostream_iterator<double>(std::cout, "\n"));

    BesselPrimeFunctor besselPrime;
    besselPrime.SetBesselOrder(order);

    const std::uintmax_t maxit = 20;                  // Limit to maximum iterations.
    int digits = std::numeric_limits<double>::digits; // Maximum possible binary digits accuracy for type T.
    // Some fraction of digits is used to control how accurate to try to make the result.
    int get_digits = digits - 3;                               // We have to have a non-zero interval at each step, so
                                                               // maximum accuracy is digits - 1.  But we also have to
                                                               // allow for inaccuracy in f(x), otherwise the last few
                                                               // iterations just thrash around.
    boost::math::tools::eps_tolerance<double> tol(get_digits); // Set the tolerance.

    std::pair<double, double> r;
    std::vector<double> extremas(numberOfRoots - 1);
    std::cout << "--- EXTREMAS ---" << std::endl;
    for (unsigned int i = 0; i < numberOfRoots - 1; ++i)
    {
        std::uintmax_t it = maxit; // Initially our chosen max iterations, but updated with actual.
        r = boost::math::tools::toms748_solve(besselPrime, roots[i], roots[i + 1], tol, it);
        extremas[i] = r.first + (r.second - r.first) / 2;

        std::cout << (roots[i] + roots[i + 1]) / 2 << " " << extremas[i] << std::endl;
    }

    return EXIT_SUCCESS;
}
