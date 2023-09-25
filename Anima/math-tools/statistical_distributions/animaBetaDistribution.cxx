#include "animaBetaDistribution.h"
#include <cmath>
#include <limits>
#include <itkMacro.h>

namespace anima
{

double BetaDistribution::GetDensity(const double &x)
{
    if (x < std::numeric_limits<double>::epsilon() || 1.0 - x < std::numeric_limits<double>::epsilon())
        return 0.0;
    return std::exp(this->GetLogDensity(x));
}

double BetaDistribution::GetLogDensity(const double &x)
{
    if (x < std::numeric_limits<double>::epsilon() || 1.0 - x < std::numeric_limits<double>::epsilon())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density of the Beta distribution is not defined for values outside (0,1).", ITK_LOCATION);
    double alphaParameter = this->GetAlphaParameter();
    double betaParameter = this->GetBetaParameter();
    double logBetaFunctionValue = std::lgamma(alphaParameter) + std::lgamma(betaParameter) - std::lgamma(alphaParameter + betaParameter);
    double resValue = (alphaParameter - 1.0) * std::log(x);
    resValue += (betaParameter - 1.0) * std::log(1.0 - x);
    resValue -= logBetaFunctionValue;
    return resValue;
}

double BetaDistribution::GetDensityDerivative(const double &x)
{
    if (x < std::numeric_limits<double>::epsilon() || 1.0 - x < std::numeric_limits<double>::epsilon())
        return 0.0;
    double alphaParameter = this->GetAlphaParameter();
    double betaParameter = this->GetBetaParameter();
    double logBetaFunctionValue = std::lgamma(alphaParameter) + std::lgamma(betaParameter) - std::lgamma(alphaParameter + betaParameter);
    double resValue = (alphaParameter - 2.0) * std::log(x);
    resValue += (betaParameter - 2.0) * std::log(1.0 - x);
    resValue -= logBetaFunctionValue;
    resValue = std::exp(resValue);
    resValue *= ((alphaParameter - 1.0) * (1.0 - x) + (betaParameter - 1.0) * x);
    return resValue;
}

void BetaDistribution::Fit(const VectorType &sample, const std::string &method)
{
    // Method-of-moments: https://en.wikipedia.org/wiki/Beta_distribution#Two_unknown_parameters
    double alphaParameter = 1.0;
    double betaParameter = 1.0;

    unsigned int numObservations = sample.size();
    
    double sampleMean = 0.0;
    for (unsigned int i = 0;i < numObservations;++i)
        sampleMean += sample[i];
    sampleMean /= numObservations;

    double sampleVariance = 0.0;
    for (unsigned int i = 0;i < numObservations;++i)
        sampleVariance += (sample[i] - sampleMean) * (sample[i] - sampleMean);
    sampleVariance /= (numObservations - 1.0);

    if (sampleVariance > std::numeric_limits<double>::epsilon() && sampleVariance < sampleMean * (1.0 - sampleMean))
    {
        double tmpValue = (sampleMean * (1.0 - sampleMean) / sampleVariance - 1.0);
        alphaParameter = sampleMean * tmpValue;
        betaParameter = (1.0 - sampleMean) * tmpValue;
    }

    this->SetAlphaParameter(alphaParameter);
    this->SetBetaParameter(betaParameter);
}

void BetaDistribution::Random(VectorType &sample, GeneratorType &generator)
{
    BetaDistributionType betaDist(this->GetAlphaParameter(), this->GetBetaParameter());
    UniformDistributionType unifDist(0.0, 1.0);

    unsigned int numObservations = sample.size();
    for (unsigned int i = 0;i < numObservations;++i)
        sample[i] = boost::math::quantile(betaDist, unifDist(generator));
}
    
} // end of namespace anima
