#include "animaBetaDistribution.h"
#include <animaGammaFunctions.h>

#include <cmath>
#include <limits>

#include <itkMacro.h>

#include <boost/math/special_functions/beta.hpp>

namespace anima
{

bool BetaDistribution::BelongsToSupport(const ValueType &x)
{
    return x >= 0.0 && x <= 1.0;
}

void BetaDistribution::SetAlphaParameter(const double val)
{
    if (val < this->GetEpsilon())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The alpha parameter of the Beta distribution should be strictly positive.", ITK_LOCATION);
    m_AlphaParameter = val;
}

void BetaDistribution::SetBetaParameter(const double val)
{
    if (val < this->GetEpsilon())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The beta parameter of the Beta distribution should be strictly positive.", ITK_LOCATION);
    m_BetaParameter = val;
}

double BetaDistribution::GetDensity(const ValueType &x)
{
    if (!this->BelongsToSupport(x))
        return 0.0;
    return std::exp(this->GetLogDensity(x));
}

double BetaDistribution::GetLogDensity(const ValueType &x)
{
    if (!this->BelongsToSupport(x))
        throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density of the Beta distribution is not defined outside [0,1].", ITK_LOCATION);
    double alphaParameter = this->GetAlphaParameter();
    double betaParameter = this->GetBetaParameter();
    double logBetaFunctionValue = std::lgamma(alphaParameter) + std::lgamma(betaParameter) - std::lgamma(alphaParameter + betaParameter);
    double resValue = (alphaParameter - 1.0) * std::log(x);
    resValue += (betaParameter - 1.0) * std::log(1.0 - x);
    resValue -= logBetaFunctionValue;
    return resValue;
}

double BetaDistribution::GetDensityDerivative(const ValueType &x)
{
    if (!this->BelongsToSupport(x))
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

double BetaDistribution::GetCumulative(const ValueType &x)
{
    /**
    * \fn double GetCumulative(const double &x)
    *
    * \author Aymeric Stamm (2024).
    *
    * \param x A numeric value specifying a point on the support of the Beta distribution.
    *
    * \return A numeric value storing the value of the cumulative distribution function of
    * the Beta distribution at point `x`. See: https://statproofbook.github.io/P/beta-cdf.
    */
    
    if (x <= 0.0)
        return 0.0;

    if (x >= 1.0)
        return 1.0;

    return boost::math::ibeta<double, double, double>(this->GetAlphaParameter(), this->GetBetaParameter(), x);
}

void BetaDistribution::Fit(const SampleType &sample, const std::string &method)
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

void BetaDistribution::Random(SampleType &sample, GeneratorType &generator)
{
    BetaDistributionType betaDist(this->GetAlphaParameter(), this->GetBetaParameter());
    UniformDistributionType unifDist(0.0, 1.0);

    unsigned int numObservations = sample.size();
    for (unsigned int i = 0;i < numObservations;++i)
        sample[i] = boost::math::quantile(betaDist, unifDist(generator));
}

double BetaDistribution::GetDistance(Self *otherDistribution)
{
    /**
    * \fn double GetDistance(BetaDistribution *otherDistribution)
    *
    * \author Aymeric Stamm (2024).
    *
    * \param otherDistribution A pointer specifying another object of class `BetaDistribution`.
    *
    * \return A numeric value storing the symmetric Kullback-Leibler divergence between the
    * current Beta distribution and the input Beta distribution. See 
    * https://math.stackexchange.com/questions/257821/kullback-liebler-divergence#comment564291_257821
    * and Rauber et al., Probabilistic distance measures of the Dirichlet and Beta distributions, 2008.
    */

    double thisAlpha = this->GetAlphaParameter();
    double thisBeta = this->GetBetaParameter();

    BetaDistribution *betaDistr = dynamic_cast<BetaDistribution *>(otherDistribution);
    double otherAlpha = betaDistr->GetAlphaParameter();
    double otherBeta = betaDistr->GetBetaParameter();

    double thisAlphaDigamma = anima::digamma(thisAlpha);
    double otherAlphaDigamma = anima::digamma(otherAlpha);
    double thisBetaDigamma = anima::digamma(thisBeta);
    double otherBetaDigamma = anima::digamma(otherBeta);
    double diffSumDigamma = anima::digamma(otherAlpha + otherBeta) - anima::digamma(thisAlpha + thisBeta);
    
    double distanceValue = (thisAlpha - otherAlpha) * (thisAlphaDigamma - otherAlphaDigamma + diffSumDigamma);
    distanceValue += (thisBeta - otherBeta) * (thisBetaDigamma - otherBetaDigamma + diffSumDigamma);
    return distanceValue;
}
    
} // end of namespace anima
