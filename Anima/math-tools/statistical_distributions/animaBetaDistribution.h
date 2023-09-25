#pragma once

#include <animaBaseDistribution.h>

#include <boost/math/distributions/beta.hpp>

namespace anima
{
    class ANIMASTATISTICALDISTRIBUTIONS_EXPORT BetaDistribution : public BaseDistribution
	{
	public:
		using UniformDistributionType = std::uniform_real_distribution<double>;
		using BetaDistributionType = boost::math::beta_distribution<double>;
		
		double GetDensity(const double &x);
		double GetLogDensity(const double &x);
		double GetDensityDerivative(const double &x);
		void Fit(const VectorType &sample, const std::string &method);
		void Random(VectorType &sample, GeneratorType &generator);
	};

} // end of namespace anima
