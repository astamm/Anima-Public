#pragma once

#include <animaBaseDistribution.h>

namespace anima
{
	class ANIMASTATISTICALDISTRIBUTIONS_EXPORT GammaDistribution : public BaseDistribution
	{
	public:
		using DistributionType = std::gamma_distribution<>;
		double GetDensity(const double &x);
		double GetLogDensity(const double &x);
		double GetDensityDerivative(const double &x);
		void Fit(const VectorType &sample, const std::string &method);
		void Random(VectorType &sample, GeneratorType &generator);
	};
    
} // end of namespace
