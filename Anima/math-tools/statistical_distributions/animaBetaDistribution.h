#pragma once

#include <animaBaseDistribution.h>

#include <boost/math/distributions/beta.hpp>

namespace anima
{
    class ANIMASTATISTICALDISTRIBUTIONS_EXPORT BetaDistribution : public BaseDistribution<double>
	{
	public:
		using UniformDistributionType = std::uniform_real_distribution<double>;
		using BetaDistributionType = boost::math::beta_distribution<double>;

		BetaDistribution()
		{
			m_AlphaParameter = 1.0;
			m_BetaParameter = 1.0;
		}
		
		bool BelongsToSupport(const ValueType &x);
		double GetDensity(const ValueType &x);
		double GetLogDensity(const ValueType &x);
		double GetDensityDerivative(const ValueType &x);
		double GetCumulative(const ValueType &x);
		void Fit(const SampleType &sample, const std::string &method);
		void Random(SampleType &sample, GeneratorType &generator);
		ValueType GetMean() { return m_AlphaParameter / (m_AlphaParameter + m_BetaParameter); }
		double GetVariance() { return m_AlphaParameter * m_BetaParameter / ((m_AlphaParameter + m_BetaParameter) * (m_AlphaParameter + m_BetaParameter) * (m_AlphaParameter + m_BetaParameter + 1.0)); }
		double GetDistance(Self *otherDistribution);

		void SetAlphaParameter(const double val);
		double GetAlphaParameter() { return m_AlphaParameter; }

		void SetBetaParameter(const double val);
		double GetBetaParameter() { return m_BetaParameter; }

	private:
		double m_AlphaParameter;
		double m_BetaParameter;
	};

} // end of namespace anima
