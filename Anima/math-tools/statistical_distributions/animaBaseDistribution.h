#pragma once

#include <string>
#include <vector>
#include <random>

#include <AnimaStatisticalDistributionsExport.h>

namespace anima
{
	class ANIMASTATISTICALDISTRIBUTIONS_EXPORT BaseDistribution
	{
	public:
		using GeneratorType = std::mt19937;
		using VectorType = std::vector<double>;

		BaseDistribution()
		{
			m_ShapeParameter = 1.0;
			m_ScaleParameter = 1.0;
			m_AlphaParameter = 1.0;
			m_BetaParameter = 1.0;
		}

		void SetShapeParameter(const double val);
		double GetShapeParameter() {return m_ShapeParameter;}

		void SetScaleParameter(const double val);
		double GetScaleParameter() {return m_ScaleParameter;}

		void SetAlphaParameter(const double val);
		double GetAlphaParameter() {return m_AlphaParameter;}

		void SetBetaParameter(const double val);
		double GetBetaParameter() {return m_BetaParameter;}

		virtual double GetDensity(const double &x) = 0;
		virtual double GetLogDensity(const double &x) = 0;
		virtual double GetDensityDerivative(const double &x) {return 0.0;}
		virtual void Fit(const VectorType &sample, const std::string &method) = 0;
		virtual void Random(VectorType &sample, GeneratorType &generator) = 0;

	private:
		double m_ShapeParameter;
		double m_ScaleParameter;
		double m_AlphaParameter;
		double m_BetaParameter;
	};
} // end of namespace
