#include "animaBaseDistribution.h"

#include <limits>
#include <itkMacro.h>

namespace anima
{

void BaseDistribution::SetShapeParameter(const double val)
{
    if (val < std::numeric_limits<double>::epsilon())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The shape parameter of a statistical distribution should be strictly positive.", ITK_LOCATION);
    m_ShapeParameter = val;
}

void BaseDistribution::SetScaleParameter(const double val)
{
    if (val < std::numeric_limits<double>::epsilon())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The scale parameter of a statistical distribution should be strictly positive.", ITK_LOCATION);
    m_ScaleParameter = val;
}

void BaseDistribution::SetAlphaParameter(const double val)
{
    if (val < std::numeric_limits<double>::epsilon())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The alpha parameter of a statistical distribution should be strictly positive.", ITK_LOCATION);
    m_AlphaParameter = val;
}

void BaseDistribution::SetBetaParameter(const double val)
{
    if (val < std::numeric_limits<double>::epsilon())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The beta parameter of a statistical distribution should be strictly positive.", ITK_LOCATION);
    m_BetaParameter = val;
}
    
} // end of namespace anima
