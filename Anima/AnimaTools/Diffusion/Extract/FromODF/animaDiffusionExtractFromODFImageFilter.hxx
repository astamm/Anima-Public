#pragma once

#include "animaDiffusionExtractFromODFImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>

namespace anima
{
template <typename TInputPixelType>
void
DiffusionExtractFromODFImageFilter<TInputPixelType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    using InputIteratorType = itk::ImageRegionConstIteratorWithIndex <TInputImage>;
    using OutputIteratorType = itk::ImageRegionIterator <TOutputImage>;

    InputIteratorType inputIt(this->GetInput(),outputRegionForThread);
    OutputIteratorType outIt(this->GetOutput(),outputRegionForThread);

    unsigned int vdim = this->GetInput()->GetNumberOfComponentsPerPixel();
    InputImagePixel tmpCoefs;

    while(!inputIt.IsAtEnd())
    {
        tmpCoefs = inputIt.Get();

        if (isZero(tmpCoefs))
        {
            ++inputIt;
            outIt.Set(0);
            ++outIt;
            continue;
        }

        double sumSquares = 0;
        for (unsigned int i = 0;i < vdim;++i)
            sumSquares += tmpCoefs[i]*tmpCoefs[i];

        outIt.Set(sqrt(1 - tmpCoefs[0]*tmpCoefs[0]/sumSquares));
        ++inputIt;
        ++outIt;
    }
}
	
} // end of namespace anima
