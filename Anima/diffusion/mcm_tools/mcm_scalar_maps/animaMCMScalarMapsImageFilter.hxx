#pragma once
#include "animaMCMScalarMapsImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

namespace anima
{

template <class TPixelType>
void
MCMScalarMapsImageFilter <TPixelType>
::DynamicThreadedGenerateData(const InputRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator <InputImageType> InputIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutputIteratorType;

    InputIteratorType inItr(this->GetInput(), outputRegionForThread);
    unsigned int numOutputs = this->GetNumberOfIndexedOutputs();
    std::vector <OutputIteratorType> outItrs(numOutputs);

    for (unsigned int i = 0;i < numOutputs;++i)
        outItrs[i] = OutputIteratorType(this->GetOutput(i), outputRegionForThread);

    InputImageType *input = const_cast <InputImageType *> (this->GetInput());
    MCModelPointer mcmPtr = input->GetDescriptionModel()->Clone();
    PixelType voxelValue;
    unsigned int numCompartments = mcmPtr->GetNumberOfCompartments();
    unsigned int numIsoCompartments = mcmPtr->GetNumberOfIsotropicCompartments();

    while (!inItr.IsAtEnd())
    {
        voxelValue = inItr.Get();
        if (isZero(voxelValue))
        {
            ++inItr;
            for (unsigned int i = 0;i < numOutputs;++i)
            {
                outItrs[i].Set(0.0);
                ++outItrs[i];
            }

            continue;
        }

        mcmPtr->SetModelVector(voxelValue);
        double anisoWeight = 0.0;
        double faCompartments = 0.0;
        double mdCompartments = 0.0;
        double parDiffCompartments = 0.0;
        double perpDiffCompartments = 0.0;
        double isoRWeight = 0.0;
        double fwWeight = 0.0;

        for (unsigned int i = numIsoCompartments;i < numCompartments;++i)
        {
            double weight = mcmPtr->GetCompartmentWeight(i);
            if (weight <= 0.0)
                continue;
            anima::BaseCompartment *currentComp = mcmPtr->GetCompartment(i);
            anisoWeight += weight;
            faCompartments += weight * currentComp->GetApparentFractionalAnisotropy();
            mdCompartments += weight * currentComp->GetApparentMeanDiffusivity();
            parDiffCompartments += weight * currentComp->GetApparentParallelDiffusivity();
            perpDiffCompartments += weight * currentComp->GetApparentPerpendicularDiffusivity();
        }

        for (unsigned int i = 0;i < numIsoCompartments;++i)
        {
            double weight = mcmPtr->GetCompartmentWeight(i);
            if (weight <= 0.0)
                continue;

            anima::BaseCompartment *currentComp = mcmPtr->GetCompartment(i);

            if (m_IncludeIsotropicWeights)
            {
                mdCompartments += weight * currentComp->GetApparentMeanDiffusivity();
                parDiffCompartments += weight * currentComp->GetApparentParallelDiffusivity();
                perpDiffCompartments += weight * currentComp->GetApparentPerpendicularDiffusivity();
            }

            if (currentComp->GetCompartmentType() == anima::FreeWater)
                fwWeight += weight;
            else if (currentComp->GetCompartmentType() == anima::PlaneSGPPulsedGradient)
                isoRWeight += weight;
        }

        if ((!m_IncludeIsotropicWeights)&&(anisoWeight > 0.0))
        {
            mdCompartments /= anisoWeight;
            faCompartments /= anisoWeight;
            parDiffCompartments /= anisoWeight;
            perpDiffCompartments /= anisoWeight;
        }

        outItrs[0].Set(fwWeight);
        outItrs[1].Set(isoRWeight);
        outItrs[2].Set(anisoWeight);
        outItrs[3].Set(faCompartments);
        outItrs[4].Set(mdCompartments);
        outItrs[5].Set(parDiffCompartments);
        outItrs[6].Set(perpDiffCompartments);

        ++inItr;
        for (unsigned int i = 0;i < numOutputs;++i)
            ++outItrs[i];
    }
}

} // end namespace anima
