#pragma once

#include <itkSingleValuedCostFunction.h>

#include <vnl/vnl_matrix.h>

#include <vector>

namespace anima
{

class DiffusionEstimateTensorsUsingNonCentralChiNoiseCostFunction :
        public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    using Self = DiffusionEstimateTensorsUsingNonCentralChiNoiseCostFunction;
    using Superclass = SingleValuedCostFunction;
    using Pointer = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer<const Self>;

    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(DiffusionEstimateTensorsUsingNonCentralChiNoiseCostFunction, SingleValuedCostFunction);

    using MeasureType = Superclass::MeasureType;
    using DerivativeType = Superclass::DerivativeType;
    using ParametersType = Superclass::ParametersType;

    virtual MeasureType GetValue(const ParametersType & parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE;

    // ItkSetMacro does not handle vector in debug
    //itkSetMacro(RawSignal, std::vector <double>);
    void SetRawSignal(std::vector <double> rawSignal ){m_RawSignal = rawSignal;}

    itkSetMacro(DesignMatrix, vnl_matrix<double>);
    itkSetMacro(Sigma, double);
    itkSetMacro(NumberOfParameters, unsigned int);
    itkSetMacro(NumberOfCoils, unsigned int);

    itkSetMacro(B0Value, double);
    itkSetMacro(UseB0Value, bool);

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE
    {
        return m_NumberOfParameters;
    }

protected:
    DiffusionEstimateTensorsUsingNonCentralChiNoiseCostFunction();

    virtual ~DiffusionEstimateTensorsUsingNonCentralChiNoiseCostFunction()
    {
    }

private:
    DiffusionEstimateTensorsUsingNonCentralChiNoiseCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector <double> m_RawSignal;
    unsigned int m_NumberOfParameters;

    vnl_matrix <double> m_DesignMatrix;
    double m_Sigma;

    unsigned int m_NumberOfCoils;

    bool m_UseB0Value;
    double m_B0Value;
};

}
