#pragma once
#include "animaMatrixOperations.h"

#include <animaVectorOperations.h>
#include <animaLinearTransformEstimationTools.h>
#include <itkSymmetricEigenAnalysis.h>

namespace anima
{

template <class VectorType>
itk::Matrix <double,3,3>
GetRotationMatrixFromVectors(const VectorType &first_direction, const VectorType &second_direction, const unsigned int dimension)
{
    // vnl_matrix<double> tmpMatrix(dimension+1,dimension+1,0);
    // anima::pairingToQuaternion(first_direction,second_direction,tmpMatrix);
    // vnl_matrix<double> AMatrix = tmpMatrix.transpose()*tmpMatrix;

    // // Needs two points, providing one on normal vector (cross product)
    // VectorType tmpVec;
    // tmpVec[0] = first_direction[1] * second_direction[2] - first_direction[2] * second_direction[1];
    // tmpVec[1] = first_direction[2] * second_direction[0] - first_direction[0] * second_direction[2];
    // tmpVec[2] = first_direction[0] * second_direction[1] - first_direction[1] * second_direction[0];

    // double normTmpVec = 0;
    // for (unsigned int i = 0;i < dimension;++i)
    //     normTmpVec += tmpVec[i] * tmpVec[i];

    // normTmpVec = std::sqrt(normTmpVec);
    // if (normTmpVec < 1.0e-8)
    // {
    //     tmpVec[0] = 0;
    //     tmpVec[1] = 0;
    //     tmpVec[2] = 1;
    // }

    // anima::pairingToQuaternion(tmpVec,tmpVec,tmpMatrix);
    // AMatrix += tmpMatrix.transpose()*tmpMatrix;

    // itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > eigenSystem(dimension+1);
    // vnl_matrix <double> eVec(dimension+1,dimension+1);
    // vnl_diag_matrix <double> eVals(dimension+1);

    // eigenSystem.SetOrderEigenValues(true);
    // eigenSystem.ComputeEigenValuesAndVectors(AMatrix, eVals, eVec);

    // vnl_matrix <double> rotationMatrix = anima::computeRotationFromQuaternion<double,double>(eVec.get_row(0));


    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    VectorType aVector, bVector, vVector;
    anima::Normalize(first_direction, aVector);
    anima::Normalize(second_direction, bVector);
    anima::ComputeCrossProduct(aVector, bVector, vVector);
    double cValue = anima::ComputeScalarProduct(aVector, bVector);
    if (cValue < -1.0)
        cValue = -1.0;
    if (cValue >  1.0)
        cValue =  1.0;

    vnl_matrix<double> rotationMatrix(3, 3, 0);
    rotationMatrix.put(0, 1, -vVector[2]);
    rotationMatrix.put(1, 0,  vVector[2]);
    rotationMatrix.put(0, 2,  vVector[1]);
    rotationMatrix.put(2, 0, -vVector[1]);
    rotationMatrix.put(1, 2, -vVector[0]);
    rotationMatrix.put(2, 1,  vVector[0]);

    rotationMatrix += rotationMatrix * rotationMatrix / (1.0 + cValue);
    for (unsigned int i = 0;i < 3;++i)
        rotationMatrix(i, i) += 1.0;

    return rotationMatrix;
}

template <class ScalarType>
itk::Matrix <double,3,3>
GetRotationMatrixFromVectors(const itk::Point<ScalarType> &first_direction, const itk::Point<ScalarType> &second_direction)
{
    unsigned int dimension = first_direction.GetPointDimension();
    return GetRotationMatrixFromVectors(first_direction, second_direction, dimension);
}

template <class ScalarType, unsigned int NDimension>
itk::Matrix <double,3,3>
GetRotationMatrixFromVectors(const itk::Vector<ScalarType,NDimension> &first_direction, const itk::Vector<ScalarType,NDimension> &second_direction)
{
    return GetRotationMatrixFromVectors(first_direction, second_direction, NDimension);
}

template <class ScalarType, unsigned int NDimension>
itk::Matrix <double,3,3>
GetRotationMatrixFromVectors(const vnl_vector_fixed<ScalarType,NDimension> &first_direction, const vnl_vector_fixed<ScalarType,NDimension> &second_direction)
{
    return GetRotationMatrixFromVectors(first_direction, second_direction, NDimension);
}

template <class ScalarType, class VectorType>
void
LowerTriangularSolver(vnl_matrix <ScalarType> &matrix, VectorType &rhs, VectorType &result, unsigned int rank)
{
    unsigned int n = matrix.rows();
    if (rank != 0)
        n = rank;

    for (unsigned int i = 0;i < n;++i)
    {
        double resValue = rhs[i];
        for (unsigned int j = 0;j < i;++j)
            resValue -= matrix.get(i,j) * result[j];

        result[i] = resValue / matrix(i,i);
    }
}

template <class ScalarType, class VectorType>
void
UpperTriangularSolver(const vnl_matrix <ScalarType> &matrix, const VectorType &rhs, VectorType &result, unsigned int rank)
{
    unsigned int n = matrix.cols();
    if (rank != 0)
        n = rank;

    for (int i = n - 1;i >= 0;--i)
    {
        double resValue = rhs[i];
        for (unsigned int j = i + 1;j < n;++j)
            resValue -= matrix.get(i,j) * result[j];

        result[i] = resValue / matrix.get(i,i);
    }
}

} // end of namespace anima
