/*=========================================================================
 *
 *  Copyright UMC Utrecht and contributors
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkCorrespondingPointsEuclideanDistancePointStackMetric_hxx
#define __itkCorrespondingPointsEuclideanDistancePointStackMetric_hxx

#include "itkCorrespondingPointsEuclideanDistancePointStackMetric.h"

namespace itk
{

/**
 * ******************* Constructor *******************
 */

template< class TFixedPointSet, class TMovingPointSet >
CorrespondingPointsEuclideanDistancePointStackMetric< TFixedPointSet, TMovingPointSet >
::CorrespondingPointsEuclideanDistancePointStackMetric()
{} // end Constructor

/**
 * ******************* GetValue *******************
 */

template< class TFixedPointSet, class TMovingPointSet >
typename CorrespondingPointsEuclideanDistancePointStackMetric< TFixedPointSet, TMovingPointSet >::MeasureType
CorrespondingPointsEuclideanDistancePointStackMetric< TFixedPointSet, TMovingPointSet >
::GetValue( const TransformParametersType & parameters ) const
{
  /** Sanity checks. */
  FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();
  if( !fixedPointSet )
  {
    itkExceptionMacro( << "Fixed point set has not been assigned" );
  }

  MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();
  if( !movingPointSet )
  {
    itkExceptionMacro( << "Moving point set has not been assigned" );
  }

  /** Initialize some variables. */
  this->m_NumberOfPointsCounted = 0;
  MeasureType     measure = NumericTraits< MeasureType >::Zero;
  InputPointType  pointB;
  OutputPointType pointA, mappedPointA, mappedPointB;
  const unsigned int lastDim = movingPointSet->PointDimension - 1;

  /** Make sure the transform parameters are up to date. */
  this->SetTransformParameters( parameters );

  /** Create iterators. */
  PointIterator pointItFixed  = fixedPointSet->GetPoints()->Begin();
  PointIterator pointItMoving = movingPointSet->GetPoints()->Begin();
  PointIterator pointEnd      = fixedPointSet->GetPoints()->End();

  /** Loop over the corresponding points. */
  while( pointItFixed != pointEnd )
  {
    /** Get the current corresponding points. */
    pointA = pointItFixed.Value();
    pointB = pointItMoving.Value();
    
    /** Calculate a linear approximation of the target point. With a single
        pair of landmarks, the target point should lie right in between point A
        and point B. */
    OutputPointType targetPoint = pointA + 0.5 * ( pointB - pointA );
    targetPoint[lastDim] = pointA[lastDim];
    mappedPointA = this->m_Transform->TransformPoint( targetPoint );
    targetPoint[lastDim] = pointB[lastDim];
    mappedPointB = this->m_Transform->TransformPoint( targetPoint );
    targetPoint += 0.5 * ( pointB - mappedPointB ) + 0.5 * ( pointA - mappedPointA);
    
    /** Now calculate the mappings of the estimated target point onto image A
        and image B. These points should 1) lie close to the given landmarks
        and 2) have the same distance vector to the landmarks. */
    targetPoint[lastDim] = pointA[lastDim];
    mappedPointA = this->m_Transform->TransformPoint( targetPoint );
    targetPoint[lastDim] = pointB[lastDim];
    mappedPointB = this->m_Transform->TransformPoint( targetPoint );
    
    /** Since we are not interested in voxel values, we do not need to check if 
        the points are inside the mask. */
    bool sampleOk = true;
    if( sampleOk )
    {
      this->m_NumberOfPointsCounted++;
      
      /** The error is the difference in distance vector from the mapped target
          points to the landmarks. */
      VnlVectorType vecB = ( pointB - mappedPointB ).GetVnlVector();
      VnlVectorType vecA = ( pointA - mappedPointA ).GetVnlVector();
      measure += (vecA - vecB).magnitude();

    } // end if sampleOk

    ++pointItFixed;
    ++pointItMoving;

  } // end loop over all corresponding points

  return measure / this->m_NumberOfPointsCounted;

} // end GetValue()


/**
 * ******************* GetDerivative *******************
 */

template< class TFixedPointSet, class TMovingPointSet >
void
CorrespondingPointsEuclideanDistancePointStackMetric< TFixedPointSet, TMovingPointSet >
::GetDerivative( const TransformParametersType & parameters,
  DerivativeType & derivative ) const
{
  /** When the derivative is calculated, all information for calculating
   * the metric value is available. It does not cost anything to calculate
   * the metric value now. Therefore, we have chosen to only implement the
   * GetValueAndDerivative(), supplying it with a dummy value variable.
   */
  MeasureType dummyvalue = NumericTraits< MeasureType >::Zero;
  this->GetValueAndDerivative( parameters, dummyvalue, derivative );

} // end GetDerivative()


/**
 * ******************* GetValueAndDerivative *******************
 */

template< class TFixedPointSet, class TMovingPointSet >
void
CorrespondingPointsEuclideanDistancePointStackMetric< TFixedPointSet, TMovingPointSet >
::GetValueAndDerivative( const TransformParametersType & parameters,
  MeasureType & value, DerivativeType & derivative ) const
{
  /** Sanity checks. */
  FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();
  if( !fixedPointSet )
  {
    itkExceptionMacro( << "Fixed point set has not been assigned" );
  }

  MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();
  if( !movingPointSet )
  {
    itkExceptionMacro( << "Moving point set has not been assigned" );
  }

  /** Initialize some variables */
  this->m_NumberOfPointsCounted = 0;
  MeasureType measure = NumericTraits< MeasureType >::Zero;
  derivative = DerivativeType( this->GetNumberOfParameters() );
  derivative.Fill( NumericTraits< DerivativeValueType >::ZeroValue() );
  NonZeroJacobianIndicesType nzji(
    this->m_Transform->GetNumberOfNonZeroJacobianIndices() );
  TransformJacobianType jacobian;

  InputPointType  pointB;
  OutputPointType pointA, mappedPointA, mappedPointB;
  const unsigned int lastDim = movingPointSet->PointDimension - 1;

  /** Call non-thread-safe stuff, such as:
   *   this->SetTransformParameters( parameters );
   *   this->GetImageSampler()->Update();
   * Because of these calls GetValueAndDerivative itself is not thread-safe,
   * so cannot be called multiple times simultaneously.
   * This is however needed in the CombinationImageToImageMetric.
   * In that case, you need to:
   * - switch the use of this function to on, using m_UseMetricSingleThreaded = true
   * - call BeforeThreadedGetValueAndDerivative once (single-threaded) before
   *   calling GetValueAndDerivative
   * - switch the use of this function to off, using m_UseMetricSingleThreaded = false
   * - Now you can call GetValueAndDerivative multi-threaded.
   */
  this->BeforeThreadedGetValueAndDerivative( parameters );

  /** Create iterators. */
  PointIterator pointItFixed  = fixedPointSet->GetPoints()->Begin();
  PointIterator pointItMoving = movingPointSet->GetPoints()->Begin();
  PointIterator pointEnd      = fixedPointSet->GetPoints()->End();

  /** Loop over the corresponding points. */
  while( pointItFixed != pointEnd )
  {
    /** Get the current corresponding points. */
    pointA = pointItFixed.Value();
    pointB = pointItMoving.Value();
    
    /** Calculate a linear approximation of the target point. With a single
        pair of landmarks, the target point should lie right in between point A
        and point B. */
    OutputPointType targetPoint = pointA + 0.5 * ( pointB - pointA );
    targetPoint[lastDim] = pointA[lastDim];
    mappedPointA = this->m_Transform->TransformPoint( targetPoint );
    targetPoint[lastDim] = pointB[lastDim];
    mappedPointB = this->m_Transform->TransformPoint( targetPoint );
    targetPoint += 0.5 * ( pointB - mappedPointB ) + 0.5 * ( pointA - mappedPointA);
    
    /** Now calculate the mappings of the estimated target point onto image A
        and image B. These points should 1) lie close to the given landmarks
        and 2) have the same distance vector to the landmarks. */
    targetPoint[lastDim] = pointA[lastDim];
    mappedPointA = this->m_Transform->TransformPoint( targetPoint );
    targetPoint[lastDim] = pointB[lastDim];
    mappedPointB = this->m_Transform->TransformPoint( targetPoint );

    /** Since we are not interested in voxel values, we do not need to check if 
        the points are inside the mask. */
    bool sampleOk = true;
    if( sampleOk )
    {
      this->m_NumberOfPointsCounted++;
      
      /** The error is the difference in distance vector from the mapped target 
          points to the landmarks. */
      VnlVectorType vecA = ( pointA - mappedPointA ).GetVnlVector();
      VnlVectorType vecB = ( pointB - mappedPointB ).GetVnlVector();
      MeasureType distance  = (vecA - vecB).magnitude();
      measure += distance;

      /** Calculate the contributions to the derivatives with respect to each parameter. */
      if( distance > std::numeric_limits< MeasureType >::epsilon() )
      {
        /** Get the TransformJacobian dT/dmu for point A. */
        targetPoint[lastDim] = pointA[lastDim];
        this->m_Transform->GetJacobian( targetPoint, jacobian, nzji );
        
        /** The mapping for point A should move towards the difference between vecA and vecB. */
        VnlVectorType diff_2 = (vecA - vecB) / distance;
        for( unsigned int i = 0; i < nzji.size(); ++i )
        {
          const unsigned int index  = nzji[ i ];
          VnlVectorType column = jacobian.get_column( i );
          derivative[ index ] -= dot_product( diff_2, column );
        }
        
        /** Get the TransformJacobian dT/dmu for point B. */
        targetPoint[lastDim] = pointB[lastDim];
        this->m_Transform->GetJacobian( targetPoint, jacobian, nzji );
        
        /** The mapping for point B should move in opposite direction. */
        for( unsigned int i = 0; i < nzji.size(); ++i )
        {
          const unsigned int index  = nzji[ i ];
          VnlVectorType column = jacobian.get_column( i );
          derivative[ index ] += dot_product( diff_2, column );
        }
      } // end if distance != 0

    } // end if sampleOk

    ++pointItFixed;
    ++pointItMoving;

  } // end loop over all corresponding points

  /** Check if enough samples were valid. */
//   this->CheckNumberOfSamples(
//     fixedPointSet->GetNumberOfPoints(), this->m_NumberOfPointsCounted );

  /** Copy the measure to value. */
  value = measure;
  if( this->m_NumberOfPointsCounted > 0 )
  {
    derivative /= (2 * this->m_NumberOfPointsCounted);
    value       = measure / this->m_NumberOfPointsCounted;
  }

} // end GetValueAndDerivative()


} // end namespace itk

#endif // end #ifndef __itkCorrespondingPointsEuclideanDistancePointStackMetric_hxx
