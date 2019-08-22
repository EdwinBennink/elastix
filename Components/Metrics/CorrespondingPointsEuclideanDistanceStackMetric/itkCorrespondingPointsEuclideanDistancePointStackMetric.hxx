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

template< class TFixedImage, class TMovingImage >
CorrespondingPointsEuclideanDistancePointStackMetric< TFixedImage, TMovingImage >
::CorrespondingPointsEuclideanDistancePointStackMetric()
{
	this->SetUseImageSampler(true);
	this->SetUseFixedImageLimiter(false);
	this->SetUseMovingImageLimiter(false);

  this->m_FixedPointSet = NULL;
  this->m_MovingPointSet = NULL;
} // end Constructor

/**
 * ******************* GetValue *******************
 */

template< class TFixedImage, class TMovingImage >
typename CorrespondingPointsEuclideanDistancePointStackMetric< TFixedImage, TMovingImage >::MeasureType
CorrespondingPointsEuclideanDistancePointStackMetric< TFixedImage, TMovingImage >
::GetValue( const TransformParametersType & parameters ) const
{
  MeasureType value = NumericTraits< MeasureType >::Zero;
  DerivativeType dummyDerivative;
  this->GetValueAndDerivativeFull( parameters, value, dummyDerivative, false );
  return value;
} // end GetValue()


/**
 * ******************* GetDerivative *******************
 */

template< class TFixedImage, class TMovingImage >
void
CorrespondingPointsEuclideanDistancePointStackMetric< TFixedImage, TMovingImage >
::GetDerivative( const TransformParametersType & parameters,
  DerivativeType & derivative ) const
{
  /** When the derivative is calculated, all information for calculating
   * the metric value is available. It does not cost anything to calculate
   * the metric value now. Therefore, we have chosen to only implement the
   * GetValueAndDerivative(), supplying it with a dummy value variable.
   */
  MeasureType dummyvalue = NumericTraits< MeasureType >::Zero;
  this->GetValueAndDerivativeFull( parameters, dummyvalue, derivative, true );

} // end GetDerivative()


/**
 * ******************* GetValueAndDerivative *******************
 */

template< class TFixedImage, class TMovingImage >
void
CorrespondingPointsEuclideanDistancePointStackMetric< TFixedImage, TMovingImage >
::GetValueAndDerivative( const TransformParametersType & parameters,
  MeasureType & value, DerivativeType & derivative ) const
{
  this->GetValueAndDerivativeFull( parameters, value, derivative, true );
} // end GetValueAndDerivative


template< class TFixedImage, class TMovingImage >
void
CorrespondingPointsEuclideanDistancePointStackMetric< TFixedImage, TMovingImage >
::GetValueAndDerivativeFull( const TransformParametersType & parameters,
  MeasureType & value, DerivativeType & derivative, bool getDerivative ) const
{
  /** Sanity checks. */
  PointSetConstPointer fixedPointSet = this->m_FixedPointSet;
  if( !fixedPointSet )
  {
    itkExceptionMacro( << "Fixed point set has not been assigned" );
  }

  PointSetConstPointer movingPointSet = this->m_MovingPointSet;
  if( !movingPointSet )
  {
    itkExceptionMacro( << "Moving point set has not been assigned" );
  }

  /** Initialize some variables */
  value = NumericTraits< MeasureType >::Zero;
  if( getDerivative )
  {
    derivative = DerivativeType( this->GetNumberOfParameters() );
    derivative.Fill( NumericTraits< DerivativeValueType >::ZeroValue() );
  }
  NonZeroJacobianIndicesType nzji(
    this->m_AdvancedTransform->GetNumberOfNonZeroJacobianIndices() );
  TransformJacobianType jacobian;
  const unsigned int lastDim     = FixedImageDimension - 1;
  const unsigned int lastDimSize = this->GetFixedImage()->GetLargestPossibleRegion().GetSize( lastDim );

  /** Calculate a normalization factor for the derivatives. */
  ImageRegion< FixedImageDimension > region         = this->GetFixedImage()->GetLargestPossibleRegion();
  NdVectorType                       spacing        = this->GetFixedImage()->GetSpacing();
  float                              diagonalSquare = 0.0;
  for( unsigned int d = 0; d < lastDim; ++d )
  {
    diagonalSquare += itk::Math::sqr( spacing[ d ] * region.GetSize( d ) );
  }
  const float normalizationFactor = 1.0 / vcl_sqrt( diagonalSquare );
 
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

  ImageSampleContainerPointer sampleContainer = this->GetImageSampler()->GetOutput();
  const unsigned int          sampleSize = sampleContainer->Size();
  const unsigned int          pointSetSize = fixedPointSet->GetNumberOfPoints();
  std::vector< float >        aSquaredDistance( pointSetSize, std::numeric_limits< float >::infinity() );
  std::vector< float >        bSquaredDistance( pointSetSize, std::numeric_limits< float >::infinity() );
  std::vector< unsigned int > aClosestPoint( pointSetSize );
  std::vector< unsigned int > bClosestPoint( pointSetSize );
  std::vector< bool >         aPointOk( pointSetSize, false );
  std::vector< bool >         bPointOk( pointSetSize, false );

  /** Iterate over all landmarks and find the closest sample point. */
  for( unsigned int s = 0; s < sampleSize; ++s )
  {
    /** Read fixed coordinates and treat the last dimension as stack index. */
    ImagePointType pMoving;
    ImagePointType pFixed = sampleContainer->GetElement( s ).m_ImageCoordinates;
    pFixed[ lastDim ] = itk::Math::Round< CoordRepType >( pFixed[ lastDim ] );

    /** Map the samplepoint to moving coordinates. */
    if( !this->TransformPoint( pFixed, pMoving ) )
    {
      continue;
    }

    for( unsigned int p = 0; p < pointSetSize; ++p )
    {
      /** Get the current landmark point. */
      InputPointType pointA = fixedPointSet->GetPoint( p );
      InputPointType pointB = movingPointSet->GetPoint( p );

      /** Check if the sample is the closest one yet to point A or B. */
      if( pMoving[ lastDim ] == pointA[ lastDim ] )
      {
        float sd = ( pointA - pMoving ).GetSquaredNorm();
        if( sd < aSquaredDistance[ p ] )
        {
          aPointOk[ p ] = true;
          aClosestPoint[ p ] = s;
          aSquaredDistance[ p ] = sd;
        }
      }
      else if( pMoving[ lastDim ] == pointB[ lastDim ] )
      {
        float sd = ( pointB - pMoving ).GetSquaredNorm();
        if( sd < bSquaredDistance[ p ] )
        {
          bPointOk[ p ] = true;
          bClosestPoint[ p ] = s;
          bSquaredDistance[ p ] = sd;
        }
      }
    }
  }

  /** Calculate local transformation matrices and estimate the position of the landmarks
      in the fixed image. */
  std::vector< NdVectorType > aFixed( pointSetSize );
  std::vector< NdVectorType > bFixed( pointSetSize );
  std::vector< NdMatrixType > aTransformFtoM( pointSetSize );
  std::vector< NdMatrixType > bTransformFtoM( pointSetSize );
  for( unsigned int p = 0; p < pointSetSize; ++p )
  {
    if( aPointOk[ p ] && bPointOk[ p ] )
    {
      /** Calculate local transformation matrix for point A. */
      ImagePointType pFixed = sampleContainer->GetElement( aClosestPoint[ p ] ).m_ImageCoordinates;
      aPointOk[ p ] = this->GetLocalLinearTransform( pFixed, aTransformFtoM[ p ] );
    }

    if( aPointOk[ p ] && bPointOk[ p ] )
    {
      /** Transform point A to fixed coordinates. */
      NdVectorType aMoving = fixedPointSet->GetPoint( p ).GetVectorFromOrigin();
      aMoving[ lastDim ] = 1.0;
      aFixed[ p ] = NdMatrixType( aTransformFtoM[ p ].GetInverse() ) * aMoving;
    }

    if( aPointOk[ p ] && bPointOk[ p ] )
    {
      /** Calculate local transformation matrix for point B. */
      ImagePointType pFixed = sampleContainer->GetElement( bClosestPoint[ p ] ).m_ImageCoordinates;
      bPointOk[ p ] = this->GetLocalLinearTransform( pFixed, bTransformFtoM[ p ] );
    }

    if( aPointOk[ p ] && bPointOk[ p ] )
    {
      /** Transform point B to fixed coordinates. */
      NdVectorType bMoving = movingPointSet->GetPoint( p ).GetVectorFromOrigin();
      bMoving[ lastDim ] = 1.0;
      bFixed[ p ] = NdMatrixType( bTransformFtoM[ p ].GetInverse() ) * bMoving;
    }
  }

  /** Calculate the errors and derivatives for the landmark pairs. */
  unsigned int numSamplesOk = 0;
  for( unsigned int p = 0; p < pointSetSize; ++p )
  {
    if( aPointOk[ p ] && bPointOk[ p ] )
    {
      numSamplesOk++;

      /** The error is given by the squared distance between the mappings of point A
      and point B in the fixed image. */
      MeasureType sd = ( aFixed[ p ] - bFixed[ p ] ).GetSquaredNorm();
      value += sd;

      if( getDerivative )
      {
        /** We now need to set the derivative at the projection of point A in fixed
            coordinates. This derivative vector points from point A in the moving 
            image to the projection of point B from fixed coordinates to moving 
            coordinates in point A's local linear space. Vice versa for point B. */

        /** Get the current landmark points. */
        InputPointType pointA = fixedPointSet->GetPoint( p );
        InputPointType pointB = movingPointSet->GetPoint( p );

        /** Point A. Divide by sqrt( norm ) to boost small gradients. */
        NdVectorType aMoving     = pointA.GetVectorFromOrigin();
        aMoving[ lastDim ]       = 1.0;
        NdVectorType aTarget     = bFixed[ p ];
        NdVectorType aDerivative = ( aTransformFtoM[ p ] * aTarget - aMoving ) * normalizationFactor;
        float        norm = aDerivative.GetNorm();
        if( norm > std::numeric_limits< MeasureType >::epsilon() )
        {
          aDerivative /= vcl_sqrt( norm );

          ImagePointType pFixed = ImagePointType( aFixed[ p ] );
          pFixed[ lastDim ] = itk::Math::Round< CoordRepType >( pointA[ lastDim ] );
          this->m_AdvancedTransform->GetJacobian( pFixed, jacobian, nzji );

          for( unsigned int i = 0; i < nzji.size(); ++i )
          {
            const unsigned int index = nzji[ i ];
            VnlVectorType column = jacobian.get_column( i );
            derivative[ index ] += dot_product( aDerivative.GetVnlVector(), column );
          }
        }

        /** Point B. Divide by sqrt( norm ) to boost small gradients. */
        NdVectorType bMoving     = pointB.GetVectorFromOrigin();
        bMoving[ lastDim ]       = 1.0;
        NdVectorType bTarget     = aFixed[ p ];
        NdVectorType bDerivative = ( bTransformFtoM[ p ] * bTarget - bMoving ) * normalizationFactor;
        norm                     = bDerivative.GetNorm();
        if( norm > std::numeric_limits< MeasureType >::epsilon() )
        {
          bDerivative /= vcl_sqrt( norm );

          ImagePointType pFixed = ImagePointType( bFixed[ p ] );
          pFixed[ lastDim ] = itk::Math::Round< CoordRepType >( pointB[ lastDim ] );
          this->m_AdvancedTransform->GetJacobian( pFixed, jacobian, nzji );

          for( unsigned int i = 0; i < nzji.size(); ++i )
          {
            const unsigned int index = nzji[ i ];
            VnlVectorType column = jacobian.get_column( i );
            derivative[ index ] += dot_product( bDerivative.GetVnlVector(), column );
          }
        }
      }
    }
  }

  /** Normalize the number of samples and the derivative. */
  value /= numSamplesOk;
  if( getDerivative )
  {
    derivative /= numSamplesOk;
  }

} // end GetValueAndDerivativeFull()


template< class TFixedImage, class TMovingImage >
bool
CorrespondingPointsEuclideanDistancePointStackMetric< TFixedImage, TMovingImage >
::GetLocalLinearTransform( const ImagePointType & point, NdMatrixType & transform ) const
{
  bool               pointOk = true;
  const unsigned int lastDim = FixedImageDimension - 1;
  const NdVectorType spacing = this->GetFixedImage()->GetSpacing();
  NdMatrixType       mFixed, mMoving;

  /** Get the fixed and moving image coordinates of nd points at 0.5 pixel spacing
      from the given point and use these to calculate a local transformation
      matrix. Break the for-loop if a transformation to moving coordinates failed. */
  for( unsigned int d = 0; d < FixedImageDimension && pointOk; ++d )
  {
    ImagePointType pFixed = point;
    ImagePointType pMoving;

    pFixed[ lastDim ] = itk::Math::Round< CoordRepType >( pFixed[ lastDim ] );
    if( d > 0 ) {
      pFixed[ d - 1 ] += spacing[ d - 1 ];
    }
    pointOk = this->TransformPoint( pFixed, pMoving );

    /** Add the coordinates of the fixed and moving image points to the matrices and
        set the last column to 1 (translation part of the transformation matrix). */
    for( unsigned int c = 0; c < lastDim; ++c )
    {
      mFixed( d, c ) = pFixed[ c ];
      mMoving( d, c ) = pMoving[ c ];
    }
    mFixed( d, lastDim ) = 1.0;
    mMoving( d, lastDim ) = 1.0;
  }

  if( pointOk )
  {
    transform = ( NdMatrixType( mFixed.GetInverse() ) * mMoving ).GetTranspose();
  }

  return pointOk;
}


} // end namespace itk

#endif // end #ifndef __itkCorrespondingPointsEuclideanDistancePointStackMetric_hxx
