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
#ifndef __itkCorrespondingPointsEuclideanDistancePointStackMetric_h
#define __itkCorrespondingPointsEuclideanDistancePointStackMetric_h

#include "itkAdvancedImageToImageMetric.h"
//#include "itkSingleValuedPointSetToPointSetMetric.h"
//#include "itkPoint.h"
#include "itkPointSet.h"
//#include "itkImage.h"

namespace itk
{

/** \class CorrespondingPointsEuclideanDistancePointStackMetric
 * \brief Computes the Euclidean distance between two point-sets in the moving 
 *  image. The last dimension is ignored, which makes this a useful metric for
 *  landmark-based groupwise registration.
 *
 *
 * \ingroup RegistrationMetrics
 */

template< class TFixedImage, class TMovingImage >
class CorrespondingPointsEuclideanDistancePointStackMetric :
  public AdvancedImageToImageMetric< TFixedImage, TMovingImage >
{
public:

  /** Standard class typedefs. */
  typedef CorrespondingPointsEuclideanDistancePointStackMetric Self;
  typedef AdvancedImageToImageMetric<
    TFixedImage, TMovingImage >                   Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CorrespondingPointsEuclideanDistancePointStackMetric,
    AdvancedImageToImageMetric );

  /** Types transferred from the base classes */
  typedef typename Superclass::FixedImageIndexType      ImageIndexType;
  typedef typename Superclass::FixedImageIndexValueType ImageIndexValueType;
  typedef typename Superclass::FixedImagePointType      ImagePointType;
  typedef typename itk::ContinuousIndex< CoordinateRepresentationType, FixedImageDimension >
    ImageContinuousIndexType;
  typedef typename Superclass::TransformType           TransformType;
  typedef typename Superclass::TransformPointer        TransformPointer;
  typedef typename Superclass::TransformParametersType TransformParametersType;
  typedef typename Superclass::TransformJacobianType   TransformJacobianType;
  typedef typename Superclass::InputPointType          InputPointType;
  typedef typename Superclass::MeasureType             MeasureType;
  typedef typename Superclass::DerivativeType          DerivativeType;
  typedef typename Superclass::DerivativeValueType     DerivativeValueType;
  typedef typename Superclass::NonZeroJacobianIndicesType NonZeroJacobianIndicesType;
  
  typedef typename ImagePointType::CoordRepType                            CoordRepType;
  typedef vnl_vector< CoordRepType >                                       VnlVectorType;
  typedef Matrix< CoordRepType, FixedImageDimension, FixedImageDimension > NdMatrixType;
  typedef Vector< CoordRepType, FixedImageDimension >                      NdVectorType;

    
  /** Typedefs for the point sets. Assuming fixed and moving pointsets are of
      equal type, which implicitly assumes that the fixed and moving image are
      of the same type.
  */
  typedef PointSet< CoordinateRepresentationType,
    TFixedImage::ImageDimension,
    DefaultStaticMeshTraits<
      CoordinateRepresentationType,
      TFixedImage::ImageDimension,
      TFixedImage::ImageDimension,
      CoordinateRepresentationType, CoordinateRepresentationType,
      CoordinateRepresentationType > >                        PointSetType;
  typedef typename PointSetType::Pointer                      PointSetPointer;
  typedef typename PointSetType::ConstPointer                 PointSetConstPointer;
  typedef typename PointSetType::PointsContainerConstIterator PointConstIterator;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & parameters ) const override;

  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType & parameters,
    DerivativeType & Derivative ) const override;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & parameters,
    MeasureType & Value, DerivativeType & Derivative ) const override;

protected:

  CorrespondingPointsEuclideanDistancePointStackMetric();
  ~CorrespondingPointsEuclideanDistancePointStackMetric() override {}

  /** Member variables. */
  PointSetPointer m_FixedPointSet;
  PointSetPointer m_MovingPointSet;

private:

  CorrespondingPointsEuclideanDistancePointStackMetric( const Self & ); // purposely not implemented
  void operator=( const Self & );                                  // purposely not implemented

  void GetValueAndDerivativeFull( const TransformParametersType & parameters,
    MeasureType & Value, DerivativeType & Derivative, bool getDerivative ) const;

  bool GetLocalLinearTransform( const ImagePointType & point, NdMatrixType & transform ) const;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCorrespondingPointsEuclideanDistancePointStackMetric.hxx"
#endif

#endif
