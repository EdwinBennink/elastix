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
#ifndef __elxCorrespondingPointsEuclideanDistanceStackMetric_H__
#define __elxCorrespondingPointsEuclideanDistanceStackMetric_H__

#include "elxIncludes.h" // include first to avoid MSVS warning
#include "itkCorrespondingPointsEuclideanDistancePointStackMetric.h"

namespace elastix
{

/**
 * \class CorrespondingPointsEuclideanDistanceStackMetric
 * \brief An metric based on the itk::CorrespondingPointsEuclideanDistancePointStackMetric.
 *
 * The parameters used in this class are:
 * \parameter Metric: Select this metric as follows:\n
 *    <tt>(Metric "CorrespondingPointsEuclideanDistanceStackMetric")</tt>
 *
 * \ingroup Metrics
 *
 */

template< class TElastix >
class CorrespondingPointsEuclideanDistanceStackMetric :
  public
  itk::CorrespondingPointsEuclideanDistancePointStackMetric<
    typename MetricBase< TElastix >::FixedImageType,
    typename MetricBase< TElastix >::MovingImageType >,
  public MetricBase< TElastix >
{
public:

  /** Standard ITK-stuff. */
  typedef CorrespondingPointsEuclideanDistanceStackMetric Self;
  typedef itk::CorrespondingPointsEuclideanDistancePointStackMetric<
    typename MetricBase< TElastix >::FixedImageType,
    typename MetricBase< TElastix >::MovingImageType > Superclass1;
  typedef MetricBase< TElastix >          Superclass2;
  typedef itk::SmartPointer< Self >       Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CorrespondingPointsEuclideanDistanceStackMetric,
    itk::CorrespondingPointsEuclideanDistancePointStackMetric );

  /** Name of this class.
   * Use this name in the parameter file to select this specific metric. \n
   * example: <tt>(Metric "CorrespondingPointsEuclideanDistanceStackMetric")</tt>\n
   */
  elxClassNameMacro( "CorrespondingPointsEuclideanDistanceStackMetric" );

  /** Typedefs from the superclass. */
  typedef typename Superclass1::CoordinateRepresentationType CoordinateRepresentationType;
  typedef typename Superclass1::PointSetType            PointSetType;
  typedef typename Superclass1::PointSetPointer         PointSetPointer;
  typedef typename Superclass1::PointSetConstPointer    PointSetConstPointer;

  typedef typename Superclass1::TransformType           TransformType;
  typedef typename Superclass1::TransformPointer        TransformPointer;
  typedef typename Superclass1::InputPointType          InputPointType;
  typedef typename Superclass1::OutputPointType         OutputPointType;
  typedef typename Superclass1::TransformParametersType TransformParametersType;
  typedef typename Superclass1::TransformJacobianType   TransformJacobianType;
  typedef typename Superclass1::FixedImageType         FixedImageType;
  typedef typename Superclass1::FixedImageMaskType     FixedImageMaskType;
  typedef typename Superclass1::FixedImageMaskPointer  FixedImageMaskPointer;
  typedef typename Superclass1::MovingImageType        MovingImageType;
  typedef typename Superclass1::MovingImageMaskType    MovingImageMaskType;
  typedef typename Superclass1::MovingImageMaskPointer MovingImageMaskPointer;
  typedef typename Superclass1::MeasureType            MeasureType;
  typedef typename Superclass1::DerivativeType         DerivativeType;
  typedef typename Superclass1::ParametersType         ParametersType;

  /** Typedefs inherited from elastix. */
  typedef typename Superclass2::ElastixType          ElastixType;
  typedef typename Superclass2::ElastixPointer       ElastixPointer;
  typedef typename Superclass2::ConfigurationType    ConfigurationType;
  typedef typename Superclass2::ConfigurationPointer ConfigurationPointer;
  typedef typename Superclass2::RegistrationType     RegistrationType;
  typedef typename Superclass2::RegistrationPointer  RegistrationPointer;
  typedef typename Superclass2::ITKBaseType          ITKBaseType;

  /** The fixed image dimension. */
  itkStaticConstMacro( FixedImageDimension, unsigned int,
    FixedImageType::ImageDimension );

  /** The moving image dimension. */
  itkStaticConstMacro( MovingImageDimension, unsigned int,
    MovingImageType::ImageDimension );

  /** Sets up a timer to measure the initialization time and calls the
   * Superclass' implementation.
   */
  void Initialize( void ) override;

  /**
   * Do some things before all:
   * \li Check and print the command line arguments fp and mp.
   *   This should be done in BeforeAllBase and not BeforeAll.
   */
  int BeforeAllBase( void ) override;

  /**
   * Do some things before registration:
   * \li Load and set the pointsets.
   */
  void BeforeRegistration( void ) override;

  /** Function to read the corresponding points. */
  unsigned int ReadLandmarks(
    const std::string & landmarkFileName,
    typename PointSetPointer & pointSet,
    const typename FixedImageType::ConstPointer image );

  /** Overwrite to silence warning. */
  void SelectNewSamples( void ) override{}

protected:

  /** The constructor. */
  CorrespondingPointsEuclideanDistanceStackMetric(){}
  /** The destructor. */
  ~CorrespondingPointsEuclideanDistanceStackMetric() override {}

private:

  /** The private constructor. */
  CorrespondingPointsEuclideanDistanceStackMetric( const Self & ); // purposely not implemented
  /** The private copy constructor. */
  void operator=( const Self & );              // purposely not implemented

};

} // end namespace elastix

#ifndef ITK_MANUAL_INSTANTIATION
#include "elxCorrespondingPointsEuclideanDistanceStackMetric.hxx"
#endif

#endif // end #ifndef __elxCorrespondingPointsEuclideanDistanceStackMetric_H__
