// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKINTERSECTIONIT_HH
#define DUNE_GRID_TEST_CHECKINTERSECTIONIT_HH

#include <cmath>

#include <dune/common/test/iteratortest.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include "checkiterators.hh"
#include "checkgeometry.hh"

/** \file
    \brief Tests for the IntersectionIterator
 */

struct CheckIntersectionIteratorErrorState
{
  unsigned int sumNormalsNonZero;

  CheckIntersectionIteratorErrorState ()
    : sumNormalsNonZero( 0 )
  {}
};


template< class Grid >
struct EnableIntersectionIteratorReverseCheck
{
  static const bool v = true;
};


// Check that normal and normal2 pointing in the same direction
template< class ctype, int dimworld, class String >
inline void checkParallel ( const Dune::FieldVector< ctype, dimworld > &normal,
                            const Dune::FieldVector< ctype, dimworld > &refNormal,
                            const String & name )
{
  using std::sqrt;
  if( (normal.two_norm()*refNormal.two_norm() - normal*refNormal) > sqrt(std::numeric_limits< ctype >::epsilon()) )
  {
    DUNE_THROW(Dune::GridError)
          << "Error: " << name << " does not point in the direction of outer normal.\n"
          << "       " << name << " = " << normal << ", outer normal = " << refNormal << std :: endl;
  }
}


// Check whether the normal is orthogonal to the intersection, i.e.,
// whether (J^-1 * n) = 0. Here J is the jacobian of the intersection
// geometry (geometry) and n is a normal.
template< class ctype, int dimworld, class JacobianInverseTransposed, class String >
inline void checkJIn ( const Dune :: FieldVector< ctype, dimworld > &normal,
                       const JacobianInverseTransposed &jit,
                       const String & name )
{
  using std::sqrt;
  const int facedim = JacobianInverseTransposed::cols;
  Dune :: FieldVector< ctype, facedim > x( ctype( 0 ) );
  jit.umtv( normal, x );
  if (x.infinity_norm() > sqrt(std::numeric_limits< ctype >::epsilon()))
  {
    const Dune::FieldMatrix< ctype, dimworld, facedim > &mjit = jit;
    DUNE_THROW(Dune::GridError)
                << "Error: The normal vector is not orthogonal to the intersection geometry (J^-1 * n != 0).\n"
                << "       Normal vector: " << name << " = " << normal << "\n"
                << "       Transposed inverse Jacobian (J^-1^T): \n" << mjit << std::endl;
  }
}


template< class GridView >
void checkIntersectionAssignment ( const GridView &view, const typename GridView::template Codim< 0 >::Iterator &l )
{
  typedef typename GridView::template Codim< 0 >::Iterator Iterator;
  typedef typename GridView::IntersectionIterator IntersectionIterator;

  Iterator l1( l );
  Iterator l2( l );
  ++l2;

  if( l2 != view.template end< 0 >() )
  {
    IntersectionIterator it1 = view.ibegin( *l1 );
    IntersectionIterator it2 = view.ibegin( *l2 );

    // verify prerequisites
    if( l1 == l2 )
      DUNE_THROW(Dune::GridError, "Error: Iterators l1 and l2 should not be equal during intersection assignment check.");
    if( it1 == it2 )
      DUNE_THROW(Dune::GridError, "Error: Intersection iterators it1 and it2 should not be equal during intersection assignment check.");
    if( it1->inside() != *l1 )
      DUNE_THROW(Dune::GridError, "Error: Intersection iterator it1 does not point to the expected inside entity.");
    if( it2->inside() != *l2 )
      DUNE_THROW(Dune::GridError, "Error: Intersection iterator it2 does not point to the expected inside entity.");

    // check assignment
    it1 = it2;
    if( it1 != it2 )
      DUNE_THROW(Dune::GridError, "Error: Assignment of intersection iterator it2 to it1 failed.");
    if( it1->inside() != *l2 )
      DUNE_THROW(Dune::GridError, "Error: After assignment, intersection iterator it1 does not point to the expected inside entity.");
  }
}


template< class Intersection >
void checkIntersection ( const Intersection &intersection, bool isCartesian = false )
{
  using std::sqrt;
  typedef typename Intersection::ctype ctype;

  typedef typename Intersection::Entity Entity;
  typedef typename Intersection::LocalGeometry LocalGeometry;
  typedef typename Intersection::Geometry Geometry;

  [[maybe_unused]] const int dimension = Entity::dimension;
  const int mydimension = Intersection::mydimension;

  // check consistency of exported types

  static_assert((std::is_same< ctype, typename Entity::Geometry::ctype >::value),
                "Type Intersection::ctype differs from Intersection::Entity::ctype.");
  static_assert((std::is_same< ctype, typename LocalGeometry::ctype >::value),
                "Type Intersection::ctype differs from Intersection::LocalGeometry::ctype.");
  static_assert((std::is_same< ctype, typename Geometry::ctype >::value),
                "Type Intersection::ctype differs from Intersection::Geometry::ctype.");

  const ctype tolerance = sqrt(std::numeric_limits< ctype >::epsilon());

  // cache some information on the intersection

  const int indexInInside = intersection.indexInInside();

  // obtain inside entity, its geometry and reference element

  const Entity inside = intersection.inside();

  const typename Entity::Geometry insideGeometry = inside.geometry();
  auto refElement = referenceElement( insideGeometry );

  // check that boundary id has positive value and that intersection is conforming

  if( intersection.boundary() )
  {
    if( !intersection.conforming() && !intersection.neighbor() )
      DUNE_THROW(Dune::GridError, "Error: Boundary intersections must be conforming or have a neighboring entity.");
  }

  // check the geometry

  const Geometry geometry = intersection.geometry();
  checkGeometry( geometry );

  if( geometry.type() != intersection.type() )
  {
    DUNE_THROW(Dune::GridError) << "Error: The type of the intersection geometry does not match the reference geometry type." << std::endl;
  }

  // check the geometryInInside

  if( !inside.type().isNone() )
  {
    const LocalGeometry geometryInInside = intersection.geometryInInside();
    checkLocalGeometry( geometryInInside, inside.type(), "geometryInInside" );

    if( geometryInInside.corners() != geometry.corners() )
      DUNE_THROW(Dune::GridError, "Error: The number of corners in geometryInInside does not match the intersection geometry.");

    // create quadrature rule as a set of test points

    const Dune::QuadratureType::Enum qt = Dune::QuadratureType::GaussLegendre;
    const Dune::QuadratureRule< ctype, mydimension > &quadrature
      = Dune::QuadratureRules< ctype, mydimension >::rule( intersection.type(), 3, qt );

    for( std::size_t i = 0; i < quadrature.size(); ++i )
    {
      const Dune::FieldVector< ctype, mydimension > &pt = quadrature[ i ].position();

      typename Geometry::GlobalCoordinate globalPos = geometry.global( pt );
      typename Geometry::GlobalCoordinate localPos = insideGeometry.global( geometryInInside.global( pt ) );

      if( (globalPos - localPos).infinity_norm() > tolerance )
      {
        std::cerr << "Error: Intersection's geometry is inconsistent with concatenation of inside entity's geometry and intersection's geometryInInside." << std::endl;

        std::cerr << "       inside()->geometry().global( intersection.geometryInInside().global( x ) ) != intersection.geometry().global( x )" << std::endl;
        std::cerr << "       x = " << pt << std::endl;
        std::cerr << "       intersection.geometry() = " << geometry.corner( 0 );
        for( int i = 1; i < geometry.corners(); ++i )
          std::cerr << " | " << geometry.corner( i );
        std::cerr << std::endl;

        std::cerr << "       intersection.geometryInInside() = " << geometryInInside.corner( 0 );
        for( int i = 1; i < geometryInInside.corners(); ++i )
          std::cerr << " | " << geometryInInside.corner( i );
        std::cerr << std::endl;

        std::cerr << "       inside()->geometry() = " << insideGeometry.corner( 0 );
        for( int i = 1; i < insideGeometry.corners(); ++i )
          std::cerr << " | " << insideGeometry.corner( i );
        std::cerr << std::endl;

        DUNE_THROW(Dune::GridError, "Intersection error");
      }
    }
  }

  // check geometryInOutside

  if( intersection.neighbor() )
  {

    const Entity outside = intersection.outside();

    if( !outside.type().isNone() )
    {
      const LocalGeometry geometryInOutside = intersection.geometryInOutside();
      checkLocalGeometry( geometryInOutside, outside.type(), "geometryInOutside" );
      if( geometryInOutside.corners() != geometry.corners() )
        DUNE_THROW(Dune::GridError, "Error: The number of corners in geometryInOutside does not match the intersection geometry.");

      if( !intersection.boundary() )
      {
        const typename Entity::Geometry outsideGeometry = outside.geometry();

        // create quadrature rule as a set of test points

        const Dune::QuadratureType::Enum qt = Dune::QuadratureType::GaussLegendre;
        const Dune::QuadratureRule< ctype, mydimension > &quadrature
          = Dune::QuadratureRules< ctype, mydimension >::rule( intersection.type(), 3, qt );

        for( std::size_t i = 0; i < quadrature.size(); ++i )
        {
          const Dune::FieldVector< ctype, mydimension > &pt = quadrature[ i ].position();

          typename Geometry::GlobalCoordinate globalPos = geometry.global( pt );
          typename Geometry::GlobalCoordinate localPos = outsideGeometry.global( geometryInOutside.global( pt ) );

          if( (globalPos - localPos).infinity_norm() > tolerance )
          {
            std::cerr << "Error: Intersection's geometry is inconsistent with concatenation of outside entity's geometry and intersection's geometryInOutside." << std::endl;

            std::cerr << "       outside()->geometry().global( intersection.geometryInOutside().global( x ) ) != intersection.geometry().global( x )" << std::endl;
            std::cerr << "       x = " << pt << std::endl;
            std::cerr << "       intersection.geometry() = " << geometry.corner( 0 );
            for( int i = 1; i < geometry.corners(); ++i )
              std::cerr << " | " << geometry.corner( i );
            std::cerr << std::endl;

            std::cerr << "       intersection.geometryInOutside() = " << geometryInOutside.corner( 0 );
            for( int i = 1; i < geometryInOutside.corners(); ++i )
              std::cerr << " | " << geometryInOutside.corner( i );
            std::cerr << std::endl;

            std::cerr << "       outside()->geometry() = " << outsideGeometry.corner( 0 );
            for( int i = 1; i < outsideGeometry.corners(); ++i )
              std::cerr << " | " << outsideGeometry.corner( i );
            std::cerr << std::endl;

            DUNE_THROW(Dune::GridError, "Intersection error");
          }
        }
      }
    }
  }

  // check normal vectors

  if( !intersection.type().isNone() )
  {

    // create quadrature rule as a set of test points

    const Dune::QuadratureType::Enum qt = Dune::QuadratureType::GaussLegendre;
    const Dune::QuadratureRule< ctype, mydimension > &quadrature
      = Dune::QuadratureRules< ctype, mydimension >::rule( intersection.type(), 3, qt );

    for( std::size_t i = 0; i < quadrature.size(); ++i )
    {
      const Dune::FieldVector< ctype, mydimension > &pt = quadrature[ i ].position();
      const typename Geometry::JacobianInverseTransposed &jit = geometry.jacobianInverseTransposed( pt );

      // independently calculate the integration outer normal for the inside entity

      typename Geometry::GlobalCoordinate refIntNormal;
      if( !inside.type().isNone() )
      {
        const LocalGeometry geometryInInside = intersection.geometryInInside();
        const typename LocalGeometry::GlobalCoordinate xInside = geometryInInside.global( pt );
        const typename LocalGeometry::GlobalCoordinate &refNormal = refElement.integrationOuterNormal( indexInInside );
        insideGeometry.jacobianInverseTransposed( xInside ).mv( refNormal, refIntNormal );

        // note: refElement.template geometry< 1 >( indexInInside ) is affine,
        //       hence we may use any point to obtain the integrationElement.
        refIntNormal *= insideGeometry.integrationElement( xInside ) * geometryInInside.integrationElement( pt )
                        / refElement.template geometry< 1 >( indexInInside ).integrationElement( pt );
      }

      // check outer normal
      const typename Intersection::GlobalCoordinate normal
        = intersection.outerNormal( pt );
      if( !inside.type().isNone() )
        checkParallel( normal, refIntNormal, "outerNormal" );

      // check normal vector is orthogonal to all vectors connecting the vertices
      // note: This only makes sense for non-curved faces. As there is no method
      //       to check curvature, we only check this for affine surface.
      if( geometry.affine() )
      {
        for( int c = 1; c < geometry.corners(); ++c )
        {
          typename Geometry::GlobalCoordinate x = geometry.corner( c-1 );
          x -= geometry.corner( c );
          if( x*normal >= tolerance )
          {
            DUNE_THROW(Dune::GridError) << "Error: The outer normal is not orthogonal to the line connecting the corners of the geometry.";
          }
        }
      }

      // check consistency with JacobianInverseTransposed (J^-1 * n) = 0,
      // where J denotes the Jacobian of the intersection's geometry and n is a
      // normal.
      checkJIn( normal, jit, "outerNormal" );

      // check integration outer normal

      const typename Intersection::GlobalCoordinate intNormal = intersection.integrationOuterNormal( pt );
      const ctype det = geometry.integrationElement( pt );
      if( std::abs( det - intNormal.two_norm() ) > tolerance )
      {
        DUNE_THROW(Dune::GridError)
                  << "Error: The length of the integrationOuterNormal does not match the integration element.\n"
                  << "       |integrationOuterNormal| = " << intNormal.two_norm()
                  << ", integrationElement = " << det << std::endl;
      }

      checkJIn( intNormal, jit, "integrationOuterNormal" );
      if( !inside.type().isNone() )
      {
        checkParallel( intNormal, refIntNormal, "integrationOuterNormal" );

        if( (intNormal - refIntNormal).two_norm() > tolerance )
        {
          std::cerr << "Error: Wrong integration outer normal (" << intNormal
                    << ", should be " << refIntNormal << ")." << std::endl;
          std::cerr << "       Intersection: " << geometry.corner( 0 );
          for( int c = 1; c < geometry.corners(); ++c )
            std::cerr << ", " << geometry.corner( c );
          std::cerr << "." << std::endl;
          DUNE_THROW(Dune::GridError, "Intersection error");
        }
      }

      // check unit outer normal

      const typename Intersection::GlobalCoordinate unitNormal = intersection.unitOuterNormal( pt );
      if( std::abs( ctype( 1 ) - unitNormal.two_norm() ) > tolerance )
      {
        DUNE_THROW(Dune::GridError)
                  << "Error: The length of the unitOuterNormal is not equal to 1.\n"
                  << "       |unitOuterNormal| = " << unitNormal.two_norm() << std::endl;
      }

      checkJIn( unitNormal, jit, "unitOuterNormal" );
      if( !inside.type().isNone() )
        checkParallel( unitNormal, refIntNormal, "unitOuterNormal" );

      // check normal for Cartesian grids

      if( isCartesian )
      {
        if( !geometry.affine() )
          DUNE_THROW( Dune::GridError, "Intersection's geometry is not affine, although isCartesian is true" );

        // check that unit normal of Cartesian grid is given in a specific way
        typename Intersection::GlobalCoordinate normal( 0 );
        normal[ indexInInside / 2 ] = 2 * (indexInInside % 2) - 1;

        if( (normal - unitNormal).infinity_norm() > tolerance )
          DUNE_THROW( Dune::GridError, "Error: The unit normal is not in Cartesian format, although isCartesian is true." );
      }
    }

    // check center unit outer normal

    auto refFace = referenceElement( geometry );

    if( (intersection.centerUnitOuterNormal() - intersection.unitOuterNormal( refFace.position( 0, 0 ) )).two_norm() > tolerance )
    {
      DUNE_THROW(Dune::GridError)
                << "Error: The centerUnitOuterNormal does not match the unitOuterNormal at the reference face position.\n"
                << "       centerUnitOuterNormal = " << intersection.centerUnitOuterNormal()
                << ", unitOuterNormal( " << refFace.position( 0, 0 ) << " ) = "
                << intersection.unitOuterNormal( refFace.position( 0, 0 ) ) << std::endl;
    }
  }
}


/** \brief Test the IntersectionIterator
 */
template< class GridViewType, class ErrorState >
void checkIntersectionIterator ( const GridViewType &view,
                                 const typename GridViewType::template Codim< 0 >::Iterator &eIt,
                                 ErrorState &errorState )
{
  using std::sqrt;
  typedef typename GridViewType::Grid GridType;
  typedef typename GridViewType::IntersectionIterator IntersectionIterator;
  typedef typename GridViewType::Intersection Intersection;

  typedef typename GridType::ctype ctype;

  typedef typename Intersection::Entity Entity;

  const bool isCartesian = Dune::Capabilities::isCartesian< GridType >::v;
  const bool checkOutside = EnableIntersectionIteratorReverseCheck< GridType >::v;

  // check consistency of exported types

  static_assert( (std::is_same< ctype, typename Intersection::ctype >::value),
                      "Type GridView::Grid::ctype differs from GridView::Intersection::ctype." );

  static_assert( (std::is_same< Intersection, typename IntersectionIterator::Intersection >::value),
                      "Type GridView::Intersection differs from GridView::IntersectionIterator::Intersection." );

  static_assert((static_cast<int>(Intersection::dimensionworld)
                      == static_cast<int>(GridType::dimensionworld)),"IntersectionIterator has wrong dimensionworld");


  // check default constructibility of intersections
  [[maybe_unused]] Intersection default_construct_intersection;

  // initialize variables for element checks

  bool hasBoundaryIntersection = false;
  typename Intersection::GlobalCoordinate sumNormal( ctype( 0 ) );

  // check whether intersection iterator is a forward iterator
  NoopFunctor< Intersection > op;
  if( 0 != testForwardIterator( view.ibegin( *eIt ), view.iend( *eIt ), op ) )
    DUNE_THROW( Dune::Exception, "IntersectionIterator does not fulfill the forward iterator concept" );

  const IntersectionIterator iend = view.iend( *eIt );
  for( IntersectionIterator iIt = view.ibegin( *eIt ); iIt != iend; ++iIt )
  {
    const Intersection &intersection = *iIt;

    // perform intersection check

    checkIntersection( intersection, isCartesian );

    // check correctness of inside entity

    if( *eIt != intersection.inside() )
    {
      DUNE_THROW(Dune::GridError, "Error: Intersection's inside entity does not equal the entity the intersection iterator was started on.");
    }

    // check if conforming() method is compatible with static information on grid view
    if( GridViewType::conforming && !intersection.conforming() )
      DUNE_THROW( Dune::GridError, "Non-conforming intersection found in a conforming grid view (static)." );

    // check if conforming() method is compatible with dynamic information on grid view
    if( view.isConforming() && !intersection.conforming() )
      DUNE_THROW( Dune::GridError, "Non-conforming intersection found in a conforming grid view (dynamic)." );

    // check symmetry of 'has-intersection-with'-relation

    if( intersection.neighbor() && checkOutside )
    {
      const Entity outside = intersection.outside();

      bool insideFound = false;

      const IntersectionIterator oiend = view.iend( outside );
      for( IntersectionIterator oiit = view.ibegin( outside ); oiit != oiend; ++oiit )
      {
        if( !oiit->neighbor() || (oiit->outside() != intersection.inside()) )
          continue;

        insideFound = true;

        if( oiit->indexInInside() != intersection.indexInOutside() )
        {
          DUNE_THROW(Dune::GridError)
                    << "Error: symmetric intersection found, but with incorrect numbering.\n"
                    << "       (from inside: indexInOutside = " << intersection.indexInOutside()
                    << ", from outside: indexInInside = " << oiit->indexInInside() << ")." << std::endl;
        }
        if( oiit->indexInOutside() != intersection.indexInInside() )
        {
          DUNE_THROW(Dune::GridError)
                    << "Error: symmetric intersection found, but with incorrect numbering.\n"
                    << "       (from inside: indexInInside = " << intersection.indexInInside()
                    << ", from outside: indexInOutside = " << oiit->indexInOutside() << ")." << std::endl;
        }

        if( oiit->boundary() != intersection.boundary() )
        {
          DUNE_THROW(Dune::GridError)
                    << "Error: symmetric intersection found, but with incorrect boundary flag.\n"
                    << "       (from inside: boundary = " << intersection.boundary()
                    << ", from outside: boundary = " << oiit->boundary() << ")." << std::endl;
        }
      }

      if( !insideFound )
        DUNE_THROW( Dune::GridError, "Could not find inside entity through intersection iterator on outside entity." );
    }
    else if( !checkOutside )
    {
      static bool called = false;
      if( !called )
      {
        Dune::derr << "Warning: Skipping reverse intersection iterator check." << std::endl;
        called = true;
      }
    }

    std::set<typename GridViewType::IndexSet::IndexType> sub_indices;

    // check consistency of conforming intersection
    if( intersection.conforming() && intersection.neighbor() && !intersection.boundary() )
    {
      const Entity outside = intersection.outside();

      const int indexInInside = intersection.indexInInside();
      const int indexInOutside = intersection.indexInOutside();

      const auto& insideFace = eIt->template subEntity<1>( indexInInside );
      const auto& outsideFace = outside.template subEntity<1>( indexInOutside );

      const typename GridViewType::IndexSet &indexSet = view.indexSet();

      if( indexSet.index(insideFace) != indexSet.index(outsideFace) )
      {
        DUNE_THROW(Dune::GridError)
                  << "Error: Index of conforming intersection differs when "
                  << "obtained from inside and outside.\n"
                  << "       inside index = " << indexSet.index(insideFace)
                  << ", outside index = " << indexSet.index(outsideFace) << std::endl;
      }

      if( indexSet.subIndex( *eIt, indexInInside, 1 ) != indexSet.subIndex( outside, indexInOutside, 1 ) )
      {
        DUNE_THROW(Dune::GridError)
                  << "Error: SubIndex of conforming intersection differs when "
                  << "obtained from inside and outside.\n"
                  << "       inside sub-index = " << indexSet.subIndex( *eIt, indexInInside, 1 )
                  << ", outside sub-index = " << indexSet.subIndex( outside, indexInOutside, 1 ) << std::endl;
      }

      for (std::size_t cc = 0; cc != GridType::dimension; ++cc) {
        for (unsigned int subSubIndex = 0; subSubIndex != insideFace.subEntities(1 + cc); ++subSubIndex)
        {
          sub_indices.insert(indexSet.subIndex(insideFace, subSubIndex, 1 + cc));
        }
        for (unsigned int subSubIndex = 0; subSubIndex != outsideFace.subEntities(1 + cc); ++subSubIndex)
        {
          if (sub_indices.erase(indexSet.subIndex(outsideFace, subSubIndex, 1 + cc)) == 0)
          {
            std::cerr << "Error: SubSubIndices of codim " << 1+cc << " of conforming intersection differs when "
                      << "obtained from inside and outside." << std::endl;
            std::cerr << "       inside sub-sub-indices = ";
            for (unsigned int subSubIndex = 0; subSubIndex != insideFace.subEntities(1 + cc); ++subSubIndex)
              std::cerr << indexSet.subIndex( *eIt, indexInInside, 1 + cc);
            std::cerr << ", outside sub-sub-index = " << indexSet.subIndex(outsideFace, subSubIndex, 1 + cc) << std::endl;
            DUNE_THROW(Dune::GridError, "Intersection error");
          }
        }
        if (not sub_indices.empty()) {
          std::cerr << "Error: SubSubIndex of conforming intersection differs when "
                    << "obtained from inside and outside." << std::endl;
            std::cerr << "       inside sub-sub-indices = ";
            for (unsigned int subSubIndex = 0; subSubIndex != insideFace.subEntities(1 + cc); ++subSubIndex)
              std::cerr << indexSet.subIndex( *eIt, indexInInside, 1 + cc);
            std::cerr << ", outside sub-sub-indices = ";
            for (unsigned int subSubIndex = 0; subSubIndex != insideFace.subEntities(1 + cc); ++subSubIndex)
              std::cerr << indexSet.subIndex(outsideFace, subSubIndex, 1 + cc);
          DUNE_THROW(Dune::GridError, "Intersection error");
        }
      }

      const typename GridType::LocalIdSet &localIdSet = view.grid().localIdSet();
      if( localIdSet.subId( *eIt, indexInInside, 1 ) != localIdSet.subId( outside, indexInOutside, 1 ) )
      {
        DUNE_THROW(Dune::GridError)
                  << "Error: Local id of conforming intersection differs when "
                  << "obtained from inside and outside.\n"
                  << "       inside id = " << localIdSet.subId( *eIt, indexInInside, 1 )
                  << ", outside id = " << localIdSet.subId( outside, indexInOutside, 1 ) << std::endl;
      }

      const typename GridType::GlobalIdSet &globalIdSet = view.grid().globalIdSet();
      if( globalIdSet.subId( *eIt, indexInInside, 1 ) != globalIdSet.subId( outside, indexInOutside, 1 ) )
      {
        DUNE_THROW(Dune::GridError)
                  << "Error: Global id of conforming intersection differs when "
                  << "obtained from inside and outside.\n"
                  << "       inside id = " << globalIdSet.subId( *eIt, indexInInside, 1 )
                  << ", outside id = " << globalIdSet.subId( outside, indexInOutside, 1 ) << std::endl;
      }
    }

    // update variables for element checks

    if( intersection.boundary() )
      hasBoundaryIntersection = true;

    const Dune::QuadratureType::Enum qt = Dune::QuadratureType::GaussLegendre;
    const Dune::QuadratureRule< ctype, Intersection::mydimension > &quadrature
      = Dune::QuadratureRules< ctype, Intersection::mydimension >::rule( intersection.type(), 3, qt );
    for( std::size_t i = 0; i < quadrature.size(); ++i )
      sumNormal.axpy( quadrature[ i ].weight(), intersection.integrationOuterNormal( quadrature[ i ].position() ) );
  }

  // check implementation of hasBoundaryIntersections
  if( hasBoundaryIntersection != eIt->hasBoundaryIntersections() )
  {
    if( !hasBoundaryIntersection )
      DUNE_THROW(Dune::GridError, "Entity::hasBoundaryIntersections has a false positive.");
    else
      DUNE_THROW( Dune::GridError, "Entity::hasBoundaryIntersections has a false negative." );
  }

  // check whether integral over the outer normals is zero
  // note: This is wrong on curved surfaces (take, e.g., the upper half sphere).
  //       Therefore we only enforce this check on affine elements.
  //       This might be also wrong for network grids where intersections maybe overlapping
  //       Therefore we only enforce this check for dim==dimworld
  // note: The errorState variable will propagate the error as a warning for the cases
  //       where this check is not enforced
  if( (sumNormal.two_norm() > sqrt(std::numeric_limits< ctype >::epsilon())) && (eIt->partitionType() != Dune::GhostEntity) )
  {
    if( eIt->geometry().affine() && int(GridViewType::dimension) == int(GridViewType::dimensionworld))
      DUNE_THROW( Dune::GridError, "Integral over outer normals on affine entity is nonzero: " << sumNormal );
    ++errorState.sumNormalsNonZero;
  }

  // check assignment operator for IntersectionIterator
  checkIntersectionAssignment( view, eIt );
}

/** \brief Test both IntersectionIterators
 */
template <class GridViewType>
void checkViewIntersectionIterator(const GridViewType& view) {

  using namespace Dune;

  typedef typename GridViewType ::template Codim<0>::Iterator
  ElementIterator;
  ElementIterator eIt    = view.template begin<0>();
  ElementIterator eEndIt = view.template end<0>();

  CheckIntersectionIteratorErrorState errorState;
  for (; eIt!=eEndIt; ++eIt)
    checkIntersectionIterator( view, eIt, errorState );

  if( errorState.sumNormalsNonZero > 0 )
  {
    std :: cerr << "Warning: Integral over outer normals is not always zero."
                << std :: endl;
    std :: cerr << "         This behaviour may be correct for entities with"
                << " nonzero curvature, or in network grids." << std :: endl;
  }
}

template <class GridType>
void checkIntersectionIterator(const GridType& grid, bool skipLevelIntersectionTest = false) {

  using namespace Dune;

  // Loop over all levels
  if(skipLevelIntersectionTest)
  {
    std::cerr<<"WARNING: skip test of LevelIntersectionIterator! \n";
  }
  else
  {
    for (int i=0; i<=grid.maxLevel(); i++)
      checkViewIntersectionIterator(grid.levelGridView(i));
  }

  // test leaf intersection iterator
  {
    checkViewIntersectionIterator(grid.leafGridView());
  }

}

#endif // DUNE_GRID_TEST_CHECKINTERSECTIONIT_HH
