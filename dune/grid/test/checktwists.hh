// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKTWISTS_HH
#define DUNE_GRID_TEST_CHECKTWISTS_HH

#include <dune/geometry/referenceelements.hh>

inline int applyTwist ( const int twist, const int i, const int numCorners )
{
  return (twist < 0 ? (2*numCorners + 1 - i + twist) : i + twist) % numCorners;
}

inline int inverseTwist ( const int twist, const int numCorners )
{
  return (twist <= 0 ? twist : numCorners - twist);
}

struct NoMapTwist
{
  int operator() ( const int twist ) const
  {
    return twist;
  }
};


struct TwistZero
{
  // overwrite twist with 0
  int operator() ( const int twist ) const
  {
    return 0;
  }
};

namespace detail {

  // for non-twist free grids this class needs to be specialized
  // to return the correct twist information
  template <class Intersection>
  struct GetTwist
  {
    static int twistInInside( const Intersection& intersection )
    {
      return 0;
    }

    static int twistInOutside( const Intersection& intersection )
    {
      return 0;
    }
  };

}

template< class Intersection, class MapTwist >
int checkTwistOnIntersection ( const Intersection &intersection, const MapTwist &mapTwist )
{
  typedef typename Intersection::Entity Entity;
  typedef typename Entity::Geometry Geometry;
  typedef typename Geometry::ctype ctype;

  const int dimension = Entity::dimension;

  typedef typename Intersection::LocalGeometry LocalGeometry;

  typedef Dune::FieldVector< typename Geometry::ctype, Geometry::coorddimension >
  WorldVector;
  typedef Dune::FieldVector< typename LocalGeometry::ctype, LocalGeometry::coorddimension >
  LocalVector;

  if( !intersection.neighbor() || intersection.boundary() || !intersection.conforming() )
    return 0;

  int errors = 0;
  const ctype tolerance = std::numeric_limits< ctype >::epsilon();

  const int tIn  = mapTwist( detail::GetTwist< Intersection >::twistInInside( intersection ) );
  const int tOut = mapTwist( detail::GetTwist< Intersection >::twistInOutside( intersection ) );

  const Entity entityIn = intersection.inside();
  const Entity entityOut = intersection.outside();

  const Geometry &geoIn = entityIn.geometry();
  const Geometry &geoOut = entityOut.geometry();

  const int nIn = intersection.indexInInside();
  const int nOut = intersection.indexInOutside();

  auto refIn = referenceElement( geoIn );
  auto refOut = referenceElement( geoOut );

  const int numCorners = refIn.size( nIn, 1, dimension );
  assert( refOut.size( nOut, 1, dimension ) == numCorners );
  for( int i = 0; i < numCorners; ++i )
  {
    const int iIn = applyTwist( inverseTwist( tIn, numCorners ), i, numCorners );
    const int iOut = applyTwist( inverseTwist( tOut, numCorners ), i, numCorners );

    WorldVector xIn = geoIn.corner( refIn.subEntity( nIn, 1, iIn, dimension ) );
    WorldVector xOut = geoOut.corner( refOut.subEntity( nOut, 1, iOut, dimension ) );

    if( (xIn - xOut).two_norm() < tolerance )
      continue;

    std::cout << "Error: corner( " << iIn << " ) = " << xIn
              << " != " << xOut << " = corner( " << iOut << " )."
              << std::endl;
    ++errors;
  }

  const LocalGeometry &lGeoIn = intersection.geometryInInside();
  const LocalGeometry &lGeoOut = intersection.geometryInOutside();

  bool twistInside = true ;
  bool twistOutside = true;
  for( int i = 0; i < numCorners; ++i )
  {
    const int gi = i;

    const int iIn = applyTwist( inverseTwist( tIn, numCorners ), i, numCorners );
    LocalVector xIn = refIn.position( refIn.subEntity( nIn, 1, iIn, dimension ), dimension );
    if( (xIn - lGeoIn.corner( gi )).two_norm() >= tolerance )
    {
      std::cout << "Error: twisted inside reference corner( " << iIn << " ) = " << xIn
                << " != " << lGeoIn.corner( gi ) << " = local corner( " << i << " )."
                << std::endl;
      twistInside = false ;
      ++errors;
    }

    const int iOut = applyTwist( inverseTwist( tOut, numCorners ), i, numCorners );
    LocalVector xOut = refOut.position( refOut.subEntity( nOut, 1, iOut, dimension ), dimension );
    if( (xOut - lGeoOut.corner( gi )).two_norm() >= tolerance )
    {
      std::cout << "Error: twisted outside reference corner( " << iOut << " ) = " << xOut
                << " != " << lGeoOut.corner( gi ) << " = local corner( " << i << " )."
                << std::endl;
      twistOutside = false;
      ++errors;
    }
  }

  // calculate inside twist
  if( ! twistInside )
  {
    for( int nTwist = numCorners-1; nTwist>= -numCorners; --nTwist )
    {
      twistInside = true ;
      for( int i = 0; i < numCorners; ++i )
      {
        const int gi = i;
        const int iIn = applyTwist( inverseTwist( nTwist, numCorners ), i, numCorners );
        LocalVector xIn = refIn.position( refIn.subEntity( nIn, 1, iIn, dimension ), dimension );
        if( (xIn - lGeoIn.corner( gi )).two_norm() >= tolerance )
        {
          twistInside = false ;
        }
      }

      if( twistInside )
      {
        std::cout << "\ninside " << nIn << "\n";
        std::cout << "twist " << tIn << " should be replaced by " << nTwist << "\n";
        break ;
      }
    }
  }

  // calculate outside twist
  if( ! twistOutside )
  {
    for( int nTwist = numCorners-1; nTwist>=-numCorners; --nTwist )
    {
      twistOutside = true ;
      for( int i = 0; i < numCorners; ++i )
      {
        const int gi = i;
        const int iOut = applyTwist( inverseTwist( nTwist, numCorners ), i, numCorners );
        LocalVector xOut = refOut.position( refOut.subEntity( nOut, 1, iOut, dimension ), dimension );
        if( (xOut - lGeoOut.corner( gi )).two_norm() >= tolerance )
        {
          twistOutside = false;
        }
      }

      if( twistOutside )
      {
        std::cout << "\noutside " << nOut << " (inside = " << nIn << ")\n";
        std::cout << "twist " << tOut << " should be replaced by " << nTwist << "\n";
        break ;
      }
    }
  }

  return errors;
}


template< class GridView, class MapTwist >
void checkTwists ( const GridView &gridView, const MapTwist &mapTwist )
{
  typedef typename GridView::template Codim< 0 >::Iterator Iterator;
  typedef typename GridView::IntersectionIterator IntersectionIterator;

  int errors = 0;

  const Iterator end = gridView.template end< 0 >();
  for( Iterator it = gridView.template begin< 0 >(); it != end; ++it )
  {
    const IntersectionIterator iend = gridView.iend( *it );
    for( IntersectionIterator iit = gridView.ibegin( *it ); iit != iend; ++iit )
      errors += checkTwistOnIntersection( *iit, mapTwist );
  }

  if( errors > 0 )
    std::cerr << "Error: Intersection twists contain errors." << std::endl;
}

template< class GridView >
void checkTwistFree ( const GridView &gridView )
{
  // twist free grids have twist 0 everywhere
  checkTwists( gridView, TwistZero() );
}

#endif // DUNE_GRID_TEST_CHECKTWISTS_HH
