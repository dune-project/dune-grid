// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
int applyTwist ( const int twist, const int i, const int numCorners )
{
  return (twist < 0 ? (2*numCorners + 1 - i + twist) : i + twist) % numCorners;
}

int inverseTwist ( const int twist, const int numCorners )
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


template< class Intersection, class MapTwist >
int checkTwistOnIntersection ( const Intersection &intersection, const MapTwist &mapTwist )
{
  typedef typename Intersection::ctype ctype;

  const int dimension = Intersection::dimension;

  typedef typename Intersection::EntityPointer EntityPointer;
  typedef typename Intersection::Entity Entity;
  typedef typename Entity::Geometry Geometry;

  typedef typename Intersection::LocalGeometry LocalGeometry;

  typedef Dune::ReferenceElement< ctype, dimension > ReferenceElement;
  typedef Dune::ReferenceElements< ctype, dimension > ReferenceElements;

  typedef FieldVector< typename Geometry::ctype, Geometry::coorddimension >
  WorldVector;
  typedef FieldVector< typename LocalGeometry::ctype, LocalGeometry::coorddimension >
  LocalVector;

  if( !intersection.neighbor() || !intersection.conforming() )
    return 0;

  int errors = 0;

  const int tIn = mapTwist( intersection.twistInSelf() );
  const int tOut = mapTwist( intersection.twistInNeighbor() );

  const EntityPointer ptrIn = intersection.inside();
  const EntityPointer ptrOut = intersection.outside();
  const Entity &entityIn = *ptrIn;
  const Entity &entityOut = *ptrOut;
  const Geometry &geoIn = entityIn.geometry();
  const Geometry &geoOut = entityOut.geometry();

  typedef GenericGeometry::MapNumberingProvider< dimension > Numbering;
  const unsigned int tidIn = GenericGeometry::topologyId( entityIn.type() );
  const int nIn = Numbering::template generic2dune< 1 >( tidIn, intersection.indexInInside() );
  const unsigned int tidOut = GenericGeometry::topologyId( entityOut.type() );
  const int nOut = Numbering::template generic2dune< 1 >( tidOut, intersection.indexInOutside() );

  const ReferenceElement &refIn = ReferenceElements::general( geoIn.type() );
  const ReferenceElement &refOut = ReferenceElements::general( geoOut.type() );

  const int numCorners = refIn.size( nIn, 1, dimension );
  assert( refOut.size( nOut, 1, dimension ) == numCorners );
  for( int i = 0; i < numCorners; ++i )
  {
    const int iIn = applyTwist( inverseTwist( tIn, numCorners ), i, numCorners );
    const int iOut = applyTwist( inverseTwist( tOut, numCorners ), i, numCorners );

    WorldVector xIn = geoIn.corner( refIn.subEntity( nIn, 1, iIn, dimension ) );
    WorldVector xOut = geoOut.corner( refOut.subEntity( nOut, 1, iOut, dimension ) );

    if( (xIn - xOut).two_norm() < 1e-12 )
      continue;

    std::cout << "Error: corner( " << iIn << " ) = " << xIn
              << " != " << xOut << " = corner( " << iOut << " )."
              << std::endl;
    ++errors;
  }

  const LocalGeometry &lGeoIn = intersection.geometryInInside();
  const LocalGeometry &lGeoOut = intersection.geometryInOutside();

  for( int i = 0; i < numCorners; ++i )
  {
    const int tid = Dune::GenericGeometry::topologyId( lGeoIn.type() );
    const int gi = Dune::GenericGeometry::MapNumberingProvider< dimension-1 >::template dune2generic< dimension-1 >( tid, i );
    //assert( lGeoIn[ i ] == lGeoIn.corner( gi ) );
    //assert( lGeoOut[ i ] == lGeoOut.corner( gi ) );

    const int iIn = applyTwist( inverseTwist( tIn, numCorners ), i, numCorners );
    LocalVector xIn = refIn.position( refIn.subEntity( nIn, 1, iIn, dimension ), dimension );
    if( (xIn - lGeoIn.corner( gi )).two_norm() >= 1e-12 )
    {
      std::cout << "Error: twisted reference corner( " << iIn << " ) = " << xIn
                << " != " << lGeoIn.corner( gi ) << " = local corner( " << i << " )."
                << std::endl;
      ++errors;
    }

    const int iOut = applyTwist( inverseTwist( tOut, numCorners ), i, numCorners );
    LocalVector xOut = refOut.position( refOut.subEntity( nOut, 1, iOut, dimension ), dimension );
    if( (xOut - lGeoOut.corner( gi )).two_norm() >= 1e-12 )
    {
      std::cout << "Error: twisted reference corner( " << iOut << " ) = " << xOut
                << " != " << lGeoOut.corner( gi ) << " = local corner( " << i << " )."
                << std::endl;
      ++errors;
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
      errors += checkTwistOnIntersection( gridView.grid().getRealIntersection( *iit ), mapTwist );
  }

  if( errors > 0 )
    std::cerr << "Error: Intersection twists contain errors." << std::endl;
}
