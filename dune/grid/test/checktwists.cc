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

#if NEW_SUBENTITY_NUMBERING
  typedef Dune::GenericReferenceElement< ctype, dimension > ReferenceElement;
  typedef Dune::GenericReferenceElements< ctype, dimension > ReferenceElements;
#else
  typedef Dune::ReferenceElement< ctype, dimension > ReferenceElement;
  typedef Dune::ReferenceElements< ctype, dimension > ReferenceElements;
#endif

  typedef Dune::FieldVector< typename Geometry::ctype, Geometry::coorddimension >
  WorldVector;
  typedef Dune::FieldVector< typename LocalGeometry::ctype, LocalGeometry::coorddimension >
  LocalVector;

  if( !intersection.neighbor() || intersection.boundary() || !intersection.conforming() )
    return 0;

  int errors = 0;

#if NEW_SUBENTITY_NUMBERING
  const int tIn = mapTwist( intersection.twistInInside() );
  const int tOut = mapTwist( intersection.twistInOutside() );
#else
  const int tIn = mapTwist( intersection.twistInSelf() );
  const int tOut = mapTwist( intersection.twistInNeighbor() );
#endif

  const EntityPointer ptrIn = intersection.inside();
  const EntityPointer ptrOut = intersection.outside();
  const Entity &entityIn = *ptrIn;
  const Entity &entityOut = *ptrOut;
  const Geometry &geoIn = entityIn.geometry();
  const Geometry &geoOut = entityOut.geometry();

#if NEW_SUBENTITY_NUMBERING
  const int nIn = intersection.indexInInside();
  const int nOut = intersection.indexInOutside();
#else
  typedef GenericGeometry::MapNumberingProvider< dimension > Numbering;
  const unsigned int tidIn = GenericGeometry::topologyId( entityIn.type() );
  const int nIn = Numbering::template generic2dune< 1 >( tidIn, intersection.indexInInside() );
  const unsigned int tidOut = GenericGeometry::topologyId( entityOut.type() );
  const int nOut = Numbering::template generic2dune< 1 >( tidOut, intersection.indexInOutside() );
#endif

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

  bool twistInside = true ;
  bool twistOutside = true;
  for( int i = 0; i < numCorners; ++i )
  {
#if NEW_SUBENTITY_NUMBERING
    const int gi = i;
#else
    const int tid = Dune::GenericGeometry::topologyId( lGeoIn.type() );
    const int gi = Dune::GenericGeometry::MapNumberingProvider< dimension-1 >::template dune2generic< dimension-1 >( tid, i );
    //assert( lGeoIn[ i ] == lGeoIn.corner( gi ) );
    //assert( lGeoOut[ i ] == lGeoOut.corner( gi ) );
#endif

    const int iIn = applyTwist( inverseTwist( tIn, numCorners ), i, numCorners );
    LocalVector xIn = refIn.position( refIn.subEntity( nIn, 1, iIn, dimension ), dimension );
    if( (xIn - lGeoIn.corner( gi )).two_norm() >= 1e-12 )
    {
      std::cout << "Error: twisted inside reference corner( " << iIn << " ) = " << xIn
                << " != " << lGeoIn.corner( gi ) << " = local corner( " << i << " )."
                << std::endl;
      twistInside = false ;
      ++errors;
    }

    const int iOut = applyTwist( inverseTwist( tOut, numCorners ), i, numCorners );
    LocalVector xOut = refOut.position( refOut.subEntity( nOut, 1, iOut, dimension ), dimension );
    if( (xOut - lGeoOut.corner( gi )).two_norm() >= 1e-12 )
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
#if NEW_SUBENTITY_NUMBERING
        const int gi = i;
#else
        const int tid = Dune::GenericGeometry::topologyId( lGeoIn.type() );
        const int gi = Dune::GenericGeometry::MapNumberingProvider< dimension-1 >::template dune2generic< dimension-1 >( tid, i );
#endif
        const int iIn = applyTwist( inverseTwist( nTwist, numCorners ), i, numCorners );
        LocalVector xIn = refIn.position( refIn.subEntity( nIn, 1, iIn, dimension ), dimension );
        if( (xIn - lGeoIn.corner( gi )).two_norm() >= 1e-12 )
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
#if NEW_SUBENTITY_NUMBERING
        const int gi = i;
#else
        const int tid = Dune::GenericGeometry::topologyId( lGeoIn.type() );
        const int gi = Dune::GenericGeometry::MapNumberingProvider< dimension-1 >::template dune2generic< dimension-1 >( tid, i );
#endif
        const int iOut = applyTwist( inverseTwist( nTwist, numCorners ), i, numCorners );
        LocalVector xOut = refOut.position( refOut.subEntity( nOut, 1, iOut, dimension ), dimension );
        if( (xOut - lGeoOut.corner( gi )).two_norm() >= 1e-12 )
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
      errors += checkTwistOnIntersection( gridView.grid().getRealIntersection( *iit ), mapTwist );
  }

  if( errors > 0 )
    std::cerr << "Error: Intersection twists contain errors." << std::endl;
}
