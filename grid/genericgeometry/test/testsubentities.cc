// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/genericgeometry/geometrytypes.hh>
#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/subentities.hh>

#ifndef GEOMETRYTYPE
#error "GEOMETRYTYPE must be one of 'simplex', 'cube', 'prism', 'pyramid'"
#endif

#ifndef DIMENSION
#error "DIMENSION not selected."
#endif

typedef Dune :: GenericGeometry :: Convert
< Dune :: GeometryType :: GEOMETRYTYPE, DIMENSION > :: type
Geometry;


template< int codim, unsigned int i, int subcodim >
struct CheckNumSubSubHelper
{
  typedef Dune :: GenericGeometry :: NumSubSubEntities
  < Geometry, codim, i, subcodim >
  NumSubSubs;

  typedef typename Dune :: GenericGeometry
  :: SubGeometry< Geometry, codim, i > :: type SubGeo;
  typedef Dune :: GenericGeometry :: NumSubEntities< SubGeo, subcodim >
  NumSubGeoSubs;

  static void check ()
  {
    std :: cerr << "SubEntity< " << codim << " > " << i << ": "
                << "size = " << NumSubSubs :: value
                << ", (" << NumSubGeoSubs :: value << " )"
                << std :: endl;
  }
};

template< int codim, unsigned int i >
struct CheckNumSubSubHelper< codim, i, 0 >
{
  typedef Dune :: GenericGeometry :: NumSubSubEntities
  < Geometry, codim, i, 0 >
  NumSubSubs;

  typedef typename Dune :: GenericGeometry
  :: SubGeometry< Geometry, codim, i > :: type SubGeo;
  typedef Dune :: GenericGeometry :: NumSubEntities< SubGeo, 0 >
  NumSubGeoSubs;

  static void check ()
  {
    std :: cerr << "SubEntity< " << codim << " > " << i << ": "
                << "size = " << NumSubSubs :: value
                << ", (" << NumSubGeoSubs :: value << " )"
                << std :: endl;
  }
};


template< int codim, unsigned int i >
struct CheckNumSubSubs
{
  static void check()
  {
    CheckNumSubSubs< codim, i-1 > :: check();
    CheckNumSubSubHelper< codim, i, Geometry :: dimension - codim > :: check();
  }
};

template< int codim >
struct CheckNumSubSubs< codim, 0 >
{
  static void check()
  {
    CheckNumSubSubHelper< codim, 0, Geometry :: dimension - codim > :: check();
  }
};


template< int codim >
struct CheckCodim
{
  typedef Dune :: GenericGeometry :: NumSubEntities< Geometry, codim > NumSubs;

  static void check()
  {
    CheckCodim< codim-1 > :: check();
    CheckNumSubSubs< codim, NumSubs :: value-1 > :: check();
  }
};

template<>
struct CheckCodim< 0 >
{
  typedef Dune :: GenericGeometry :: NumSubEntities< Geometry, 0 > NumSubs;

  static void check()
  {
    CheckNumSubSubs< 0, NumSubs :: value-1 > :: check();
  }
};

int main ()
{
  std :: cout << "Generic geometry type: " << Geometry :: name() << std :: endl;

  CheckCodim< Geometry :: dimension > :: check();
}
