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

bool verbose = false;
unsigned int errors = 0;

template< int codim, unsigned int i, int subcodim, unsigned int j >
struct CheckSubSubNumbering
{
  typedef Dune :: GenericGeometry
  :: SubEntityNumber< Geometry, codim, i, subcodim, j > SubNumber;
  typedef Dune :: GenericGeometry
  :: NumSubEntities< Geometry, codim + subcodim > NumSubs;

  static void check ()
  {
    CheckSubSubNumbering< codim, i, subcodim, j-1 > :: check();

    bool error = false;
    error |= ((codim == 0) && (SubNumber :: value != j));
    error |= ((subcodim == 0) && (SubNumber :: value != i));
    error |= ((unsigned int)SubNumber :: value >= (unsigned int)NumSubs :: value);

    if( verbose || error )
    {
      std :: cerr << "SubEntityNumber< " << codim << ", " << i
                  << ", " << subcodim << ", " << j << " > = "
                  << SubNumber :: value << std :: endl;
    }
    if( error )
      ++errors;
  }
};

template< int codim, unsigned int i, int subcodim >
struct CheckSubSubNumbering< codim, i, subcodim, 0 >
{
  typedef Dune :: GenericGeometry
  :: SubEntityNumber< Geometry, codim, i, subcodim, 0 > SubNumber;
  typedef Dune :: GenericGeometry
  :: NumSubEntities< Geometry, codim + subcodim > NumSubs;

  static void check ()
  {
    bool error = false;
    error |= ((codim == 0) && (SubNumber :: value != 0));
    error |= ((subcodim == 0) && (SubNumber :: value != i));
    error |= ((unsigned int)SubNumber :: value >= (unsigned int)NumSubs :: value);

    if( verbose || error )
    {
      std :: cerr << "SubEntityNumber< " << codim << ", " << i
                  << ", " << subcodim << ", " << 0 << " > = "
                  << SubNumber :: value << std :: endl;
    }
    if( error )
      ++errors;
  }
};


template< int codim, unsigned int i, int subcodim >
struct CheckNumSubSubHelper
{
  typedef typename Dune :: GenericGeometry
  :: SubGeometry< Geometry, codim, i > :: type SubGeo;
  typedef Dune :: GenericGeometry :: NumSubEntities< SubGeo, subcodim >
  NumSubSubs;

  static void check ()
  {
    CheckNumSubSubHelper< codim, i, subcodim-1 > :: check();

    if( verbose )
    {
      std :: cerr << "SubEntity< " << codim << " > " << i
                  << ": size< " << subcodim << " > = " << NumSubSubs :: value
                  << std :: endl;
    }

    CheckSubSubNumbering< codim, i, subcodim, NumSubSubs :: value - 1 > :: check();
  }
};

template< int codim, unsigned int i >
struct CheckNumSubSubHelper< codim, i, 0 >
{
  typedef typename Dune :: GenericGeometry
  :: SubGeometry< Geometry, codim, i > :: type SubGeo;
  typedef Dune :: GenericGeometry :: NumSubEntities< SubGeo, 0 >
  NumSubSubs;

  static void check ()
  {
    if( verbose )
    {
      std :: cerr << "SubEntity< " << codim << " > " << i
                  << ": size< 0 > = " << NumSubSubs :: value
                  << std :: endl;
    }

    CheckSubSubNumbering< codim, i, 0, NumSubSubs :: value - 1 > :: check();
  }
};


template< int codim, unsigned int i >
struct CheckNumSubSubs
{
  typedef typename Dune :: GenericGeometry
  :: SubGeometry< Geometry, codim, i > :: type SubGeo;

  static void check()
  {
    CheckNumSubSubs< codim, i-1 > :: check();

    if( verbose )
    {
      std :: cerr << "SubEntity< " << codim << " > " << i
                  << ": type = " << SubGeo :: name()
                  << std :: endl;
    }

    CheckNumSubSubHelper< codim, i, Geometry :: dimension - codim > :: check();
  }
};

template< int codim >
struct CheckNumSubSubs< codim, 0 >
{
  typedef typename Dune :: GenericGeometry
  :: SubGeometry< Geometry, codim, 0 > :: type SubGeo;

  static void check()
  {
    if( verbose )
    {
      std :: cerr << "SubEntity< " << codim << " > 0"
                  << ": type = " << SubGeo :: name()
                  << std :: endl;
    }

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

int main ( int argc, char **argv )
{
  for( int i = 1; i < argc; ++i )
  {
    verbose |= (strcmp( argv[ i ], "-v" ) == 0);
  }

  std :: cerr << "Generic geometry type: " << Geometry :: name() << std :: endl;

  CheckCodim< Geometry :: dimension > :: check();

  std :: cerr << "Number of errors: " << errors << std :: endl;
  return (errors > 0 ? 1 : 0);
}
