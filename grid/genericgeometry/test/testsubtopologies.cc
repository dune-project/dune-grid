// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/subtopologies.hh>

#ifndef GEOMETRYTYPE
#error "GEOMETRYTYPE must be one of 'simplex', 'cube', 'prism', 'pyramid'"
#endif

#ifndef DIMENSION
#error "DIMENSION not selected."
#endif

typedef Dune :: GenericGeometry :: Convert
< Dune :: GeometryType :: GEOMETRYTYPE, DIMENSION > :: type
Topology;

bool verbose = false;
unsigned int errors = 0;

typedef Dune :: GenericGeometry :: SubTopologyNumbering< Topology >
SubTopologyNumbering;



template< int codim >
struct CheckCodim
{
  template< int i >
  struct CheckSub;

  static void apply();
};

template< int codim >
template< int i >
struct CheckCodim< codim > :: CheckSub
{
  template< int subcodim >
  struct CheckSubCodim;

  static void apply ();
};

template< int codim >
template< int i >
template< int subcodim >
struct CheckCodim< codim > :: CheckSub< i > :: CheckSubCodim
{
  template< int j >
  struct CheckSubSub;

  static void apply ();
};

template< int codim >
template< int i >
template< int subcodim >
template< int j >
struct CheckCodim< codim > :: CheckSub< i > :: CheckSubCodim< subcodim >
:: CheckSubSub
{
  template< int subsubcodim >
  struct CheckSubSubCodim;

  static void apply ();
};

template< int codim >
template< int i >
template< int subcodim >
template< int j >
template< int subsubcodim >
struct CheckCodim< codim > :: CheckSub< i > :: CheckSubCodim< subcodim >
:: CheckSubSub< j > :: CheckSubSubCodim
{
  template< int k >
  struct CheckSubSubSub;

  static void apply ();
};

template< int codim >
template< int i >
template< int subcodim >
template< int j >
template< int subsubcodim >
template< int k >
struct CheckCodim< codim > :: CheckSub< i > :: CheckSubCodim< subcodim >
:: CheckSubSub< j > :: CheckSubSubCodim< subsubcodim > :: CheckSubSubSub
{
  static void apply ();
};




template< int codim >
void CheckCodim< codim > :: apply ()
{
  typedef Dune :: GenericGeometry :: Size< Topology, codim > NumSubs;
  Dune :: GenericGeometry :: ForLoop< CheckSub, 0, NumSubs :: value-1 > :: apply();
}

template< int codim >
template< int i >
void CheckCodim< codim > :: CheckSub< i > :: apply ()
{
  typedef typename Dune :: GenericGeometry
  :: SubTopology< Topology, codim, i > :: type SubTopology;

  if( verbose )
  {
    std :: cerr << "SubTopology< " << codim << " > " << i
                << ": type = " << SubTopology :: name()
                << std :: endl;
  }

  Dune :: GenericGeometry :: ForLoop
  < CheckSubCodim, 0, SubTopology :: dimension > :: apply();
}

template< int codim >
template< int i >
template< int subcodim >
void CheckCodim< codim > :: CheckSub< i > :: CheckSubCodim< subcodim > :: apply ()
{
  typedef typename Dune :: GenericGeometry
  :: SubTopology< Topology, codim, i > :: type SubTopology;
  typedef Dune :: GenericGeometry :: Size< SubTopology, subcodim > NumSubSubs;
  typedef Dune :: GenericGeometry :: SubTopologySize< Topology, codim, subcodim >
  SubTopologySize;

  bool error = (NumSubSubs :: value != SubTopologySize :: size( i ));
  if( verbose || error )
  {
    std :: cerr << "SubTopology< " << codim << " > " << i
                << ": size< " << subcodim << " > = " << NumSubSubs :: value
                << " (" << SubTopologySize :: size( i ) << ")"
                << std :: endl;
  }
  if( error )
    ++errors;

  Dune :: GenericGeometry :: ForLoop
  < CheckSubSub, 0, NumSubSubs :: value - 1 > :: apply();
}

template< int codim >
template< int i >
template< int subcodim >
template< int j >
void CheckCodim< codim > :: CheckSub< i > :: CheckSubCodim< subcodim >
:: CheckSubSub< j > :: apply ()
{
  typedef Dune :: GenericGeometry
  :: SubTopologyNumber< Topology, codim, i, subcodim, j > SubNumber;
  typedef Dune :: GenericGeometry
  :: Size< Topology, codim + subcodim > NumSubs;

  const unsigned int subEntity
    = SubTopologyNumbering :: template subEntity< codim, subcodim >( i, j );

  bool error = false;
  error |= ((codim == 0) && (SubNumber :: value != j));
  error |= ((subcodim == 0) && (SubNumber :: value != i));
  error |= ((unsigned int)SubNumber :: value >= (unsigned int)NumSubs :: value);
  error |= (SubNumber :: value != subEntity);

  if( verbose || error )
  {
    std :: cerr << "SubTopologyNumber< " << codim << ", " << i
                << ", " << subcodim << ", " << j << " > = "
                << SubNumber :: value << std :: endl;
  }
  if( error )
    ++errors;

  Dune :: GenericGeometry :: ForLoop
  < CheckSubSubCodim, 0, Topology :: dimension - codim - subcodim > :: apply();
}

template< int codim >
template< int i >
template< int subcodim >
template< int j >
template< int subsubcodim >
void CheckCodim< codim > :: CheckSub< i > :: CheckSubCodim< subcodim >
:: CheckSubSub< j > :: CheckSubSubCodim< subsubcodim > :: apply ()
{
  typedef Dune :: GenericGeometry
  :: SubTopologyNumber< Topology, codim, i, subcodim, j > SubNumber;
  typedef typename Dune :: GenericGeometry
  :: SubTopology< Topology, codim + subcodim, SubNumber :: value > :: type
  SubSubTopology;
  typedef Dune :: GenericGeometry
  :: Size< SubSubTopology, subsubcodim > NumSubSubSubs;

  Dune :: GenericGeometry :: ForLoop
  < CheckSubSubSub, 0, NumSubSubSubs :: value-1 > :: apply();
}

template< int codim >
template< int i >
template< int subcodim >
template< int j >
template< int subsubcodim >
template< int k >
void CheckCodim< codim > :: CheckSub< i > :: CheckSubCodim< subcodim >
:: CheckSubSub< j > :: CheckSubSubCodim< subsubcodim > :: CheckSubSubSub< k >
:: apply ()
{
  typedef typename Dune :: GenericGeometry
  :: SubTopology< Topology, codim, i > :: type SubTopology;
  typedef Dune :: GenericGeometry :: SubTopologyNumber
  < Topology, codim, i, subcodim, j > Nij;
  typedef Dune :: GenericGeometry :: SubTopologyNumber
  < SubTopology, subcodim, j, subsubcodim, k > Njk;

  typedef Dune :: GenericGeometry :: SubTopologyNumber
  < Topology, codim + subcodim, Nij :: value, subsubcodim, k > Nij_k;
  typedef Dune :: GenericGeometry :: SubTopologyNumber
  < Topology, codim, i, subcodim + subsubcodim, Njk :: value > Ni_jk;

  if( (unsigned int)Nij_k :: value != (unsigned int)Ni_jk :: value )
  {
    std :: cerr << "Error< " << codim << ", " << subcodim << ", "
                << subsubcodim << " >: "
                << "Sub( Sub( " << i << ", " << j << "), " << k << " ) = "
                << Nij_k :: value << " != " << Ni_jk :: value
                << " = Sub( " << i << ", Sub( " << j << ", " << k << " ) )"
                << std :: endl;
    ++errors;
  }
}



int main ( int argc, char **argv )
{
  for( int i = 1; i < argc; ++i )
  {
    verbose |= (strcmp( argv[ i ], "-v" ) == 0);
  }

  std :: cerr << "Generic geometry type: " << Topology :: name() << std :: endl;

  Dune :: GenericGeometry :: ForLoop< CheckCodim, 0, Topology :: dimension > :: apply();

  std :: cerr << "Number of errors: " << errors << std :: endl;
  return (errors > 0 ? 1 : 0);
}
