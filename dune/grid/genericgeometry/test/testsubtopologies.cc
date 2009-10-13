// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <cstring>

#include <dune/common/forloop.hh>

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
  static void apply ();
};



template< int codim >
void CheckCodim< codim > :: apply ()
{
  typedef Dune :: GenericGeometry :: Size< Topology, codim > NumSubs;
  Dune::ForLoop< CheckSub, 0, NumSubs::value-1 >::apply();
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

  Dune::ForLoop< CheckSubCodim, 0, SubTopology::dimension >::apply();
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

  Dune::ForLoop< CheckSubSub, 0, NumSubSubs::value - 1 >::apply();
}

template< int codim >
template< int i >
template< int subcodim >
template< int j >
void CheckCodim< codim > :: CheckSub< i > :: CheckSubCodim< subcodim >
:: CheckSubSub< j > :: apply ()
{
  /*
     typedef Dune :: GenericGeometry
     :: SubTopologyNumber< Topology, codim, i, subcodim, j > SubNumber;
   */
  typedef Dune :: GenericGeometry
  :: Size< Topology, codim + subcodim > NumSubs;

  const unsigned int number = Dune :: GenericGeometry
                              :: SubTopologyNumbering< Topology, codim, subcodim > :: number( i, j );
  const unsigned int genericNumber = Dune :: GenericGeometry
                                     :: GenericSubTopologyNumbering< Topology, codim, subcodim > :: number( i, j );

  bool error = false;
  error |= ((codim == 0) && (number != j));
  error |= ((subcodim == 0) && (number != i));
  error |= (number >= (unsigned int)NumSubs :: value);
  //error |= (number != SubNumber :: value);
  error |= (number != genericNumber);

  if( verbose || error )
  {
    std :: cerr << "SubTopologyNumber< " << codim << ", " << i
                << ", " << subcodim << ", " << j << " > = " << number
                << std :: endl;
  }
  if( error )
    ++errors;
}

// *********************************************
typedef Dune::GenericGeometry::Pyramid< Topology > YTopology;
typedef Dune::GenericGeometry::Pyramid< YTopology > XTopology;
template <class T,unsigned int i>
struct ForSubTest
{
  const static unsigned int codimension = XTopology::dimension-T::dimension;
  const static unsigned int size = Dune::GenericGeometry::Size< XTopology, codimension >::value;
  static void apply ( )
  {
    std::cout << codimension << " " << i << " -> " << T::name()
              << " ( " << size << " ) " << std::endl;
  }
};

int main ( int argc, char **argv )
{
  for( int i = 1; i < argc; ++i )
  {
    verbose |= (strcmp( argv[ i ], "-v" ) == 0);
  }

  std :: cerr << "Generic geometry type: " << Topology :: name() << std :: endl;

  Dune::ForLoop< CheckCodim, 0, Topology::dimension >::apply();

  std :: cerr << "Number of errors: " << errors << std :: endl;

  // ForSubTopology does not exist anymore
#if 0
  Dune :: GenericGeometry :: ForSubTopology< ForSubTest, XTopology >::apply( );
#endif

  return (errors > 0 ? 1 : 0);
}
