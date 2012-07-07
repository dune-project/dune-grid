// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_SGRIDINTERNAL_GENERIC2DUNE_HH
#define DUNE_GRID_SGRIDINTERNAL_GENERIC2DUNE_HH

/** \file
 *  \brief Renumber from the Dune subentity numbering to the SGrid-internal one
 *
 * Once upon a time SGrid used a certain system to number its element subentities.
 * That system worked, but the code that implemented it wasn't documented.
 * So in order to understand that hidden numbering system, a second system was
 * implemented additionally, and it was voted that the new system should be the
 * official Dune numbering system.  For a transitional period code was implemented
 * that computed from one system to the other.  All code was eventually updated
 * to use the new system, except for SGrid.  And that is why the transformation
 * code from Dune numbering to SGrid numbering is still here.
 */

#include <dune/common/static_assert.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/genericgeometry/subtopologies.hh>

namespace Dune
{

  namespace SGridInternal
  {

    // MapNumbering
    // ------------

    template< class Topology >
    struct MapNumbering;


    struct MapNumberingIdentical
    {
      template< unsigned int codim >
      static unsigned int generic2dune ( unsigned int i )
      {
        return i;
      }
    };

    // MapNumbering for Point
    template<>
    struct MapNumbering< GenericGeometry::Point >
      : public MapNumberingIdentical
    {};

    // MapNumbering for Line
    template<>
    struct MapNumbering< GenericGeometry::Prism< GenericGeometry::Point > >
      : public MapNumberingIdentical
    {};

    template<>
    struct MapNumbering< GenericGeometry::Pyramid< GenericGeometry::Point > >
      : public MapNumberingIdentical
    {};

    // MapNumbering for Triangle
    struct MapNumberingTriangle
    {
      template< unsigned int codim >
      static unsigned int generic2dune ( unsigned int i )
      {
        return (codim == 1 ? 2 - i : i);
      }
    };

    template<>
    struct MapNumbering< GenericGeometry::Pyramid< GenericGeometry::Pyramid< GenericGeometry::Point > > >
      : public MapNumberingTriangle
    {};

    template<>
    struct MapNumbering< GenericGeometry::Pyramid< GenericGeometry::Prism< GenericGeometry::Point > > >
      : public MapNumberingTriangle
    {};


    // MapNumbering for Quadrilateral
    template<>
    struct MapNumbering< GenericGeometry::Prism< GenericGeometry::Pyramid< GenericGeometry::Point > > >
      : public MapNumberingIdentical
    {};

    template<>
    struct MapNumbering< GenericGeometry::Prism< GenericGeometry::Prism< GenericGeometry::Point > > >
      : public MapNumberingIdentical
    {};

    // MapNumbering for Tetrahedron
    struct MapNumberingTetrahedron
    {
      template< unsigned int codim >
      static unsigned int generic2dune ( unsigned int i )
      {
        static unsigned int edge[ 6 ] = { 0, 2, 1, 3, 4, 5 };
        return (codim == 1 ? 3 - i : (codim == 2 ? edge[ i ] : i));
      }
    };

    template<>
    struct MapNumbering< GenericGeometry::Pyramid< GenericGeometry::Pyramid< GenericGeometry::Pyramid< GenericGeometry::Point > > > >
      : public MapNumberingTetrahedron
    {};

    template<>
    struct MapNumbering< GenericGeometry::Pyramid< GenericGeometry::Pyramid< GenericGeometry::Prism< GenericGeometry::Point > > > >
      : public MapNumberingTetrahedron
    {};

    // MapNumbering for Cube
    struct MapNumberingCube
    {
      template< unsigned int codim >
      static unsigned int generic2dune ( unsigned int i )
      {
        static unsigned int edge[ 12 ] = { 0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11 };
        return (codim == 2 ? edge[ i ] : i);
      }
    };

    template<>
    struct MapNumbering< GenericGeometry::Prism< GenericGeometry::Prism< GenericGeometry::Pyramid< GenericGeometry::Point > > > >
      : public MapNumberingCube
    {};

    template<>
    struct MapNumbering< GenericGeometry::Prism< GenericGeometry::Prism< GenericGeometry::Prism< GenericGeometry::Point > > > >
      : public MapNumberingCube
    {};

    // MapNumbering for 4D-Cube
    struct MapNumbering4DCube
    {
      template< unsigned int codim >
      static unsigned int generic2dune ( unsigned int i )
      {
        static unsigned int codim2[ 24 ] =
        { 0, 1, 2, 3, 4, 5, 12, 13, 6, 7, 14, 15,
          8, 9, 16, 17, 20, 21, 10, 11, 18, 19, 22, 23 };
        static unsigned int codim3[ 32 ] =
        { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 24, 25,
          18, 19, 26, 27, 12, 13, 14, 15, 20, 21, 28, 29, 22, 23, 30, 31 };
        if (codim == 2)
          return codim2[i];
        else if (codim == 3)
          return codim3[i];
        else
          return i;
      }
    };

    template<>
    struct MapNumbering< GenericGeometry::Prism< GenericGeometry::Prism< GenericGeometry::Prism< GenericGeometry::Pyramid< GenericGeometry::Point > > > > >
      : public MapNumbering4DCube
    {};

    template<>
    struct MapNumbering< GenericGeometry::Prism< GenericGeometry::Prism< GenericGeometry::Prism< GenericGeometry::Prism< GenericGeometry::Point > > > > >
      : public MapNumbering4DCube
    {};

    // MapNumbering for Pyramid
    struct MapNumberingPyramid
    {
      template< unsigned int codim >
      static unsigned int generic2dune ( unsigned int i )
      {
        static unsigned int vertex[ 5 ] = { 0, 1, 3, 2, 4 };
        static unsigned int edge[ 8 ] = { 3, 1, 0, 2, 4, 5, 7, 6 };
        static unsigned int face[ 5 ] = { 0, 4, 2, 1, 3 };

        if( codim == 3 )
          return vertex[ i ];
        else if( codim == 2 )
          return edge[ i ];
        else if( codim == 1 )
          return face[ i ];
        else
          return i;
      }
    };

    template<>
    struct MapNumbering< GenericGeometry::Pyramid< GenericGeometry::Prism< GenericGeometry::Pyramid< GenericGeometry::Point > > > >
      : public MapNumberingPyramid
    {};

    template<>
    struct MapNumbering< GenericGeometry::Pyramid< GenericGeometry::Prism< GenericGeometry::Prism< GenericGeometry::Point > > > >
      : public MapNumberingPyramid
    {};

    // MapNumbering for Prism
    struct MapNumberingPrism
    {
      template< unsigned int codim >
      static unsigned int generic2dune ( unsigned int i )
      {
        static unsigned int edge[ 9 ] = { 3, 4, 5, 0, 2, 1, 6, 8, 7 };
        static unsigned int face[ 5 ] = { 1, 3, 2, 0, 4 };

        if( codim == 2 )
          return edge[ i ];
        else if( codim == 1 )
          return face[ i ];
        else
          return i;
      }
    };

    template<>
    struct MapNumbering< GenericGeometry::Prism< GenericGeometry::Pyramid< GenericGeometry::Pyramid< GenericGeometry::Point > > > >
      : public MapNumberingPrism
    {};

    template<>
    struct MapNumbering< GenericGeometry::Prism< GenericGeometry::Pyramid< GenericGeometry::Prism< GenericGeometry::Point > > > >
      : public MapNumberingPrism
    {};



    // MapNumberingProvider
    // --------------------

    template< unsigned int dim >
    struct MapNumberingProvider
    {
      static const unsigned int dimension = dim;
      static const unsigned int numTopologies = (1 << dimension);

    private:
      template< int i >
      struct Builder;

      typedef std :: vector< unsigned int > Map;

      Map generic2dune_[ numTopologies ][ dimension+1 ];

      MapNumberingProvider ()
      {
        ForLoop< Builder, 0, (1 << dim)-1 >::apply( generic2dune_ );
      }

      static const MapNumberingProvider &instance ()
      {
        static MapNumberingProvider inst;
        return inst;
      }

    public:
      static unsigned int
      generic2dune ( unsigned int topologyId, unsigned int i, unsigned int codim )
      {
        assert( (topologyId < numTopologies) && (codim <= dimension) );
        return instance().generic2dune_[ topologyId ][ codim ][ i ];
      }
    };


    template< unsigned int dim >
    template< int topologyId >
    struct MapNumberingProvider< dim >::Builder
    {
      typedef typename GenericGeometry::Topology< topologyId, dimension >::type Topology;
      typedef SGridInternal::MapNumbering< Topology > MapNumbering;

      template< int codim >
      struct Codim;

      static void apply ( Map (&generic2dune)[ numTopologies ][ dimension+1 ] )
      {
        ForLoop< Codim, 0, dimension >::apply( generic2dune[ topologyId ] );
      }
    };

    template< unsigned int dim >
    template< int i >
    template< int codim >
    struct MapNumberingProvider< dim >::Builder< i >::Codim
    {
      static void apply ( Map (&generic2dune)[ dimension+1 ] )
      {
        const unsigned int size = GenericGeometry::Size< Topology, codim >::value;

        Map &g2d = generic2dune[ codim ];
        g2d.resize( size );
        for( unsigned int j = 0; j < size; ++j )
          g2d[ j ] = MapNumbering::template generic2dune< codim >( j );
      }
    };

  }

}

#endif // #ifndef DUNE_GRID_SGRIDINTERNAL_GENERIC2DUNE_HH
