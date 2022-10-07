// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_MISC_HH
#define DUNE_ALBERTA_MISC_HH

#include <cassert>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/albertagrid/albertaheader.hh>

#if HAVE_ALBERTA

// should the coordinates be cached in a vector (required for ALBERTA 2.0)?
#ifndef DUNE_ALBERTA_CACHE_COORDINATES
#define DUNE_ALBERTA_CACHE_COORDINATES 1
#endif

namespace Dune
{

  // Exceptions
  // ----------

  class AlbertaError
    : public Exception
  {};

  class AlbertaIOError
    : public IOError
  {};



  namespace Alberta
  {

    // Import Types
    // ------------

    static const int dimWorld = DIM_OF_WORLD;

    typedef ALBERTA REAL Real;
    typedef ALBERTA REAL_B LocalVector; // in barycentric coordinates
    typedef ALBERTA REAL_D GlobalVector;
    typedef ALBERTA REAL_DD GlobalMatrix;
    typedef ALBERTA AFF_TRAFO AffineTransformation;
    typedef ALBERTA MESH Mesh;
    typedef ALBERTA EL Element;

    static const int meshRefined = MESH_REFINED;
    static const int meshCoarsened = MESH_COARSENED;

    static const int InteriorBoundary = INTERIOR;
    static const int DirichletBoundary = DIRICHLET;
    typedef ALBERTA BNDRY_TYPE BoundaryId;

    typedef U_CHAR ElementType;

    typedef ALBERTA FE_SPACE DofSpace;



    // Memory Manipulation Functions
    // -----------------------------

    template< class Data >
    inline Data *memAlloc ( size_t size )
    {
      return MEM_ALLOC( size, Data );
    }

    template< class Data >
    inline Data *memCAlloc ( size_t size )
    {
      return MEM_CALLOC( size, Data );
    }

    template< class Data >
    inline Data *memReAlloc ( Data *ptr, size_t oldSize, size_t newSize )
    {
      return MEM_REALLOC( ptr, oldSize, newSize, Data );
    }

    template< class Data >
    inline void memFree ( Data *ptr, size_t size )
    {
      return MEM_FREE( ptr, size, Data );
    }



    // GlobalSpace
    // -----------

    class GlobalSpace
    {
      typedef GlobalSpace This;

    public:
      typedef GlobalMatrix Matrix;
      typedef GlobalVector Vector;

    private:
      Matrix identityMatrix_;
      Vector nullVector_;

      GlobalSpace ()
      {
        for( int i = 0; i < dimWorld; ++i )
        {
          for( int j = 0; j < dimWorld; ++j )
            identityMatrix_[ i ][ j ] = Real( 0 );
          identityMatrix_[ i ][ i ] = Real( 1 );
          nullVector_[ i ] = Real( 0 );
        }
      }

      static This &instance ()
      {
        static This theInstance;
        return theInstance;
      }

    public:
      static const Matrix &identityMatrix ()
      {
        return instance().identityMatrix_;
      }

      static const Vector &nullVector ()
      {
        return instance().nullVector_;
      }
    };



    // NumSubEntities
    // --------------

    template< int dim, int codim >
    struct NumSubEntities;

    template< int dim >
    struct NumSubEntities< dim, 0 >
    {
      static const int value = 1;
    };

    template< int dim >
    struct NumSubEntities< dim, dim >
    {
      static const int value = dim+1;
    };

    template<>
    struct NumSubEntities< 0, 0 >
    {
      static const int value = 1;
    };

    template<>
    struct NumSubEntities< 2, 1 >
    {
      static const int value = 3;
    };

    template<>
    struct NumSubEntities< 3, 1 >
    {
      static const int value = 4;
    };

    template<>
    struct NumSubEntities< 3, 2 >
    {
      static const int value = 6;
    };



    // CodimType
    // ---------

    template< int dim, int codim >
    struct CodimType;

    template< int dim >
    struct CodimType< dim, 0 >
    {
      static const int value = CENTER;
    };

    template< int dim >
    struct CodimType< dim, dim >
    {
      static const int value = VERTEX;
    };

    template<>
    struct CodimType< 2, 1 >
    {
      static const int value = EDGE;
    };

    template<>
    struct CodimType< 3, 1 >
    {
      static const int value = FACE;
    };

    template<>
    struct CodimType< 3, 2 >
    {
      static const int value = EDGE;
    };



    // FillFlags
    // ---------

    template< int dim >
    struct FillFlags
    {
      typedef ALBERTA FLAGS Flags;

      static const Flags nothing = FILL_NOTHING;

      static const Flags coords = FILL_COORDS;

      static const Flags neighbor = FILL_NEIGH;

      static const Flags orientation = (dim == 3 ? FILL_ORIENTATION : FILL_NOTHING);

      static const Flags projection = FILL_PROJECTION;

      static const Flags elementType = FILL_NOTHING;

      static const Flags boundaryId = FILL_MACRO_WALLS;

      static const Flags nonPeriodic = FILL_NON_PERIODIC;

      static const Flags all = coords | neighbor | boundaryId | nonPeriodic
                               | orientation | projection | elementType;

      static const Flags standardWithCoords = all & ~nonPeriodic & ~projection;

#if DUNE_ALBERTA_CACHE_COORDINATES
      static const Flags standard = standardWithCoords & ~coords;
#else
      static const Flags standard = standardWithCoords;
#endif
    };



    // RefinementEdge
    // --------------

    template< int dim >
    struct RefinementEdge
    {
      static const int value = 0;
    };

    template<>
    struct RefinementEdge< 2 >
    {
      static const int value = 2;
    };



    // Dune2AlbertaNumbering
    // ---------------------

    template< int dim, int codim >
    struct Dune2AlbertaNumbering
    {
      static int apply ( const int i )
      {
        assert( (i >= 0) && (i < NumSubEntities< dim, codim >::value) );
        return i;
      }
    };

    template<>
    struct Dune2AlbertaNumbering< 3, 2 >
    {
      static const int numSubEntities = NumSubEntities< 3, 2 >::value;

      static int apply ( const int i )
      {
        assert( (i >= 0) && (i < numSubEntities) );
        static int dune2alberta[ numSubEntities ] = { 0, 3, 1, 2, 4, 5 };
        return dune2alberta[ i ];
      }
    };



    // Generic2AlbertaNumbering
    // ------------------------

    template< int dim, int codim >
    struct Generic2AlbertaNumbering
    {
      static int apply ( const int i )
      {
        assert( (i >= 0) && (i < NumSubEntities< dim, codim >::value) );
        return i;
      }
    };

    template< int dim >
    struct Generic2AlbertaNumbering< dim, 1 >
    {
      static int apply ( const int i )
      {
        assert( (i >= 0) && (i < NumSubEntities< dim, 1 >::value) );
        return dim - i;
      }
    };

    template<>
    struct Generic2AlbertaNumbering< 1, 1 >
    {
      static int apply ( const int i )
      {
        assert( (i >= 0) && (i < NumSubEntities< 1, 1 >::value) );
        return i;
      }
    };

    template<>
    struct Generic2AlbertaNumbering< 3, 2 >
    {
      static const int numSubEntities = NumSubEntities< 3, 2 >::value;

      static int apply ( const int i )
      {
        assert( (i >= 0) && (i < numSubEntities) );
        static int generic2alberta[ numSubEntities ] = { 0, 1, 3, 2, 4, 5 };
        return generic2alberta[ i ];
      }
    };



    // NumberingMap
    // ------------

    template< int dim, template< int, int > class Numbering = Generic2AlbertaNumbering >
    class NumberingMap
    {
      typedef NumberingMap< dim, Numbering > This;

      template< int codim >
      struct Initialize;

      int *dune2alberta_[ dim+1 ];
      int *alberta2dune_[ dim+1 ];
      int numSubEntities_[ dim+1 ];

      NumberingMap ( const This & );
      This &operator= ( const This & );

    public:
      NumberingMap ()
      {
        Hybrid::forEach( std::make_index_sequence< dim+1 >{}, [ & ]( auto i ){ Initialize< i >::apply( *this ); } );
      }

      ~NumberingMap ()
      {
        for( int codim = 0; codim <= dim; ++codim )
        {
          delete[]( dune2alberta_[ codim ] );
          delete[]( alberta2dune_[ codim ] );
        }
      }

      int dune2alberta ( int codim, int i ) const
      {
        assert( (codim >= 0) && (codim <= dim) );
        assert( (i >= 0) && (i < numSubEntities( codim )) );
        return dune2alberta_[ codim ][ i ];
      }

      int alberta2dune ( int codim, int i ) const
      {
        assert( (codim >= 0) && (codim <= dim) );
        assert( (i >= 0) && (i < numSubEntities( codim )) );
        return alberta2dune_[ codim ][ i ];
      }

      int numSubEntities ( int codim ) const
      {
        assert( (codim >= 0) && (codim <= dim) );
        return numSubEntities_[ codim ];
      }
    };



    // NumberingMap::Initialize
    // ------------------------

    template< int dim, template< int, int > class Numbering >
    template< int codim >
    struct NumberingMap< dim, Numbering >::Initialize
    {
      static const int numSubEntities = NumSubEntities< dim, codim >::value;

      static void apply ( NumberingMap< dim, Numbering > &map )
      {
        map.numSubEntities_[ codim ] = numSubEntities;
        map.dune2alberta_[ codim ] = new int[ numSubEntities ];
        map.alberta2dune_[ codim ] = new int[ numSubEntities ];

        for( int i = 0; i < numSubEntities; ++i )
        {
          const int j = Numbering< dim, codim >::apply( i );
          map.dune2alberta_[ codim ][ i ] = j;
          map.alberta2dune_[ codim ][ j ] = i;
        }
      }
    };



    // MapVertices
    // -----------

    template< int dim, int codim >
    struct MapVertices;

    template< int dim >
    struct MapVertices< dim, 0 >
    {
      static int apply ( int subEntity, int vertex )
      {
        assert( subEntity == 0 );
        assert( (vertex >= 0) && (vertex <= NumSubEntities< dim, dim >::value) );
        return vertex;
      }
    };

    template<>
    struct MapVertices< 2, 1 >
    {
      static int apply ( int subEntity, int vertex )
      {
        assert( (subEntity >= 0) && (subEntity < 3) );
        assert( (vertex >= 0) && (vertex < 2) );
        //static const int map[ 3 ][ 2 ] = { {1,2}, {2,0}, {0,1} };
        static const int map[ 3 ][ 2 ] = { {1,2}, {0,2}, {0,1} };
        return map[ subEntity ][ vertex ];
      }
    };

    template<>
    struct MapVertices< 3, 1 >
    {
      static int apply ( int subEntity, int vertex )
      {
        assert( (subEntity >= 0) && (subEntity < 4) );
        assert( (vertex >= 0) && (vertex < 3) );
        //static const int map[ 4 ][ 3 ] = { {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1} };
        static const int map[ 4 ][ 3 ] = { {1,2,3}, {0,2,3}, {0,1,3}, {0,1,2} };
        return map[ subEntity ][ vertex ];
      }
    };

    template<>
    struct MapVertices< 3, 2 >
    {
      static int apply ( int subEntity, int vertex )
      {
        assert( (subEntity >= 0) && (subEntity < 6) );
        assert( (vertex >= 0) && (vertex < 2) );
        static const int map[ 6 ][ 2 ] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
        return map[ subEntity ][ vertex ];
      }
    };

    template< int dim >
    struct MapVertices< dim, dim >
    {
      static int apply ( int subEntity, int vertex )
      {
        assert( (subEntity >= 0) && (subEntity < NumSubEntities< dim, 1 >::value) );
        assert( vertex == 0 );
        return subEntity;
      }
    };



    // Twist
    // -----

    // ******************************************************************
    // Meaning of the twist (same as in ALU)
    // -------------------------------------
    //
    // Consider a fixed ordering of the vertices v_1, ... v_n of a face
    // (here, we assume their indices to be increasing). Denote by k the
    // local number of a vertex v within the element and by t the twist.
    // Then, v = v_j, where j is computed by the following formula:
    //
    //        / (2n + 1 - k + t) % n, if t < 0
    //   j = <
    //        \ (k + t) % n,          if t >= 0
    //
    //  Note: We use the order of the 0-th vertex dof to assign the twist.
    //        This is ok for two reasons:
    //        - ALBERTA preserves the relative order of the dofs during
    //          dof compression.
    //        - ALBERTA enforces the first vertex dof admin to be periodic.
    // ******************************************************************

    template< int dim, int subdim >
    struct Twist
    {
      static const int numSubEntities = NumSubEntities< dim, dim-subdim >::value;

      static const int minTwist = 0;
      static const int maxTwist = 0;

      static int twist ( [[maybe_unused]] const Element *element,
                         [[maybe_unused]] int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < numSubEntities) );
        return 0;
      }
    };

    template< int dim >
    struct Twist< dim, 1 >
    {
      static const int numSubEntities = NumSubEntities< dim, dim-1 >::value;

      static const int minTwist = 0;
      static const int maxTwist = 1;

      static int twist ( const Element *element, int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < numSubEntities) );
        const int numVertices = NumSubEntities< 1, 1 >::value;
        int dof[ numVertices ];
        for( int i = 0; i < numVertices; ++i )
        {
          const int j = MapVertices< dim, dim-1 >::apply( subEntity, i );
          dof[ i ] = element->dof[ j ][ 0 ];
        }
        return (dof[ 0 ] < dof[ 1 ] ? 0 : 1);
      }
    };


    template<>
    struct Twist< 1, 1 >
    {
      static const int minTwist = 0;
      static const int maxTwist = 0;

      static int twist ( [[maybe_unused]] const Element *element,
                         [[maybe_unused]] int subEntity )
      {
        assert( subEntity == 0 );
        return 0;
      }
    };


    template< int dim >
    struct Twist< dim, 2 >
    {
      static const int numSubEntities = NumSubEntities< dim, dim-2 >::value;

      static const int minTwist = -3;
      static const int maxTwist = 2;

      static int twist ( const Element *element, int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < numSubEntities) );
        const int numVertices = NumSubEntities< 2, 2 >::value;
        int dof[ numVertices ];
        for( int i = 0; i < numVertices; ++i )
        {
          const int j = MapVertices< dim, dim-2 >::apply( subEntity, i );
          dof[ i ] = element->dof[ j ][ 0 ];
        }

        const int twist[ 8 ] = { -2, 1, 666, -1, 2, 666, -3, 0 };
        const int k = int( dof[ 0 ] < dof[ 1 ] )
                      | (int( dof[ 0 ] < dof[ 2 ] ) << 1)
                      | (int( dof[ 1 ] < dof[ 2 ] ) << 2);
        assert( twist[ k ] != 666 );
        return twist[ k ];
      }
    };


    template<>
    struct Twist< 2, 2 >
    {
      static const int minTwist = 0;
      static const int maxTwist = 0;

      static int twist ( [[maybe_unused]] const Element *element,
                         [[maybe_unused]] int subEntity )
      {
        assert( subEntity == 0 );
        return 0;
      }
    };



    template< int dim >
    inline int applyTwist ( int twist, int i )
    {
      const int numCorners = NumSubEntities< dim, dim >::value;
      return (twist < 0 ? (2*numCorners + 1 - i + twist) : i + twist) % numCorners;
    }

    template< int dim >
    inline int applyInverseTwist ( int twist, int i )
    {
      const int numCorners = NumSubEntities< dim, dim >::value;
      return (twist < 0 ? (2*numCorners + 1 - i + twist) : numCorners + i - twist) % numCorners;
    }

  }

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_MISC_HH
