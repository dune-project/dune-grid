// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_MISC_HH
#define DUNE_ALBERTA_MISC_HH

#include <cassert>

#include <dune/common/exceptions.hh>

#include <dune/grid/genericgeometry/codimtable.hh>
#include <dune/grid/genericgeometry/misc.hh>

#include <dune/grid/albertagrid/albertaheader.hh>

#if HAVE_ALBERTA

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

    using GenericGeometry::ForLoop;
    using GenericGeometry::CodimTable;

    static const int dimWorld = DIM_OF_WORLD;
#if DUNE_ALBERTA_VERSION < 0x200
    static const int dimGrid = DIM;
#endif // #if DUNE_ALBERTA_VERSION < 0x200

    typedef ALBERTA REAL Real;
    typedef ALBERTA REAL_D GlobalVector;

    typedef ALBERTA MESH Mesh;
    typedef ALBERTA MACRO_EL MacroElement;
    typedef ALBERTA EL Element;

    static const int InteriorBoundary = INTERIOR;
    static const int DirichletBoundary = DIRICHLET;
#if DUNE_ALBERTA_VERSION >= 0x201
    typedef ALBERTA BNDRY_TYPE BoundaryId;
#else
    typedef S_CHAR BoundaryId;
#endif
#if DUNE_ALBERTA_VERSION < 0x200
    typedef ALBERTA BOUNDARY Boundary;
#endif

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

#if (DUNE_ALBERTA_VERSION >= 0x200) || (DIM == 3)
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
#endif



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
      static int apply ( int i )
      {
        assert( (i >= 0) && (i < NumSubEntities< dim, codim >::value) );
        return i;
      }
    };

    template<>
    struct Dune2AlbertaNumbering< 3, 2 >
    {
      static const int numSubEntities = NumSubEntities< 3, 2 >::value;

      static int apply ( int i )
      {
        assert( (i >= 0) && (i < numSubEntities) );
        static int dune2alberta[ numSubEntities ] = { 0, 3, 1, 2, 4, 5 };
        return dune2alberta[ i ];
      }
    };



    // NumberingMap
    // ------------

    template< int dim >
    class NumberingMap
    {
      typedef NumberingMap< dim > This;

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
        ForLoop< Initialize, 0, dim >::apply( *this );
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



    // NumberingMap::Setup
    // -------------------

    template< int dim >
    template< int codim >
    struct NumberingMap< dim >::Initialize
    {
      static const int numSubEntities = NumSubEntities< dim, codim >::value;

      static void apply ( NumberingMap< dim > &map )
      {
        map.numSubEntities_[ codim ] = numSubEntities;
        map.dune2alberta_[ codim ] = new int[ numSubEntities ];
        map.alberta2dune_[ codim ] = new int[ numSubEntities ];

        for( int i = 0; i < numSubEntities; ++i )
        {
          const int j = Dune2AlbertaNumbering< dim, codim >::apply( i );
          map.dune2alberta_[ codim ][ i ] = j;
          map.alberta2dune_[ codim ][ j ] = i;
        }
      }
    };

  }

}

#endif // #if HAVE_ALBERTA

#endif
