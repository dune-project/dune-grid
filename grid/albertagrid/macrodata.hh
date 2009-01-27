// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_MACRODATA_HH
#define DUNE_ALBERTA_MACRODATA_HH

#include <dune/common/fvector.hh>

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/albertaheader.hh>
#include <dune/grid/albertagrid/referencetopo.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    template< int dim >
    class MacroData
    {
      typedef MacroData< dim > This;

      typedef ALBERTA MACRO_DATA Data;

      static const int dimension = dim;
      static const int numVertices = NumSubEntities< dimension, dimension >::value;
      static const int numEdges = NumSubEntities< dimension, dimension-1 >::value;

      static const int initialSize = 4096;

    public:
      typedef int ElementId[ numVertices ];

      static const int supportPeriodicity = (DUNE_ALBERTA_VERSION >= 0x201);

    private:
      Data *data_;
      int vertexCount_;
      int elementCount_;

    public:
      MacroData ()
        : data_( NULL ),
          vertexCount_( -1 ),
          elementCount_( -1 )
      {}

      operator Data * () const
      {
        return data_;
      }

      int vertexCount () const
      {
        return (vertexCount_ < 0 ? data_->n_total_vertices : vertexCount_);
      }

      int elementCount () const
      {
        return (elementCount_ < 0 ? data_->n_macro_elements : elementCount_);
      }

      ElementId &element ( int i ) const;
      GlobalVector &vertex ( int i ) const;
      int &neighbor ( int element, int i ) const;
      BoundaryId &boundaryId ( int element, int i ) const;

      /** \brief create a new macro data structure
       *
       *  A new macro data structure is created and put into insert mode.
       */
      void create ();

      /** \brief compress macro data structure
       *
       *  Compress the macro data structure to its minimum size and leave
       *  insert mode.
       *
       *  \note This method may always be called. It does nothing outside of
       *        insert mode.
       */
      void finalize ();

      /** \brief mark the longest edge of all elements as refinement edges
       *
       *  This is a postprocessing step and should be done after finalizing
       *  the triangulation.
       *
       *  \note Thouth it is possible to call markLongestEdge in insert mode,
       *        you must make sure that all required vertices have been set.
       */
      void markLongestEdge ();

      /** \brief release the macro data structure */
      void release ()
      {
        if( data_ != NULL )
        {
          ALBERTA free_macro_data( data_ );
          data_ = NULL;
        }
        vertexCount_ = elementCount_ = -1;
      }

      /** \brief insert element
       *
       *  Insert an element into the macro data structure. This may only be
       *  done in insert mode.
       */
      int insertElement ( const ElementId &id )
      {
        assert( elementCount_ >= 0 );
        if( elementCount_ >= data_->n_macro_elements )
          resizeElements( 2*elementCount_ );

        ElementId &e = element( elementCount_ );
        for( int i = 0; i < numVertices; ++i )
        {
          e[ i ] = id[ i ];
          boundaryId( elementCount_, i ) = InteriorBoundary;
        }

        return elementCount_++;
      }

      /** \brief insert vertex
       *
       *  Insert a vertex into the macro data structure. This may only be
       *  done in insert mode.
       */
      int insertVertex ( const GlobalVector &coords )
      {
        assert( vertexCount_ >= 0 );
        if( vertexCount_ >= data_->n_total_vertices )
          resizeVertices( 2*vertexCount_ );
        copy( coords, vertex( vertexCount_ ) );
        return vertexCount_++;
      }

      /** \brief insert vertex
       *
       *  Insert a vertex into the macro data structure. This may only be
       *  done in insert mode.
       */
      int insertVertex ( const FieldVector< Real, dimWorld > &coords )
      {
        assert( vertexCount_ >= 0 );
        if( vertexCount_ >= data_->n_total_vertices )
          resizeVertices( 2*vertexCount_ );
        copy( coords, vertex( vertexCount_ ) );
        return vertexCount_++;
      }

      void insertWallTrafo ( const GlobalMatrix &m, const GlobalVector &t );
      void insertWallTrafo ( const FieldMatrix< Real, dimWorld, dimWorld > &matrix,
                             const FieldVector< Real, dimWorld > &shift );

      void read ( const std::string &filename, bool binary = false );

      bool write ( const std::string &filename, bool binary = false ) const
      {
        if( binary )
          return ALBERTA write_macro_data_xdr( data_, filename.c_str() );
        else
          return ALBERTA write_macro_data( data_, filename.c_str() );
      }

    private:
      Real edgeLength ( const ElementId &e, int edge ) const;
      int longestEdge ( const ElementId &e ) const;

      template< class Vector >
      void copy ( const Vector &x, GlobalVector &y )
      {
        for( int i = 0; i < dimWorld; ++i )
          y[ i ] = x[ i ];
      }

#if DUNE_ALBERTA_VERSION >= 0x200
      template< class Type >
      void rotate ( Type *array, int i, int shift );
#else
      template< class Type >
      void rotate ( Type (*array)[ numVertices ], int i, int shift );
#endif

      void resizeElements ( const int newSize );

      void resizeVertices ( const int newSize )
      {
        const int oldSize = data_->n_total_vertices;
        data_->n_total_vertices = newSize;
        data_->coords = memReAlloc< GlobalVector >( data_->coords, oldSize, newSize );
        assert( data_->coords != NULL );
      }
    };


#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    inline typename MacroData< dim >::ElementId &
    MacroData< dim >::element ( int i ) const
    {
      assert( (i >= 0) && (i < data_->n_macro_elements) );
      const int offset = i * numVertices;
      return *reinterpret_cast< ElementId * >( data_->mel_vertices + offset );
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x200

#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    inline typename MacroData< dim >::ElementId &
    MacroData< dim >::element ( int i ) const
    {
      assert( (i >= 0) && (i < data_->n_macro_elements) );
      return data_->mel_vertices[ i ];
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x200


    template< int dim >
    inline GlobalVector &MacroData< dim >::vertex ( int i ) const
    {
      assert( (i >= 0) && (i < data_->n_total_vertices) );
      return data_->coords[ i ];
    }


#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    inline int &MacroData< dim >::neighbor ( int element, int i ) const
    {
      assert( (element >= 0) && (element < data_->n_macro_elements) );
      assert( (i >= 0) && (i < numVertices) );
      return data_->neigh[ element*numVertices + i ];
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x200

#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    inline int &MacroData< dim >::neighbor ( int element, int i ) const
    {
      assert( (element >= 0) && (element < data_->n_macro_elements) );
      assert( (i >= 0) && (i < numVertices) );
      return data_->neigh[ element ][ i ];
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x200


#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    inline BoundaryId &MacroData< dim >::boundaryId ( int element, int i ) const
    {
      assert( (element >= 0) && (element < data_->n_macro_elements) );
      assert( (i >= 0) && (i < numVertices) );
      return data_->boundary[ element*numVertices + i ];
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x200

#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    inline BoundaryId &MacroData< dim >::boundaryId ( int element, int i ) const
    {
      assert( (element >= 0) && (element < data_->n_macro_elements) );
      assert( (i >= 0) && (i < numVertices) );
      return data_->boundary[ element ][ i ];
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x200


#if DUNE_ALBERTA_VERSION >= 0x201
    template< int dim >
    inline void MacroData< dim >::create ()
    {
      release();
      data_ = ALBERTA alloc_macro_data( dim, initialSize, initialSize );
      data_->boundary = memAlloc< BoundaryId >( initialSize*numVertices );
      vertexCount_ = elementCount_ = 0;
      elementCount_ = 0;
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x201

#if DUNE_ALBERTA_VERSION == 0x200
    template< int dim >
    inline void MacroData< dim >::create ()
    {
      release();
      data_ = ALBERTA alloc_macro_data( dim, initialSize, initialSize, 0 );
      data_->boundary = memAlloc< BoundaryId >( initialSize*numVertices );
      vertexCount_ = elementCount_ = 0;
      elementCount_ = 0;
    }
#endif // #if DUNE_ALBERTA_VERSION == 0x200

#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    inline void MacroData< dim >::create ()
    {
      dune_static_assert( dimension == dimGrid,
                          "Wrong grid dimension used for ALBERTA 1.2." );
      release();
      data_ = ALBERTA alloc_macro_data( initialSize, initialSize, 0 );
      data_->boundary = memAlloc< BoundaryId[ numVertices ] >( initialSize );
      vertexCount_ = elementCount_ = 0;
      elementCount_ = 0;
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x200


#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    inline void MacroData< dim >::finalize ()
    {
      if( (vertexCount_ >= 0) && (elementCount_ >= 0) )
      {
        resizeVertices( vertexCount_ );
        resizeElements( elementCount_ );
        ALBERTA compute_neigh_fast( data_ );

        // assign default boundary id (if none is assigned)
        for( int element = 0; element < elementCount_; ++element )
        {
          for( int i = 0; i < numVertices; ++i )
          {
            BoundaryId &id = boundaryId( element, i );
            if( neighbor( element, i ) >= 0 )
            {
              assert( id == InteriorBoundary );
              id = InteriorBoundary;
            }
            else
              id = (id == InteriorBoundary ? DirichletBoundary : id);
          }
        }

        vertexCount_ = elementCount_ = -1;
      }
      assert( (vertexCount_ < 0) && (elementCount_ < 0) );
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x200

#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    inline void MacroData< dim >::finalize ()
    {
      if( (vertexCount_ >= 0) && (elementCount_ >= 0) )
      {
        resizeVertices( vertexCount_ );
        resizeElements( elementCount_ );

        std::cerr << "Warning: GridFactory for ALBERTA 1.2 does not support "
                  << "boundary ids, yet." << std::endl << std::endl;
        memFree( data_->boundary, elementCount_ );
        data_->boundary = NULL;

        vertexCount_ = elementCount_ = -1;
      }
      assert( (vertexCount_ < 0) && (elementCount_ < 0) );
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x200


    template< int dim >
    inline void MacroData< dim >::markLongestEdge ()
    {
      assert( data_ != NULL );
      if( dimension == 1 )
        return;

      const int count = elementCount();
      for( int i = 0; i < count; ++i )
      {
        const int refEdge = RefinementEdge< dimension >::value;
        const int edge = longestEdge( element( i ) );
        if( edge == refEdge )
          continue;
        std::cout << "Longest edge: " << edge << std::endl;

        // shift vertices such that the refinement edge is the longest edge
        const int shift = edge + (numVertices - refEdge);

        // rotate necessary fields
        rotate( data_->mel_vertices, i, shift );

        // correct opposite vertices
#if DUNE_ALBERTA_VERSION >= 0x201
        if( data_->opp_vertex != NULL )
        {
          assert( data_->neigh != NULL );
          for( int j = 0; j < numVertices; ++j )
          {
            const int nb = data_->neigh[ i*numVertices + j ];
            if( nb < 0 )
              continue;
            const int ov = data_->opp_vertex[ i*numVertices + j ];
            assert( data_->neigh[ nb*numVertices + ov ] == i );
            assert( data_->opp_vertex[ nb*numVertices + ov ] == j );
            data_->opp_vertex[ nb*numVertices + ov ] = (j+shift) % numVertices;
          }
          rotate( data_->opp_vertex, i, shift );
        }
#endif

        // correct neighbors and boundaries
        rotate( data_->neigh, i, shift );
        rotate( data_->boundary, i, shift );
      }
    }


#if DUNE_ALBERTA_VERSION >= 0x201
    template< int dim >
    inline void MacroData< dim >
    ::insertWallTrafo ( const GlobalMatrix &matrix, const GlobalVector &shift )
    {
      typedef ALBERTA AFF_TRAFO AffineTrafo;

      int &count = data_->n_wall_trafos;
      AffineTrafo *&array = data_->wall_trafos;

      // resize wall trafo array
      array = memReAlloc< AffineTrafo >( array, count, count+1 );
      assert( data_->wall_trafos != NULL );

      // copy matrix and shift
      for( int i = 0; i < dimWorld; ++i )
        copy( matrix[ i ], array[ count ].M[ i ] );
      copy( shift, array[ count ].t );
      ++count;
    }

    template< int dim >
    inline void MacroData< dim >
    ::insertWallTrafo ( const FieldMatrix< Real, dimWorld, dimWorld > &matrix,
                        const FieldVector< Real, dimWorld > &shift )
    {
      typedef ALBERTA AFF_TRAFO AffineTrafo;

      int &count = data_->n_wall_trafos;
      AffineTrafo *&array = data_->wall_trafos;

      // resize wall trafo array
      array = memReAlloc< AffineTrafo >( array, count, count+1 );
      assert( data_->wall_trafos != NULL );

      // copy matrix and shift
      for( int i = 0; i < dimWorld; ++i )
        copy( matrix[ i ], array[ count ].M[ i ] );
      copy( shift, array[ count ].t );
      ++count;
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x201

#if DUNE_ALBERTA_VERSION <= 0x200
    template< int dim >
    inline void MacroData< dim >
    ::insertWallTrafo ( const GlobalMatrix &m, const GlobalVector &t )
    {
      DUNE_THROW( AlbertaError,
                  "Periodic grids are only supported in ALBERTA 2.1 or higher." );
    }

    template< int dim >
    inline void MacroData< dim >
    ::insertWallTrafo ( const FieldMatrix< Real, dimWorld, dimWorld > &matrix,
                        const FieldVector< Real, dimWorld > &shift )
    {
      DUNE_THROW( AlbertaError,
                  "Periodic grids are only supported in ALBERTA 2.1 or higher." );
    }
#endif // #if DUNE_ALBERTA_VERSION <= 0x200


#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    inline void MacroData< dim >::read ( const std::string &filename, bool binary )
    {
      release();
      if( binary )
        data_ = ALBERTA read_macro_xdr( filename.c_str() );
      else
        data_ = ALBERTA read_macro( filename.c_str() );
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x200

#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    inline void MacroData< dim >::read ( const std::string &filename, bool binary )
    {
      release();
      DUNE_THROW( NotImplemented, "In ALBERTA 1.2, macro data cannot be read." );
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x200


    template< int dim >
    inline Real MacroData< dim >::edgeLength ( const ElementId &e, int edge ) const
    {
      const int i = ALBERTA AlbertHelp::MapVertices< 1, 3 >::mapVertices( edge, 0 );
      assert( (vertexCount_ < 0) || (e[ i ] < vertexCount_) );
      const GlobalVector &x = vertex( e[ i ] );

      const int j = ALBERTA AlbertHelp::MapVertices< 1, 3 >::mapVertices( edge, 1 );
      assert( (vertexCount_ < 0) || (e[ j ] < vertexCount_) );
      const GlobalVector &y = vertex( e[ j ] );

      Real sum = (y[ 0 ] - x[ 0 ]) * (y[ 0 ] - x[ 0 ]);
      for( int i = 1; i < dimWorld; ++i )
        sum += (y[ i ] - x[ i ]) * (y[ i ] - x[ i ]);
      return sqrt( sum );
    }


    template< int dim >
    inline int MacroData< dim >::longestEdge ( const ElementId &e ) const
    {
      int maxEdge = 0;
      Real maxLength = edgeLength( e, 0 );
      for( int i = 1; i < numEdges; ++i )
      {
        const Real length = edgeLength( e, i );
        if( length <= maxLength )
          continue;
        maxEdge = i;
        maxLength = length;
      }
      return maxEdge;
    }


#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    template< class Type >
    inline void MacroData< dim >::rotate ( Type *array, int i, int shift )
    {
      assert( (i >= 0) && (i < data_->n_macro_elements) );
      if( array == NULL )
        return;

      const int offset = i*numVertices;
      Type old[ numVertices ];
      for( int j = 0; j < numVertices; ++j )
        old[ j ] = array[ offset + j ];
      for( int j = 0; j < numVertices; ++j )
        array[ offset + j ] = old[ (j+shift) % numVertices ];
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x200

#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    template< class Type >
    inline void MacroData< dim >::rotate ( Type (*array)[ numVertices ], int i, int shift )
    {
      assert( (i >= 0) && (i < data_->n_macro_elements) );
      if( array == NULL )
        return;

      Type old[ numVertices ];
      for( int j = 0; j < numVertices; ++j )
        old[ j ] = array[ i ][ j ];
      for( int j = 0; j < numVertices; ++j )
        array[ i ][ j ] = old[ (j+shift) % numVertices ];
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x200


#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    inline void MacroData< dim >::resizeElements ( const int newSize )
    {
      const int oldSize = data_->n_macro_elements;
      data_->n_macro_elements = newSize;
      data_->mel_vertices = memReAlloc( data_->mel_vertices, oldSize*numVertices, newSize*numVertices );
      data_->boundary = memReAlloc( data_->boundary, oldSize*numVertices, newSize*numVertices );
      assert( data_->mel_vertices != NULL );
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x200

#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    inline void MacroData< dim >::resizeElements ( const int newSize )
    {
      const int oldSize = data_->n_macro_elements;
      data_->n_macro_elements = newSize;
      data_->mel_vertices = memReAlloc( data_->mel_vertices, oldSize, newSize );
      assert( data_->mel_vertices != NULL );
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x200

  }

}

#endif // #if HAVE_ALBERTA

#endif
