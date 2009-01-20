// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_MACRODATA_HH
#define DUNE_ALBERTA_MACRODATA_HH

#include <dune/grid/albertagrid/albertaheader.hh>
#include <dune/grid/albertagrid/misc.hh>

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

      /** \brief create a new macro data structure
       *
       *  A new macro data structure is created and put into insert mode.
       */
      void create ()
      {
        release();
        data_ = ALBERTA alloc_macro_data( dim, 4096, 4096 );
        vertexCount_ = elementCount_ = 0;
        elementCount_ = 0;
      }

      /** \brief compress macro data structure
       *
       *  Compress the macro data structure to its minimum size and leave
       *  insert mode.
       */
      void finalize ()
      {
        assert( (vertexCount_ >= 0) && (elementCount_ >= 0) );
        resizeVertices( vertexCount_ );
        resizeElements( elementCount_ );
        vertexCount_ = elementCount_ = -1;
      }

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
      int insertElement ( const int (&element)[ numVertices ] )
      {
        assert( elementCount_ >= 0 );
        if( elementCount_ >= data_->n_macro_elements )
          resizeElements( 2*elementCount_ );
        const int offset = elementCount_*numVertices;
        for( int i = 0; i < numVertices; ++i )
          data_->mel_vertices[ offset + i ] = element[ i ];
        return elementCount_++;
      }

      /** \brief insert vertex
       *
       *  Insert a vertex into the macro data structure. This may only be
       *  done in insert mode.
       */
      int insertVertex ( const GlobalVector &v )
      {
        assert( vertexCount_ >= 0 );
        if( vertexCount_ >= data_->n_total_vertices )
          resizeVertices( 2*vertexCount_ );
        for( int i = 0; i < dimWorld; ++i )
          data_->coords[ vertexCount_ ][ i ] = v[ i ];
        return vertexCount_++;
      }

      void read ( const std::string &filename, bool binary = false )
      {
        release();
        if( binary )
          data_ = ALBERTA read_macro_xdr( filename.c_str() );
        else
          data_ = ALBERTA read_macro( filename.c_str() );
      }

      bool write ( const std::string &filename, bool binary = false ) const
      {
        if( binary )
          return ALBERTA write_macro_data_xdr( data_, filename.c_str() );
        else
          return ALBERTA write_macro_data( data_, filename.c_str() );
      }

    private:
      void resizeElements ( const int newSize )
      {
        const int oldSize = data_->n_macro_elements;
        data_->n_macro_elements = newSize;
        data_->mel_vertices = MEM_REALLOC( data_->mel_vertices, oldSize*numVertices, newSize*numVertices, int );
        assert( data_->mel_vertices != NULL );
      }

      void resizeVertices ( const int newSize )
      {
        const int oldSize = data_->n_total_vertices;
        data_->n_total_vertices = newSize;
        data_->coords = MEM_REALLOC( data_->coords, oldSize, newSize, GlobalVector );
        assert( data_->coords != NULL );
      }
    };

  }

}

#endif
