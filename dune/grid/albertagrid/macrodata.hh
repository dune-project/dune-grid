// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_MACRODATA_HH
#define DUNE_ALBERTA_MACRODATA_HH

/** \file
 *  \author Martin Nolte
 *  \brief  provides a wrapper for ALBERTA's macro_data structure
 */

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/algebra.hh>
#include <dune/grid/albertagrid/albertaheader.hh>

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
      template< int >
      struct Library;

      template< int > friend struct InstantiateMacroDataLibrary;

    public:
      typedef int ElementId[ numVertices ];

      static const int supportPeriodicity = 1;

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
       *  \note Though it is possible to call markLongestEdge in insert mode,
       *        you must make sure that all required vertices have been set.
       */
      void markLongestEdge ();

      /** \brief set the orientation of all elements
       *
       *  This is a postprocessing step and should be done after finalizing
       *  the triangulation.
       *
       *  \note Though it is possible to call setOrientation in insert mode,
       *        you must make sure that all required vertices have been set.
       */
      void setOrientation ( const Real orientation );

      /** \brief check the neighbor information
       *
       *  This method allows the verification of neighbor information in a
       *  finalized (and possibly postprecessed) macro triangulation.
       *
       *  \note On unfinalized macro triangulations there is no neighbor
       *        information. Hence this check will succeed in this case.
       *
       *  \returns true, if all generated neighbor information is correct.
       */
      bool checkNeighbors () const;

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
      int insertElement ( const ElementId &id );

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

      void checkCycles ();

      void read ( const std::string &filename, bool binary = false );

      bool write ( const std::string &filename, bool binary = false ) const
      {
        if( binary )
          return ALBERTA write_macro_data_xdr( data_, filename.c_str() );
        else
          return ALBERTA write_macro_data( data_, filename.c_str() );
      }

    private:
      template< class Vector >
      void copy ( const Vector &x, GlobalVector &y )
      {
        for( int i = 0; i < dimWorld; ++i )
          y[ i ] = x[ i ];
      }

      void resizeElements ( const int newSize );

      void resizeVertices ( const int newSize )
      {
        const int oldSize = data_->n_total_vertices;
        data_->n_total_vertices = newSize;
        data_->coords = memReAlloc< GlobalVector >( data_->coords, oldSize, newSize );
        assert( (data_->coords != NULL) || (newSize == 0) );
      }

    private:
      Data *data_;
      int vertexCount_;
      int elementCount_;
    };



    // MacroData::Library
    // ------------------

    template< int dim >
    template< int >
    struct MacroData< dim >::Library
    {
      typedef Alberta::MacroData< dim > MacroData;

      static bool checkNeighbors ( const MacroData &macroData );
      static void markLongestEdge ( MacroData &macroData );
      static void setOrientation ( [[maybe_unused]] MacroData &macroData,
                                   [[maybe_unused]] const Real orientation );

    private:
      static Real edgeLength ( const MacroData &macroData, const ElementId &e, int edge );
      static int longestEdge ( const MacroData &macroData, const ElementId &e );

      template< class Type >
      static void rotate ( Type *array, int i, int shift );

      static void rotate ( MacroData &macroData, int i, int shift );
      static void swap ( MacroData &macroData, int el, int v1, int v2 );
    };



    // Implementation of MacroData
    // ---------------------------

    template< int dim >
    inline typename MacroData< dim >::ElementId &
    MacroData< dim >::element ( int i ) const
    {
      assert( (i >= 0) && (i < data_->n_macro_elements) );
      const int offset = i * numVertices;
      return *reinterpret_cast< ElementId * >( data_->mel_vertices + offset );
    }


    template< int dim >
    inline GlobalVector &MacroData< dim >::vertex ( int i ) const
    {
      assert( (i >= 0) && (i < data_->n_total_vertices) );
      return data_->coords[ i ];
    }


    template< int dim >
    inline int &MacroData< dim >::neighbor ( int element, int i ) const
    {
      assert( (element >= 0) && (element < data_->n_macro_elements) );
      assert( (i >= 0) && (i < numVertices) );
      return data_->neigh[ element*numVertices + i ];
    }


    template< int dim >
    inline BoundaryId &MacroData< dim >::boundaryId ( int element, int i ) const
    {
      assert( (element >= 0) && (element < data_->n_macro_elements) );
      assert( (i >= 0) && (i < numVertices) );
      return data_->boundary[ element*numVertices + i ];
    }


    template< int dim >
    inline void MacroData< dim >::create ()
    {
      release();
      data_ = ALBERTA alloc_macro_data( dim, initialSize, initialSize );
      data_->boundary = memAlloc< BoundaryId >( initialSize*numVertices );
      if( dim == 3 )
        data_->el_type = memAlloc< ElementType >( initialSize );
      vertexCount_ = elementCount_ = 0;
      elementCount_ = 0;
    }


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


    template< int dim >
    inline void MacroData< dim >::markLongestEdge ()
    {
      Library< dimWorld >::markLongestEdge( *this );
    }


    template< int dim >
    inline void MacroData< dim >::setOrientation ( const Real orientation )
    {
      Library< dimWorld >::setOrientation( *this, orientation );
    }


    template< int dim >
    inline bool MacroData< dim >::checkNeighbors () const
    {
      return Library< dimWorld >::checkNeighbors( *this );
    }


    template< int dim >
    inline int MacroData< dim >::insertElement ( const ElementId &id )
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
      if( dim == 3 )
        data_->el_type[ elementCount_ ] = 0;

      return elementCount_++;
    }


    template< int dim >
    inline void MacroData< dim >
    ::insertWallTrafo ( const GlobalMatrix &matrix, const GlobalVector &shift )
    {
      int &count = data_->n_wall_trafos;
      AffineTransformation *&array = data_->wall_trafos;

      // resize wall trafo array
      array = memReAlloc< AffineTransformation >( array, count, count+1 );
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
      int &count = data_->n_wall_trafos;
      AffineTransformation *&array = data_->wall_trafos;

      // resize wall trafo array
      array = memReAlloc< AffineTransformation >( array, count, count+1 );
      assert( data_->wall_trafos != NULL );

      // copy matrix and shift
      for( int i = 0; i < dimWorld; ++i )
        copy( matrix[ i ], array[ count ].M[ i ] );
      copy( shift, array[ count ].t );
      ++count;
    }


    template< int dim >
    inline void MacroData< dim >::checkCycles ()
    {
      // ensure that the macro data has been finalized
      finalize();
      ALBERTA macro_test( data_, NULL );
    }


    template< int dim >
    inline void MacroData< dim >::read ( const std::string &filename, bool binary )
    {
      release();
      if( binary )
        data_ = ALBERTA read_macro_xdr( filename.c_str() );
      else
        data_ = ALBERTA read_macro( filename.c_str() );
    }


    template< int dim >
    inline void MacroData< dim >::resizeElements ( const int newSize )
    {
      const int oldSize = data_->n_macro_elements;
      data_->n_macro_elements = newSize;
      data_->mel_vertices = memReAlloc( data_->mel_vertices, oldSize*numVertices, newSize*numVertices );
      data_->boundary = memReAlloc( data_->boundary, oldSize*numVertices, newSize*numVertices );
      if( dim == 3 )
        data_->el_type = memReAlloc( data_->el_type, oldSize, newSize );
      assert( (newSize == 0) || (data_->mel_vertices != NULL) );
    }

  }

}

#endif // #if HAVE_ALBERTA

#endif
