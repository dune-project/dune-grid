// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

/** \file
 *  \author Martin Nolte
 *  \brief  provides a wrapper for ALBERTA's macro_data structure
 */

#include <dune/common/forloop.hh>

#if HAVE_ALBERTA

#include <dune/grid/albertagrid/macrodata.hh>

namespace Dune
{

  namespace Alberta
  {

    template< int dim >
    template< int dimWorld >
    bool MacroData< dim >::Library< dimWorld >
    ::checkNeighbors ( const MacroData &macroData )
    {
      assert( macroData.data_ != NULL );
      if( macroData.data_->neigh == NULL )
        return true;

#if DUNE_ALBERTA_VERSION >= 0x300
      const bool hasOppVertex = (macroData.data_->opp_vertex != NULL);
#else
      const bool hasOppVertex = false;
#endif

      const int count = macroData.elementCount();
      for( int i = 0; i < count; ++i )
      {
        for( int j = 0; j <= dimension; ++j )
        {
          const int nb = macroData.data_->neigh[ i*numVertices + j ];
          if( nb < 0 )
            continue;
          if( nb >= macroData.elementCount() )
            return false;

#if DUNE_ALBERTA_VERSION >= 0x300
          if( hasOppVertex )
          {
            const int ov = macroData.data_->opp_vertex[ i*numVertices + j ];
            if( ov > dimension )
              return false;
            if( macroData.data_->neigh[ nb*numVertices + ov ] != i )
              return false;
            if( macroData.data_->opp_vertex[ nb*numVertices + ov ] != j )
              return false;
          }
#endif // #if DUNE_ALBERTA_VERSION >= 0x300

          if( !hasOppVertex )
          {
            bool foundSelf = false;
            for( int k = 0; k <= dimension; ++k )
              foundSelf |= (macroData.data_->neigh[ nb*numVertices + k ] == i);
            if( !foundSelf )
              return false;
          }
        }
      }
      return true;
    }


    template< int dim >
    template< int dimWorld >
    void MacroData< dim >::Library< dimWorld >::markLongestEdge ( MacroData &macroData )
    {
      assert( macroData.data_ != NULL );
      if( dimension == 1 )
        return;

      const int count = macroData.elementCount();
      for( int i = 0; i < count; ++i )
      {
        const int refEdge = RefinementEdge< dimension >::value;
        const int edge = longestEdge( macroData, macroData.element( i ) );
        if( edge == refEdge )
          continue;

        // shift vertices such that the refinement edge is the longest edge
        const int shift = edge + (numVertices - refEdge);

        // rotate necessary fields
        rotate( macroData.data_->mel_vertices, i, shift );

        // correct opposite vertices
#if DUNE_ALBERTA_VERSION >= 0x300
        if( macroData.data_->opp_vertex != NULL )
        {
          assert( macroData.data_->neigh != NULL );
          const int shiftBack = numVertices - (shift % numVertices);
          for( int j = 0; j < numVertices; ++j )
          {
            const int nb = macroData.data_->neigh[ i*numVertices + j ];
            if( nb < 0 )
              continue;
            const int ov = macroData.data_->opp_vertex[ i*numVertices + j ];
            assert( macroData.data_->neigh[ nb*numVertices + ov ] == i );
            assert( macroData.data_->opp_vertex[ nb*numVertices + ov ] == j );
            macroData.data_->opp_vertex[ nb*numVertices + ov ] = (j+shiftBack) % numVertices;
          }
          rotate( macroData.data_->opp_vertex, i, shift );
        }
#endif // #if DUNE_ALBERTA_VERSION >= 0x300

        // correct neighbors and boundaries
        rotate( macroData.data_->neigh, i, shift );
        rotate( macroData.data_->boundary, i, shift );
      }
    }


    template< int dim >
    template< int dimWorld >
    void MacroData< dim >::Library< dimWorld >
    ::setOrientation ( MacroData &macroData, const Real orientation )
    {}

    template<>
    template<>
    void MacroData< dimWorld >::Library< dimWorld >
    ::setOrientation ( MacroData &macroData, const Real orientation )
    {
      assert( macroData.data_ != NULL );

      const int count = macroData.elementCount();
      for( int i = 0; i < count; ++i )
      {
        FieldMatrix< Real, dimWorld, dimWorld > jacobian;
        ElementId &id = macroData.element( i );

        const GlobalVector &x = macroData.vertex( id[ 0 ] );
        for( int j = 0; j < dimWorld; ++j )
        {
          const GlobalVector &y = macroData.vertex( id[ j+1 ] );
          for( int k = 0; k < dimWorld; ++k )
            jacobian[ j ][ k ] = y[ k ] - x[ k ];
        }
        if( determinant( jacobian ) * orientation < 0 )
          swap( macroData, i, (dimWorld == 3 ? 2 : 0), (dimWorld == 3 ? 3 : 1) );
      }
    }



    template< int dim >
    template< int dimWorld >
    inline Real MacroData< dim >::Library< dimWorld >
    ::edgeLength ( const MacroData &macroData, const ElementId &e, int edge )
    {
      const int i = MapVertices< dim, dim-1 >::apply( edge, 0 );
      assert( (macroData.vertexCount_ < 0) || (e[ i ] < macroData.vertexCount_) );
      const GlobalVector &x = macroData.vertex( e[ i ] );

      const int j = MapVertices< dim, dim-1 >::apply( edge, 1 );
      assert( (macroData.vertexCount_ < 0) || (e[ j ] < macroData.vertexCount_) );
      const GlobalVector &y = macroData.vertex( e[ j ] );

      Real sum = (y[ 0 ] - x[ 0 ]) * (y[ 0 ] - x[ 0 ]);
      for( int i = 1; i < dimWorld; ++i )
        sum += (y[ i ] - x[ i ]) * (y[ i ] - x[ i ]);
      return sqrt( sum );
    }


    template< int dim >
    template< int dimWorld >
    inline int MacroData< dim >::Library< dimWorld >
    ::longestEdge ( const MacroData &macroData, const ElementId &e )
    {
      int maxEdge = 0;
      Real maxLength = edgeLength( macroData, e, 0 );
      for( int i = 1; i < numEdges; ++i )
      {
        const Real length = edgeLength( macroData, e, i );
        if( length <= maxLength )
          continue;
        maxEdge = i;
        maxLength = length;
      }
      return maxEdge;
    }


    template< int dim >
    template< int dimWorld >
    template< class Type >
    inline void MacroData< dim >::Library< dimWorld >::rotate ( Type *array, int i, int shift )
    {
      if( array == NULL )
        return;

      const int offset = i*numVertices;
      Type old[ numVertices ];
      for( int j = 0; j < numVertices; ++j )
        old[ j ] = array[ offset + j ];
      for( int j = 0; j < numVertices; ++j )
        array[ offset + j ] = old[ (j+shift) % numVertices ];
    }


    template< int dim >
    template< int dimWorld >
    inline void MacroData< dim >::Library< dimWorld >
    ::swap ( MacroData &macroData, int el, int v1, int v2 )
    {
      std::swap( macroData.element( el )[ v1 ], macroData.element( el )[ v2 ] );
#if DUNE_ALBERTA_VERSION >= 0x300
      if( macroData.data_->opp_vertex != NULL )
      {
        assert( macroData.data_->neigh != NULL );

        const int nb1 = macroData.neighbor( el, v1 );
        if( nb1 >= 0 )
        {
          const int ov = macroData.data_->opp_vertex[ el*numVertices + v1 ];
          assert( macroData.neighbor( nb1, ov ) == el );
          assert( macroData.data_->opp_vertex[ nb1*numVertices + ov ] == v1 );
          macroData.data_->opp_vertex[ nb1*numVertices + ov ] = v2;
        }

        const int nb2 = macroData.neighbor( el, v2 );
        if( nb2 >= 0 )
        {
          const int ov = macroData.data_->opp_vertex[ el*numVertices + v2 ];
          assert( macroData.neighbor( nb2, ov ) == el );
          assert( macroData.data_->opp_vertex[ nb2*numVertices + ov ] == v2 );
          macroData.data_->opp_vertex[ nb2*numVertices + ov ] = v1;
        }

        std::swap( macroData.data_->opp_vertex[ el*numVertices + v1 ],
                   macroData.data_->opp_vertex[ el*numVertices + v2 ] );
      }
#endif // #if DUNE_ALBERTA_VERSION >= 0x300

      if( macroData.data_->neigh != NULL )
        std::swap( macroData.neighbor( el, v1 ), macroData.neighbor( el, v2 ) );
      if( macroData.data_->boundary != NULL )
        std::swap( macroData.boundaryId( el, v1 ), macroData.boundaryId( el, v2 ) );
    }



    // Instantiation
    // -------------

    template struct Dune::Alberta::MacroData< 1 >::Library< dimWorld >;
#if ALBERTA_DIM >= 2
    template struct Dune::Alberta::MacroData< 2 >::Library< dimWorld >;
#endif // #if ALBERTA_DIM >= 2
#if ALBERTA_DIM >= 3
    template struct Dune::Alberta::MacroData< 3 >::Library< dimWorld >;
#endif // #if ALBERTA_DIM >= 3

  } // end namespace Alberta

} // end namespace Dune

#endif // #if HAVE_ALBERTA
