// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

/** \file
 *  \author Martin Nolte
 *  \brief  provides a wrapper for ALBERTA's el_info structure
 */

#if HAVE_ALBERTA

#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/albertagrid/meshpointer.hh>

namespace Dune
{

  namespace Alberta
  {

    // Implementation of ElementInfo
    // -----------------------------

    template< int dim >
    template< int dimWorld >
    inline int ElementInfo< dim >::Library< dimWorld >
    ::macroNeighbor ( const ElementInfo &element, int face, ElementInfo &neighbor )
    {
      assert( (face >= 0) && (face < numFaces) );
      const MacroElement *const macroElement = element.elInfo().macro_el;
      const MacroElement *const macroNeighbor = macroElement->neigh[ face ];
      if( macroNeighbor != NULL )
      {
        neighbor = ElementInfo( element.mesh(), *macroNeighbor, element.elInfo().fill_flag );
        return macroElement->opp_vertex[ face ];
      }
      else
        return -1;
    }


    template<>
    template<>
    int ElementInfo< 1 >::Library< dimWorld >
    ::leafNeighbor ( const ElementInfo &element, const int face, ElementInfo &neighbor )
    {
      static const int neighborInFather[ 2 ][ numFaces ] = { {-1, 1}, {0, -1} };

      assert( !!element );

      int faceInNeighbor;
      if( element.level() > 0 )
      {
        assert( (face >= 0) && (face < numFaces) );

        const int myIndex = element.indexInFather();
        const int nbInFather = neighborInFather[ myIndex ][ face ];
        if( nbInFather >= 0 )
          return leafNeighbor( element.father(), nbInFather, neighbor );
        else
        {
          neighbor = element.father().child( 1-myIndex );
          faceInNeighbor = 1-myIndex;
        }
      }
      else
        faceInNeighbor = macroNeighbor( element, face, neighbor );

      if( faceInNeighbor >= 0 )
      {
        // refine until we are on the leaf level (faceInNeighbor < 2 is always true)
        while( !neighbor.isLeaf() )
          neighbor = neighbor.child( 1-faceInNeighbor );
        assert( neighbor.el() == element.elInfo().neigh[ face ] );
      }
      return faceInNeighbor;
    }


    template<>
    template<>
    int ElementInfo< 2 >::Library< dimWorld >
    ::leafNeighbor ( const ElementInfo &element, const int face, ElementInfo &neighbor )
    {
      static const int neighborInFather[ 2 ][ numFaces ] = { {2, -1, 1}, {-1, 2, 0} };

      assert( !!element );

      int faceInNeighbor;
      if( element.level() > 0 )
      {
        assert( (face >= 0) && (face < numFaces) );

        const int myIndex = element.indexInFather();
        const int nbInFather = neighborInFather[ myIndex ][ face ];
        if( nbInFather >= 0 )
        {
          faceInNeighbor = leafNeighbor( element.father(), nbInFather, neighbor );

          // handle a common face of in refinement patch
          if( (faceInNeighbor >= 0) && (nbInFather >= 2) )
          {
            assert( faceInNeighbor >= 2 );

            int childIndex = myIndex;
            if( element.father().el()->dof[ 0 ][ 0 ] != neighbor.el()->dof[ 0 ][ 0 ] )
            {
              assert( element.father().el()->dof[ 0 ][ 0 ] == neighbor.el()->dof[ 1 ][ 0 ] );
              childIndex = 1-myIndex;
            }
            neighbor = neighbor.child( childIndex );
            faceInNeighbor = childIndex;
          }
        }
        else
        {
          neighbor = element.father().child( 1-myIndex );
          faceInNeighbor = myIndex;
        }
      }
      else
        faceInNeighbor = macroNeighbor( element, face, neighbor );

      if( faceInNeighbor >= 0 )
      {
        // refine until we share a refinement face of the neighbor
        if( !neighbor.isLeaf() && (faceInNeighbor < 2) )
        {
          neighbor = neighbor.child( 1-faceInNeighbor );
          faceInNeighbor = dimension;
        }
        assert( neighbor.el() == element.elInfo().neigh[ face ] );
      }
      return faceInNeighbor;
    }


    template<>
    template<>
    int ElementInfo< 3 >::Library< dimWorld >
    ::leafNeighbor ( const ElementInfo &element, const int face, ElementInfo &neighbor )
    {
      // father.neigh[ neighborInFather[ child[ i ].el_type ][ i ][ j ] == child[ i ].neigh[ j ]
      static const int neighborInFather[ 3 ][ 2 ][ numFaces ]
        = { { {-1, 2, 3, 1}, {-1, 2, 3, 0} },
            { {-1, 2, 3, 1}, {-1, 3, 2, 0} },
            { {-1, 2, 3, 1}, {-1, 2, 3, 0} } };

      assert( !!element );

      int faceInNeighbor;
      if( element.level() > 0 )
      {
        assert( (face >= 0) && (face < numFaces) );

        const int myIndex = element.indexInFather();
        const int nbInFather = neighborInFather[ element.type() ][ myIndex ][ face ];
        if( nbInFather >= 0 )
        {
          faceInNeighbor = leafNeighbor( element.father(), nbInFather, neighbor );

          // handle a common face of in refinement patch
          if( (faceInNeighbor >= 0) && (nbInFather >= 2) )
          {
            assert( faceInNeighbor >= 2 );

            int childIndex = myIndex;
            if( element.father().el()->dof[ 0 ][ 0 ] != neighbor.el()->dof[ 0 ][ 0 ] )
            {
              assert( element.father().el()->dof[ 0 ][ 0 ] == neighbor.el()->dof[ 1 ][ 0 ] );
              childIndex = 1-myIndex;
            }

            const int oppDof = neighbor.el()->dof[ faceInNeighbor ][ 0 ];
            neighbor = neighbor.child( childIndex );
            faceInNeighbor = (oppDof == neighbor.el()->dof[ 1 ][ 0 ] ? 1 : 2);
            assert( oppDof == neighbor.el()->dof[ faceInNeighbor ][ 0 ] );
          }
        }
        else
        {
          neighbor = element.father().child( 1-myIndex );
          faceInNeighbor = 0;
        }
      }
      else
        faceInNeighbor = macroNeighbor( element, face, neighbor );

      if( faceInNeighbor >= 0 )
      {
        // refine until we share a refinement face of the neighbor
        if( !neighbor.isLeaf() && (faceInNeighbor < 2) )
        {
          neighbor = neighbor.child( 1-faceInNeighbor );
          faceInNeighbor = dimension;
        }
        assert( neighbor.el() == element.elInfo().neigh[ face ] );
      }
      return faceInNeighbor;
    }

  }

}

#else
#error "Library for AlbertaGrid can only be compiled if ALBERTA has been found by configure."
#endif // #if HAVE_ALBERTA
