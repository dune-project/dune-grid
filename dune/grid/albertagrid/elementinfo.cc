// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
      const MacroElement &macroElement = element.macroElement();
      const MacroElement *macroNeighbor = macroElement.neighbor( face );
      if( macroNeighbor != 0 )
      {
        neighbor = ElementInfo( element.mesh(), *macroNeighbor, element.elInfo().fill_flag );
        return macroElement.opp_vertex[ face ];
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


    template<>
    template<>
    int ElementInfo< 1 >::Library< dimWorld >
    ::levelNeighbors ( const ElementInfo &element, const int face,
                       ElementInfo (&neighbor)[ maxLevelNeighbors ],
                       int (&faceInNeighbor)[ maxLevelNeighbors ] )
    {
      static const int neighborInFather[ 2 ][ numFaces ] = { {-1, 1}, {0, -1} };

      assert( !!element );

      int numNeighbors; // number of neighbors if grid is sufficiently fine
      if( element.level() > 0 )
      {
        assert( (face >= 0) && (face < numFaces) );

        const int myIndex = element.indexInFather();
        const int nbInFather = neighborInFather[ myIndex ][ face ];
        if( nbInFather >= 0 )
        {
          numNeighbors = levelNeighbors( element.father(), nbInFather, neighbor, faceInNeighbor );
          if( numNeighbors >= 0 )
          {
            if( !neighbor[ 0 ].isLeaf() )
              neighbor[ 0 ] = neighbor[ 0 ].child( 1-faceInNeighbor[ 0 ] );
            else
            {
              faceInNeighbor[ 0 ] = -1;
              numNeighbors = 0;
            }
          }
        }
        else
        {
          // the neighbor is the other child of our father
          neighbor[ 0 ] = element.father().child( 1-myIndex );
          faceInNeighbor[ 0 ] = 1-myIndex;
          numNeighbors = 1;
        }
      }
      else
      {
        // find macro level neighbors
        faceInNeighbor[ 0 ] = macroNeighbor( element, face, neighbor[ 0 ] );
        numNeighbors = (faceInNeighbor[ 0 ] >= 0);
      }

      return numNeighbors;
    }


    template<>
    template<>
    int ElementInfo< 2 >::Library< dimWorld >
    ::levelNeighbors ( const ElementInfo &element, const int face,
                       ElementInfo (&neighbor)[ maxLevelNeighbors ],
                       int (&faceInNeighbor)[ maxLevelNeighbors ] )
    {
      static const int neighborInFather[ 2 ][ numFaces ] = { {2, -1, 1}, {-1, 2, 0} };

      assert( !!element );

      int numNeighbors; // number of neighbors if grid is sufficiently fine
      if( element.level() > 0 )
      {
        assert( (face >= 0) && (face < numFaces) );

        const int myIndex = element.indexInFather();
        const int nbInFather = neighborInFather[ myIndex ][ face ];
        if( nbInFather >= 0 )
        {
          numNeighbors = levelNeighbors( element.father(), nbInFather, neighbor, faceInNeighbor );

          if( numNeighbors >= 0 )
          {
            if( nbInFather >= 2 )
            {
              // handle a refinement edge in inside
              if( faceInNeighbor[ 0 ] >= 2 )
              {
                // handle a refinement edge in outside (common refinement edge)
                assert( numNeighbors < 2 );

                int childIndex = myIndex;
                if( element.father().el()->dof[ 0 ][ 0 ] != neighbor[ 0 ].el()->dof[ 0 ][ 0 ] )
                {
                  assert( element.father().el()->dof[ 0 ][ 0 ] != neighbor[ 0 ].el()->dof[ 1 ][ 0 ] );
                  childIndex = 1-myIndex;
                }
                neighbor[ 0 ] = neighbor[ 0 ].child( childIndex );
                faceInNeighbor[ 0 ] = childIndex;
              }
              else
              {
                // handle a non-refinement edge in outside
                if( numNeighbors >= 2 )
                {
                  // drop the neighbor for the other child
                  neighbor[ 0 ] = neighbor[ myIndex ];
                  faceInNeighbor[ 0 ] = faceInNeighbor[ myIndex ];
                  numNeighbors = 1;
                }
                neighbor[ 0 ] = neighbor[ 0 ].child( 1-faceInNeighbor[ 0 ] );
                faceInNeighbor[ 0 ] = 2;
              }
            }
            else
            {
              // handle non-refinement edge in inside
              if( faceInNeighbor[ 0 ] >= 2 )
              {
                // handle refinement edge in outside
                assert( numNeighbors < 2 );
                if( !neighbor[ 0 ].isLeaf() )
                {
                  if( element.father().el()->dof[ 2 ][ 0 ] != neighbor[ 0 ].el()->dof[ faceInNeighbor[ 0 ] ][ 0 ] )
                  {
                    assert( element.father().el()->dof[ 2 ][ 0 ] == neighbor[ 0 ].el()->dof[ 1-faceInNeighbor[ 0 ] ][ 0 ] );
                    faceInNeighbor[ 0 ] = 0;
                  }
                  else
                    faceInNeighbor[ 0 ] = 1;

                  faceInNeighbor[ 1 ] = 1 - faceInNeighbor[ 0 ];
                  neighbor[ 1 ] = neighbor[ 0 ].child( faceInNeighbor[ 1 ] );
                  neighbor[ 0 ] = neighbor[ 0 ].child( faceInNeighbor[ 0 ] );
                  numNeighbors = 2;
                }
                else
                  numNeighbors = 0;
              }
              else
              {
                // handle non-refinement edge in outside
                int realNumNeighbors = 0;
                for( int i = 0; i < numNeighbors; ++i )
                {
                  assert( faceInNeighbor[ i ] < 2 );
                  if( faceInNeighbor[ i ]  < 0 )
                    continue;

                  if( !neighbor[ i ].isLeaf() )
                  {
                    neighbor[ i ] = neighbor[ i ].child( 1-faceInNeighbor[ i ] );
                    faceInNeighbor[ i ] = 2;
                    ++realNumNeighbors;
                  }
                  else
                    faceInNeighbor[ i ] = -1;
                }
                numNeighbors = (realNumNeighbors > 0 ? numNeighbors : 0);
              }
            }
          }
        }
        else
        {
          // the neighbor is the other child of our father
          neighbor[ 0 ] = element.father().child( 1-myIndex );
          faceInNeighbor[ 0 ] = myIndex;
          faceInNeighbor[ 1 ] = -1;
          numNeighbors = 1;
        }
      }
      else
      {
        // find macro level neighbors
        faceInNeighbor[ 0 ] = macroNeighbor( element, face, neighbor[ 0 ] );
        faceInNeighbor[ 1 ] = -1;
        numNeighbors = (faceInNeighbor[ 0 ] >= 0);
      }

      return numNeighbors;
    }


    template<>
    template<>
    int ElementInfo< 3 >::Library< dimWorld >
    ::levelNeighbors ( const ElementInfo &element, const int face,
                       ElementInfo (&neighbor)[ maxLevelNeighbors ],
                       int (&faceInNeighbor)[ maxLevelNeighbors ] )
    {
      assert( !!element );

      int numNeighbors;
      if( element.level() > 0 )
      {
        // we support level neighbors only on the macro level for now
        numNeighbors = 0;
      }
      else
      {
        // find macro level neighbors
        faceInNeighbor[ 0 ] = macroNeighbor( element, face, neighbor[ 0 ] );
        numNeighbors = (faceInNeighbor[ 0 ] >= 0);
      }

      return numNeighbors;
    }

  }

}

#else
#error "Library for AlbertaGrid can only be compiled if ALBERTA has been found by configure."
#endif // #if HAVE_ALBERTA
