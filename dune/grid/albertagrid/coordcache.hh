// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_COORDCACHE_HH
#define DUNE_ALBERTA_COORDCACHE_HH

#include <dune/grid/albertagrid/meshpointer.hh>
#include <dune/grid/albertagrid/dofadmin.hh>
#include <dune/grid/albertagrid/dofvector.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    // CoordCache
    // ----------

    template< int dim >
    class CoordCache
    {
      typedef DofVectorPointer< GlobalVector > CoordVectorPointer;
      typedef Alberta::DofAccess< dim, dim > DofAccess;

      class LocalCaching;
      struct Interpolation;

    public:
      static const int dimension = dim;

      typedef Alberta::ElementInfo< dimension > ElementInfo;
      typedef Alberta::MeshPointer< dimension > MeshPointer;
      typedef HierarchyDofNumbering< dimension > DofNumbering;

      GlobalVector &operator() ( const Element *element, int vertex ) const
      {
        assert( !(!coords_) );
        GlobalVector *array = (GlobalVector *)coords_;
        return array[ dofAccess_( element, vertex ) ];
      }

      GlobalVector &operator() ( const ElementInfo &elementInfo, int vertex ) const
      {
        return (*this)( elementInfo.el(), vertex );
      }

      void create ( const DofNumbering &dofNumbering )
      {
        MeshPointer mesh = dofNumbering.mesh();
        const DofSpace *dofSpace = dofNumbering.dofSpace( dimension );

        coords_.create( dofSpace, "Coordinate Cache" );
        LocalCaching localCaching( coords_ );
        mesh.hierarchicTraverse( localCaching, FillFlags< dimension >::coords );
        coords_.template setupInterpolation< Interpolation >();

        dofAccess_ = DofAccess( dofSpace );
      }

      void release ()
      {
        coords_.release();
      }

    private:
      CoordVectorPointer coords_;
      DofAccess dofAccess_;
    };



    // CoordCache::LocalCaching
    // ------------------------

    template< int dim >
    class CoordCache< dim >::LocalCaching
    {
      CoordVectorPointer coords_;
      DofAccess dofAccess_;

    public:
      explicit LocalCaching ( const CoordVectorPointer &coords )
        : coords_( coords ),
          dofAccess_( coords.dofSpace() )
      {}

      void operator() ( const ElementInfo &elementInfo ) const
      {
        GlobalVector *array = (GlobalVector *)coords_;
        for( int i = 0; i < DofAccess::numSubEntities; ++i )
        {
          const GlobalVector &x = elementInfo.coordinate( i );
          GlobalVector &y = array[ dofAccess_( elementInfo.el(), i ) ];
          for( int i = 0; i < dimWorld; ++i )
            y[ i ] = x[ i ];
        }
      }
    };



    // CoordCache::Interpolation
    // -------------------------

    template< int dim >
    struct CoordCache< dim >::Interpolation
    {
      static const int dimension = dim;

      typedef Alberta::Patch< dimension > Patch;

      static void
      interpolateVector ( const CoordVectorPointer &dofVector, const Patch &patch )
      {
        DofAccess dofAccess( dofVector.dofSpace() );
        GlobalVector *array = (GlobalVector *)dofVector;

        const Element *element = patch[ 0 ];

        // new vertex is always the last one
        assert( element->child[ 0 ] != NULL );
        GlobalVector &newCoord = array[ dofAccess( element->child[ 0 ], dimension ) ];

        if( element->new_coord != NULL )
        {
          for( int j = 0; j < dimWorld; ++j )
            newCoord[ j ] = element->new_coord[ j ];
        }
        else
        {
          // new coordinate is the average of of old ones on the same edge
          // refinement edge is always between vertices 0 and 1
          const GlobalVector &coord0 = array[ dofAccess( element, 0 ) ];
          const GlobalVector &coord1 = array[ dofAccess( element, 1 ) ];
          for( int j = 0; j < dimWorld; ++j )
            newCoord[ j ] = 0.5 * (coord0[ j ] + coord1[ j ]);
        }
      }
    };

  }

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_COORDCACHE_HH
