// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_DOFADMIN_HH
#define DUNE_ALBERTA_DOFADMIN_HH

#include <utility>

#include <dune/common/hybridutilities.hh>

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/elementinfo.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    // External Forward Declarations
    // -----------------------------

    template< int dim >
    class MeshPointer;



    // DofAccess
    // ---------

    template< int dim, int codim >
    class DofAccess
    {
      static const int codimtype = CodimType< dim, codim >::value;

    public:
      static const int numSubEntities = NumSubEntities< dim, codim >::value;

      static const int dimension = dim;
      static const int codimension = codim;

      typedef Alberta::ElementInfo< dimension > ElementInfo;

      DofAccess ()
        : node_( -1 )
      {}

      explicit DofAccess ( const DofSpace *dofSpace )
      {
        assert( dofSpace );
        node_ = dofSpace->admin->mesh->node[ codimtype ];
        index_ = dofSpace->admin->n0_dof[ codimtype ];
      }

      int operator() ( const Element *element, int subEntity, int i ) const
      {
        assert( element );
        assert( node_ != -1 );
        assert( subEntity < numSubEntities );
        return element->dof[ node_ + subEntity ][ index_ + i ];
      }

      int operator() ( const Element *element, int subEntity ) const
      {
        return (*this)( element, subEntity, 0 );
      }

      int operator() ( const ElementInfo &elementInfo, int subEntity, int i ) const
      {
        return (*this)( elementInfo.el(), subEntity, i );
      }

      int operator() ( const ElementInfo &elementInfo, int subEntity ) const
      {
        return (*this)( elementInfo.el(), subEntity );
      }

    private:
      int node_;
      int index_;
    };



    // HierarchyDofNumbering
    // ---------------------

    template< int dim >
    class HierarchyDofNumbering
    {
      typedef HierarchyDofNumbering< dim > This;

    public:
      static const int dimension = dim;

      typedef Alberta::MeshPointer< dimension > MeshPointer;
      typedef Alberta::ElementInfo< dimension > ElementInfo;

    private:
      static const int nNodeTypes = N_NODE_TYPES;

      template< int codim >
      struct CreateDofSpace;

      template< int codim >
      struct CacheDofSpace;

      typedef std::pair< int, int > Cache;

    public:
      HierarchyDofNumbering ()
      {}

    private:
      HierarchyDofNumbering ( const This & );
      This &operator= ( const This & );

    public:
      ~HierarchyDofNumbering ()
      {
        release();
      }

      int operator() ( const Element *element, int codim, unsigned int subEntity ) const
      {
        assert( !(*this) == false );
        assert( (codim >= 0) && (codim <= dimension) );
        const Cache &cache = cache_[ codim ];
        return element->dof[ cache.first + subEntity ][ cache.second ];
      }

      int operator() ( const ElementInfo &element, int codim, unsigned int subEntity ) const
      {
        return (*this)( element.el(), codim, subEntity );
      }

      explicit operator bool () const
      {
        return (bool)mesh_;
      }

      const DofSpace *dofSpace ( int codim ) const
      {
        assert( *this );
        assert( (codim >= 0) && (codim <= dimension) );
        return dofSpace_[ codim ];
      }

      const DofSpace *emptyDofSpace () const
      {
        assert( *this );
        return emptySpace_;
      }

      const MeshPointer &mesh () const
      {
        return mesh_;
      }

      int size ( int codim ) const
      {
        return dofSpace( codim )->admin->size;
      }

      void create ( const MeshPointer &mesh );

      void release ()
      {
        if( *this )
        {
          for( int codim = 0; codim <= dimension; ++codim )
            freeDofSpace( dofSpace_[ codim ] );
          freeDofSpace( emptySpace_ );
          mesh_ = MeshPointer();
        }
      }

    private:
      static const DofSpace *createEmptyDofSpace ( const MeshPointer &mesh );
      static const DofSpace *createDofSpace ( const MeshPointer &mesh,
                                              const std::string &name,
                                              const int (&ndof)[ nNodeTypes ],
                                              const bool periodic = false );
      static void freeDofSpace ( const DofSpace *dofSpace );

      MeshPointer mesh_;
      const DofSpace *emptySpace_;
      const DofSpace *dofSpace_[ dimension+1 ];
      Cache cache_[ dimension+1 ];
    };



    template< int dim >
    inline void
    HierarchyDofNumbering< dim >::create ( const MeshPointer &mesh )
    {
      release();

      if( !mesh )
        return;

      mesh_ = mesh;

      Hybrid::forEach( std::make_index_sequence< dimension+1 >{}, [ & ]( auto i ){ CreateDofSpace< i >::apply( mesh_, dofSpace_ ); } );
      Hybrid::forEach( std::make_index_sequence< dimension+1 >{}, [ & ]( auto i ){ CacheDofSpace< i >::apply( dofSpace_, cache_ ); } );

      emptySpace_ = createEmptyDofSpace( mesh_ );
      for( int i = 0; i < nNodeTypes; ++i )
        assert( emptySpace_->admin->n_dof[ i ] == 0 );
    }



    template< int dim >
    inline const DofSpace *
    HierarchyDofNumbering< dim >::createEmptyDofSpace ( const MeshPointer &mesh )
    {
      int ndof[ nNodeTypes ];
      for( int i = 0; i < nNodeTypes; ++i )
        ndof[ i ] = 0;
      std::string name = "Empty";
      return createDofSpace( mesh, name, ndof );
    }


    template< int dim >
    inline const DofSpace *
    HierarchyDofNumbering< dim >::createDofSpace ( const MeshPointer &mesh,
                                                   const std::string &name,
                                                   const int (&ndof)[ nNodeTypes ],
                                                   const bool periodic )
    {
      const ALBERTA FLAGS flags
        = ADM_PRESERVE_COARSE_DOFS | (periodic ? ADM_PERIODIC : 0);
      return ALBERTA get_dof_space ( mesh, name.c_str(), ndof, flags );
    }


    template< int dim >
    inline void
    HierarchyDofNumbering< dim >::freeDofSpace ( const DofSpace *dofSpace )
    {
      ALBERTA free_fe_space( dofSpace );
    }



    // HierarchyDofNumbering::CreateDofSpace
    // -------------------------------------

    template< int dim >
    template< int codim >
    struct HierarchyDofNumbering< dim >::CreateDofSpace
    {
      static void apply ( const MeshPointer &mesh, const DofSpace *(&dofSpace)[ dim+1 ] )
      {
        int ndof[ nNodeTypes ];
        for( int i = 0; i < nNodeTypes; ++i )
          ndof[ i ] = 0;
        ndof[ CodimType< dim, codim >::value ] = 1;

        std::string name = "Codimension ";
        name += (char)(codim + '0');

        dofSpace[ codim ] = createDofSpace( mesh, name, ndof );
        assert( dofSpace[ codim ] );
      }
    };



    // HierarchyDofNumbering::CacheDofSpace
    // ------------------------------------

    template< int dim >
    template< int codim >
    struct HierarchyDofNumbering< dim >::CacheDofSpace
    {
      static void apply ( const DofSpace *(&dofSpace)[ dim+1 ], Cache (&cache)[ dim+1 ] )
      {
        assert( dofSpace[ codim ] );
        const int codimtype = CodimType< dim, codim >::value;
        cache[ codim ].first = dofSpace[ codim ]->mesh->node[ codimtype ];
        cache[ codim ].second = dofSpace[ codim ]->admin->n0_dof[ codimtype ];
      }
    };

  } // namespace Alberta

} // namespace Dune

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_DOFADMIN_HH
