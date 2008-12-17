// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_DOFADMIN_HH
#define DUNE_ALBERTA_DOFADMIN_HH

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



    // HierarchyDofNumbering
    // ---------------------

    template< int dim >
    class HierarchyDofNumbering
    {
      typedef HierarchyDofNumbering< dim > This;

    public:
      typedef Alberta::MeshPointer< dim > MeshPointer;
      typedef Alberta::ElementInfo< dim > ElementInfo;

    private:
      template< int codim >
      struct CreateDofSpace;

      template< int codim >
      struct CacheDofSpace;

      typedef ALBERTA FE_SPACE DofSpace;
      typedef std::pair< int, int > Cache;

      MeshPointer mesh_;
      const DofSpace *dofSpace_[ dim+1 ];
      Cache cache_[ dim+1 ];

    public:
      HierarchyDofNumbering ( const MeshPointer &mesh )
        : mesh_( mesh )
      {
        ForLoop< CreateDofSpace, 0, dim >::apply( mesh_, dofSpace_ );
        ForLoop< CacheDofSpace, 0, dim >::apply( dofSpace_, cache_ );
      }

    private:
      HierarchyDofNumbering ( const This & );
      This &operator= ( const This & );

    public:
      ~HierarchyDofNumbering ()
      {
        for( int codim = 0; codim <= dim; ++codim )
          freeDofSpace( dofSpace_[ codim ] );
      }

      int operator() ( Element *element, int codim, unsigned int subEntity )
      {
        Cache &cache = cache_[ codim ];
        return element->dof[ cache.first + subEntity ][ cache.second ];
      }

      int operator() ( const ElementInfo &element, int codim, unsigned int subEntity )
      {
        return number( element.el(), codim, subEntity );
      }

    private:
      static void freeDofSpace ( const DofSpace *dofSpace );
    };



#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    inline void
    HierarchyDofNumbering< dim >::freeDofSpace ( const DofSpace *dofSpace )
    {
      ALBERTA free_fe_space( dofSpace );
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x200

#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    inline void
    HierarchyDofNumbering< dim >::freeDofSpace ( const DofSpace *dofSpace )
    {
      if( dofSpace->name != NULL )
        free( (char *)(dofSpace->name) );
      ALBERTA MEM_FREE( dofSpace, 1, DofSpace );
    }
#endif



    // HierarchyDofNumbering::CreateDofSpace
    // -------------------------------------

#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    template< int codim >
    struct HierarchyDofNumbering< dim >::CreateDofSpace
    {
      static void apply ( const MeshPointer &mesh, const DofSpace *(&dofSpace)[ dim+1 ] )
      {
        int ndof[ N_NODE_TYPES ];
        for( int i = 0; i < N_NODE_TYPES; ++i )
          ndof[ i ] = 0;
        ndof[ CodimType< dim, codim >::value ] = 1;

        std::string name = "Codimension ";
        name += (char)(codim + '0');

#if DUNE_ALBERTA_VERSION >= 0x201
        const ALBERTA FLAGS flags = ADM_PRESERVE_COARSE_DOFS;
        dofSpace[ codim ] = ALBERTA get_dof_space ( mesh, name.c_str(), ndof, flags );
#else
        dofSpace[ codim ] = ALBERTA get_fe_space ( mesh, name.c_str(), ndof, NULL, 1 );
#endif
        assert( dofSpace[ codim ] != NULL );
      }
    };
#endif // #if DUNE_ALBERTA_VERSION >= 0x200

#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    template< int codim >
    struct HierarchyDofNumbering< dim >::CreateDofSpace
    {
      static void apply ( const MeshPointer &mesh, const DofSpace *(&dofSpace)[ dim+1 ] )
      {
        int ndof[ DIM+1 ];
        for( int i = 0; i <= DIM; ++i )
          ndof[ i ] = 0;
        ndof[ CodimType< dim, codim >::value ] = 1;

        std::string name = "Codimension ";
        name += (char)(codim + '0');

        dofSpace[ codim ] = ALBERTA get_fe_space ( mesh, name.c_str(), ndof, NULL );
        assert( dofSpace[ codim ] != NULL );
      }
    };
#endif // #if DUNE_ALBERTA_VERSION < 0x200



    // HierarchyDofNumbering::CacheDofSpace
    // -------------------------------------

    template< int dim >
    template< int codim >
    struct HierarchyDofNumbering< dim >::CacheDofSpace
    {
      static void apply ( const DofSpace *(&dofSpace)[ dim+1 ], Cache (&cache)[ dim+1 ] )
      {
        const int codimtype = CodimType< dim, codim >::value;
        cache[ codim ].first = dofSpace[ codim ]->mesh->node[ codimtype ];
        cache[ codim ].second = dofSpace[ codim ]->admin->n0_dof[ codimtype ];
      }
    };

  }

}

#endif // #if HAVE_ALBERTA

#endif
