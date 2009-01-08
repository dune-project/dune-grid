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
      static const int dimension = dim;

      typedef Alberta::MeshPointer< dimension > MeshPointer;
      typedef Alberta::ElementInfo< dimension > ElementInfo;

    private:
      template< int codim >
      struct CreateDofSpace;

      template< int codim >
      struct CacheDofSpace;

      typedef std::pair< int, int > Cache;

      MeshPointer mesh_;
      const DofSpace *dofSpace_[ dimension+1 ];
      Cache cache_[ dimension+1 ];

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
        Cache &cache = cache_[ codim ];
        return element->dof[ cache.first + subEntity ][ cache.second ];
      }

      int operator() ( const ElementInfo &element, int codim, unsigned int subEntity ) const
      {
        return (*this)( element.el(), codim, subEntity );
      }

      bool operator! () const
      {
        return !mesh_;
      }

      const DofSpace *dofSpace ( int codim ) const
      {
        assert( !(*this) == false );
        assert( (codim >= 0) && (codim <= dimension) );
        return dofSpace_[ codim ];
      }

      void create ( const MeshPointer &mesh );

      void release ()
      {
        if( !(*this) )
          return;

        for( int codim = 0; codim <= dimension; ++codim )
          freeDofSpace( dofSpace_[ codim ] );
        mesh_ = MeshPointer();
      }

    private:
      static void freeDofSpace ( const DofSpace *dofSpace );
    };



    template< int dim >
    inline void
    HierarchyDofNumbering< dim >::create ( const MeshPointer &mesh )
    {
      release();

      if( !mesh )
        return;

      mesh_ = mesh;
      ForLoop< CreateDofSpace, 0, dimension >::apply( mesh_, dofSpace_ );
      ForLoop< CacheDofSpace, 0, dimension >::apply( dofSpace_, cache_ );
    }


#if DUNE_ALBERTA_VERSION >= 0x201
    template< int dim >
    inline void
    HierarchyDofNumbering< dim >::freeDofSpace ( const DofSpace *dofSpace )
    {
      ALBERTA free_fe_space( dofSpace );
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x201

#if DUNE_ALBERTA_VERSION == 0x200
    template< int dim >
    inline void
    HierarchyDofNumbering< dim >::freeDofSpace ( const DofSpace *dofSpace )
    {
      ALBERTA free_fe_space( const_cast< DofSpace * >( dofSpace ) );
    }
#endif // #if DUNE_ALBERTA_VERSION == 0x200

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

    private:
      int node_;
      int index_;
#ifndef NDEBUG
      int count_;
#endif

    public:
      DofAccess ()
        : node_( -1 )
      {}

      explicit DofAccess ( const DofSpace *dofSpace )
      {
        node_ = dofSpace->admin->mesh->node[ codimtype ];
        index_ = dofSpace->admin->n0_dof[ codimtype ];
#ifndef NDEBUG
        count_ = dofSpace->admin->n_dof[ codimtype ];
#endif
      }

      int operator() ( const Element *element, int subEntity, int i ) const
      {
#ifndef NDEBUG
        assert( node_ != -1 );
        assert( subEntity < numSubEntities );
        assert( i < count_ );
#endif
        return element->dof[ node_ + subEntity ][ index_ + i ];
      }

      int operator() ( const Element *element, int subEntity ) const
      {
        return (*this)( element, subEntity, 0 );
      }
    };

  }

}

#endif // #if HAVE_ALBERTA

#endif
