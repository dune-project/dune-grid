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
      static const int numSubEntities = NumSubEntities< dim, codim >::value;

    public:
      static const int dimension = dim;
      static const int codimension = codim;

    private:
      int node_;
#ifndef NDEBUG
      int count_;
#endif
      int index_;

    public:
      explicit DofAccess ( const DofSpace *dofSpace )
        : node_( dofSpace->admin->mesh->node[ codimtype ] ),
#ifndef NDEBUG
          count_( dofSpace->admin->n_dof[ codimtype ] ),
#endif
          index_( dofSpace->admin->n0_dof[ codimtype ] )
      {}

      int operator() ( const Element *element, int subEntity, int i ) const
      {
        assert( subEntity < numSubEntities );
#ifndef NDEBUG
        assert( i < count_ );
#endif
        return element->dof[ node_ + subEntity ][ index_ + i ];
      }

      int operator() ( const Element *element, int subEntity ) const
      {
        return (*this)( element, subEntity, 0 );
      }
    };



    // DofVectorProvider
    // -----------------

    template< class Dof >
    struct DofVectorProvider;

    template<>
    struct DofVectorProvider< int >
    {
      typedef ALBERTA DOF_INT_VEC DofVector;

      static DofVector *get ( const DofSpace *dofSpace, const std::string &name )
      {
        return ALBERTA get_dof_int_vec( name.c_str(), dofSpace );
      }

      static void free ( DofVector *dofVector )
      {
        ALBERTA free_dof_int_vec( dofVector );
      }

      static DofVector *read ( const std::string &filename, Mesh *mesh, DofSpace *dofSpace )
      {
        return ALBERTA read_dof_int_vec_xdr( filename.c_str(), mesh, dofSpace );
      }

      static bool write ( const DofVector *dofVector, const std::string &filename )
      {
        int success = ALBERTA write_dof_int_vec_xdr( dofVector, filename.c_str() );
        return (success == 0);
      }
    };

    template<>
    struct DofVectorProvider< signed char >
    {
      typedef ALBERTA DOF_SCHAR_VEC DofVector;

      static DofVector *get ( const DofSpace *dofSpace, const std::string &name )
      {
        return ALBERTA get_dof_schar_vec( name.c_str(), dofSpace );
      }

      static void free ( DofVector *dofVector )
      {
        ALBERTA free_dof_schar_vec( dofVector );
      }

      static DofVector *read ( const std::string &filename, Mesh *mesh, DofSpace *dofSpace )
      {
        return ALBERTA read_dof_schar_vec_xdr( filename.c_str(), mesh, dofSpace );
      }

      static bool write ( const DofVector *dofVector, const std::string &filename )
      {
        int success = ALBERTA write_dof_schar_vec_xdr( dofVector, filename.c_str() );
        return (success == 0);
      }
    };

    template<>
    struct DofVectorProvider< unsigned char >
    {
      typedef ALBERTA DOF_UCHAR_VEC DofVector;

      static DofVector *get ( const DofSpace *dofSpace, const std::string &name )
      {
        return ALBERTA get_dof_uchar_vec( name.c_str(), dofSpace );
      }

      static void free ( DofVector *dofVector )
      {
        ALBERTA free_dof_uchar_vec( dofVector );
      }

      static DofVector *read ( const std::string &filename, Mesh *mesh, DofSpace *dofSpace )
      {
        return ALBERTA read_dof_uchar_vec_xdr( filename.c_str(), mesh, dofSpace );
      }

      static bool write ( const DofVector *dofVector, const std::string &filename )
      {
        int success = ALBERTA write_dof_uchar_vec_xdr( dofVector, filename.c_str() );
        return (success == 0);
      }
    };

    template<>
    struct DofVectorProvider< Real >
    {
      typedef ALBERTA DOF_REAL_VEC DofVector;

      static DofVector *get ( const DofSpace *dofSpace, const std::string &name )
      {
        return ALBERTA get_dof_real_vec( name.c_str(), dofSpace );
      }

      static void free ( DofVector *dofVector )
      {
        ALBERTA free_dof_real_vec( dofVector );
      }

      static DofVector *read ( const std::string &filename, Mesh *mesh, DofSpace *dofSpace )
      {
        return ALBERTA read_dof_real_vec_xdr( filename.c_str(), mesh, dofSpace );
      }

      static bool write ( const DofVector *dofVector, const std::string &filename )
      {
        int success = ALBERTA write_dof_real_vec_xdr( dofVector, filename.c_str() );
        return (success == 0);
      }
    };

    template<>
    struct DofVectorProvider< GlobalVector >
    {
      typedef ALBERTA DOF_REAL_D_VEC DofVector;

      static DofVector *get ( const DofSpace *dofSpace, const std::string &name )
      {
        return ALBERTA get_dof_real_d_vec( name.c_str(), dofSpace );
      }

      static void free ( DofVector *dofVector )
      {
        ALBERTA free_dof_real_d_vec( dofVector );
      }

      static DofVector *read ( const std::string &filename, Mesh *mesh, DofSpace *dofSpace )
      {
        return ALBERTA read_dof_real_d_vec_xdr( filename.c_str(), mesh, dofSpace );
      }

      static bool write ( const DofVector *dofVector, const std::string &filename )
      {
        int success = ALBERTA write_dof_real_d_vec_xdr( dofVector, filename.c_str() );
        return (success == 0);
      }
    };



    // DofVector
    // ---------

    template< class Dof >
    class DofVectorPointer
    {
      typedef DofVectorPointer< Dof > This;

      typedef Alberta::DofVectorProvider< Dof > DofVectorProvider;

    public:
      typedef typename DofVectorProvider::DofVector DofVector;

    private:
      DofVector *dofVector_;

    public:
      DofVectorPointer ()
        : dofVector_( NULL )
      {}

      explicit DofVectorPointer ( DofVector *dofVector )
        : dofVector_( dofVector )
      {}

      operator DofVector * () const
      {
        return dofVector_;
      }

      operator Dof * () const
      {
        Dof *ptr = NULL;
        GET_DOF_VEC( ptr, dofVector_ );
        return ptr;
      }

      bool operator! () const
      {
        return (dofVector_ == NULL);
      }

      const DofSpace *dofSpace () const
      {
        return dofVector_->fe_space;
      }

      std::string name () const
      {
        if( dofVector_ != NULL )
          return dofVector_->name;
        else
          return std::string();
      }

      void create ( const DofSpace *dofSpace, const std::string &name = "" )
      {
        release();
        dofVector_ = DofVectorProvider::get( dofSpace, name );
      }

      template< int dim >
      void read ( const std::string &filename, const MeshPointer< dim > &meshPointer )
      {
        release();
        dofVector_ = DofVectorProvider::read( filename, meshPointer, NULL );
      }

      bool write ( const std::string &filename ) const
      {
        return DofVectorProvider::write( dofVector_, filename );
      }

      void release ()
      {
        if( dofVector_ != NULL )
        {
          DofVectorProvider::free( dofVector_ );
          dofVector_ = NULL;
        }
      }

      template< class Functor >
      void forEach ( Functor &functor )
      {
        Dof *array = (Dof *)(*this);
        FOR_ALL_DOFS( dofSpace()->admin, functor( array[ dof ] ) );
      }

      void initialize ( const Dof &value )
      {
        Dof *array = (Dof *)(*this);
        FOR_ALL_DOFS( dofSpace()->admin, array[ dof ] = value );
      }

      template< class Interpolation >
      void setupInterpolation ()
      {
        dofVector_->refine_interpol = &refineInterpolate< Interpolation >;
      }

      template< class Restriction >
      void setupRestriction ()
      {
        dofVector_->coarse_restrict = &coarsenRestrict< Restriction >;
      }

    private:
      template< class Interpolation >
      static void refineInterpolate ( DofVector *dofVector, RC_LIST_EL *list, int n )
      {
        const This dofVectorPointer( dofVector );
        Interpolation interpolation( dofVectorPointer );
        for( int i = 0; i < n; ++i )
        {
#if DUNE_ALBERTA_VERSION < 0x200
          const Element *element = list[ i ].el;
#else
          const Element *element = list[ i ].el_info.el;
#endif
          interpolation( element );
        }
      }

      template< class Restriction >
      static void coarsenRestrict ( DofVector *dofVector, RC_LIST_EL *list, int n )
      {
        const This dofVectorPointer( dofVector );
        Restriction restriction( dofVectorPointer );
        for( int i = 0; i < n; ++i )
        {
#if DUNE_ALBERTA_VERSION < 0x200
          const Element *element = list[ i ].el;
#else
          const Element *element = list[ i ].el_info.el;
#endif
          restriction( element );
        }
      }
    };



    inline void abs ( const DofVectorPointer< int > &dofVector )
    {
      int *array = (int *)dofVector;
      FOR_ALL_DOFS( dofVector.dofSpace()->admin,
                    array[ dof ] = std::abs( array[ dof ] ) );
    }


    inline int maxAbs ( const DofVectorPointer< int > &dofVector )
    {
      int *array = (int *)dofVector;
      int result = 0;
      FOR_ALL_DOFS( dofVector.dofSpace()->admin,
                    result = std::max( result, std::abs( array[ dof ] ) ) );
      return result;
    }

  }

}

#endif // #if HAVE_ALBERTA

#endif
