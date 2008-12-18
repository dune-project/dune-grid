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
      HierarchyDofNumbering ( const MeshPointer &mesh )
        : mesh_( mesh )
      {
        ForLoop< CreateDofSpace, 0, dimension >::apply( mesh_, dofSpace_ );
        ForLoop< CacheDofSpace, 0, dimension >::apply( dofSpace_, cache_ );
      }

    private:
      HierarchyDofNumbering ( const This & );
      This &operator= ( const This & );

    public:
      ~HierarchyDofNumbering ()
      {
        for( int codim = 0; codim <= dimension; ++codim )
          freeDofSpace( dofSpace_[ codim ] );
      }

      int operator() ( Element *element, int codim, unsigned int subEntity )
      {
        assert( (codim >= 0) && (codim <= dimension) );
        Cache &cache = cache_[ codim ];
        return element->dof[ cache.first + subEntity ][ cache.second ];
      }

      int operator() ( const ElementInfo &element, int codim, unsigned int subEntity )
      {
        return number( element.el(), codim, subEntity );
      }

      const DofSpace *dofSpace ( int codim ) const
      {
        assert( (codim >= 0) && (codim <= dimension) );
        return dofSpace_[ codim ];
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

      void initialize ( const Dof &value )
      {
        Dof *array = NULL;
        GET_DOF_VEC( array, dofVector_ );
        FOR_ALL_DOFS( dofSpace()->admin, array[ dof ] = value );
      }

    };

  }

}

#endif // #if HAVE_ALBERTA

#endif
