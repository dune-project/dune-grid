// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_DOFVECTOR_HH
#define DUNE_ALBERTA_DOFVECTOR_HH

#include <limits>

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



    // Patch
    // -----

    class Patch
    {
      typedef ALBERTA RC_LIST_EL ElementList;

      ElementList *list_;
      int count_;

      template< int dim, int codim >
      struct ForEachInternalSubChild;

    public:
      Patch ( ElementList *list, int count )
        : list_( list ),
          count_( count )
      {
        assert( count > 0 );
      }

      Element *operator[] ( int i ) const;

      int count () const
      {
        return count_;
      }

      int elementType ( int i ) const;
      bool hasNeighbor ( int i, int neighbor ) const;
      int neighborIndex ( int i, int neighbor ) const;

      template< class Functor >
      void forEach ( Functor &functor ) const
      {
        for( int i = 0; i < count(); ++i )
          functor( (*this)[ i ] );
      }

      template< class Functor >
      void forEachInternalSubChild ( Functor &functor ) const
      {
        const int dim = Functor::dimension;
        const int codim = Functor::codimension;
        dune_static_assert( ((dim >= 1) && (dim <= 3)),
                            "Alberta supports only dimensions 1, 2, 3" );
        ForEachInternalSubChild< dim, codim >::apply( functor, *this );
      }
    };


#if DUNE_ALBERTA_VERSION < 0x200
    Element *Patch::operator[] ( int i ) const
    {
      assert( (i >= 0) && (i < count()) );
      return list_[ i ].el;
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x200

#if DUNE_ALBERTA_VERSION >= 0x200
    Element *Patch::operator[] ( int i ) const
    {
      assert( (i >= 0) && (i < count()) );
      return list_[ i ].el_info.el;
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x200


#if DUNE_ALBERTA_VERSION < 0x200
    int Patch::elementType ( int i ) const
    {
      assert( (i >= 0) && (i < count()) );
#if DIM == 3
      return list_[ i ].el_type;
#else
      return 0;
#endif
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x200

#if DUNE_ALBERTA_VERSION >= 0x200
    int Patch::elementType ( int i ) const
    {
      assert( (i >= 0) && (i < count()) );
      return list_[ i ].el_info.el_type;
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x200


#if (DUNE_ALBERTA_VERSION >= 0x200) || (DIM == 3)
    bool Patch::hasNeighbor ( int i, int neighbor ) const
    {
      return (list_[ i ].neigh[ neighbor ] != NULL);
    }

    int Patch::neighborIndex ( int i, int neighbor ) const
    {
      assert( hasNeighbor( i, neighbor ) );
      return (list_[ i ].neigh[ neighbor ]->no);
    }
#endif // #if (DUNE_ALBERTA_VERSION >= 0x200) || (DIM == 3)



    // Patch::ForEachInternalSubEntity
    // -------------------------------

    template< int dim >
    struct Patch::ForEachInternalSubChild< dim, 0 >
    {
      template< class Functor >
      static void apply ( Functor &functor, const Patch &patch )
      {
        for( int i = 0; i < patch.count(); ++i )
        {
          Element *const father = patch[ i ];
          functor( father->child[ 0 ], 0 );
          functor( father->child[ 1 ], 0 );
        }
      }
    };

    template< int dim >
    struct Patch::ForEachInternalSubChild< dim, dim >
    {
      template< class Functor >
      static void apply ( Functor &functor, const Patch &patch )
      {
        functor( patch[ 0 ]->child[ 0 ], dim );
      }
    };

    template<>
    struct Patch::ForEachInternalSubChild< 2, 1 >
    {
      template< class Functor >
      static void apply ( Functor &functor, const Patch &patch )
      {
        // see alberta/src/2d/lagrange_2_2d.c for details
        Element *const firstFather = patch[ 0 ];

        Element *const firstChild = firstFather->child[ 0 ];
        functor( firstChild, 0 );
        functor( firstChild, 1 );

        functor( firstFather->child[ 1 ], 1 );

        if( patch.count() > 1 )
        {
          Element *const father = patch[ 1 ];
          functor( father->child[ 0 ], 1 );
        }
      }
    };

    template<>
    struct Patch::ForEachInternalSubChild< 3, 1 >
    {
      template< class Functor >
      static void apply ( Functor &functor, const Patch &patch )
      {
        // see alberta/src/3d/lagrange_3_3d.c for details
        Element *const firstFather = patch[ 0 ];

        Element *const firstChild = firstFather->child[ 0 ];
        functor( firstChild, 0 );
        functor( firstChild, 1 );
        functor( firstChild, 2 );

        Element *const secondChild = firstFather->child[ 1 ];
        functor( secondChild, 1 );
        functor( secondChild, 2 );

        for( int i = 1; i < patch.count(); ++i )
        {
          Element *const father = patch[ i ];
          const int type = patch.elementType( i );

          int lr_set = 0;
          if( patch.hasNeighbor( i, 0 ) && (patch.neighborIndex( i, 0 ) < i) )
            lr_set = 1;
          if( patch.hasNeighbor( i, 1 ) && (patch.neighborIndex( i, 1 ) < i) )
            lr_set += 2;
          assert( lr_set != 0 );

          functor( father->child[ 0 ], 0 );
          switch( lr_set )
          {
          case 1 :
            functor( father->child[ 0 ], 2 );
            functor( father->child[ 1 ], (type == 0 ? 1 : 2) );
            break;

          case 2 :
            functor( father->child[ 0 ], 1 );
            functor( father->child[ 1 ], (type == 0 ? 2 : 1) );
            break;
          }
        }
      }
    };

    template<>
    struct Patch::ForEachInternalSubChild< 3, 2 >
    {
      template< class Functor >
      static void apply ( Functor &functor, const Patch &patch )
      {
        // see alberta/src/3d/lagrange_2_3d.c for details
        Element *const firstFather = patch[ 0 ];

        Element *const firstChild = firstFather->child[ 0 ];
        functor( firstChild, 2 );
        functor( firstChild, 4 );
        functor( firstChild, 5 );

        functor( firstFather->child[ 1 ], 2 );

        for( int i = 1; i < patch.count(); ++i )
        {
          Element *const father = patch[ i ];

          int lr_set = 0;
          if( patch.hasNeighbor( i, 0 ) && (patch.neighborIndex( i, 0 ) < i) )
            lr_set = 1;
          if( patch.hasNeighbor( i, 1 ) && (patch.neighborIndex( i, 1 ) < i) )
            lr_set += 2;
          assert( lr_set != 0 );

          switch( lr_set )
          {
          case 1 :
            functor( father->child[ 0 ], 4 );
            break;

          case 2 :
            functor( father->child[ 0 ], 5 );
            break;
          }
        }
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

      explicit DofVectorPointer ( const DofSpace *dofSpace,
                                  const std::string &name = "" )
        : dofVector_ ( DofVectorProvider::get( dofSpace, name ) )
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
        const Patch patch( list, n );
        Interpolation::interpolateVector( dofVectorPointer, patch );
      }

      template< class Restriction >
      static void coarsenRestrict ( DofVector *dofVector, RC_LIST_EL *list, int n )
      {
        const This dofVectorPointer( dofVector );
        const Patch patch( list, n );
        Restriction::restrictVector( dofVectorPointer, patch );
      }
    };



    inline void abs ( const DofVectorPointer< int > &dofVector )
    {
      assert( !dofVector == false );
      int *array = (int *)dofVector;
      FOR_ALL_DOFS( dofVector.dofSpace()->admin,
                    array[ dof ] = std::abs( array[ dof ] ) );
    }


    inline int max ( const DofVectorPointer< int > &dofVector )
    {
      assert( !dofVector == false );
      int *array = (int *)dofVector;
      int result = std::numeric_limits< int >::min();
      FOR_ALL_DOFS( dofVector.dofSpace()->admin,
                    result = std::max( result, array[ dof ] ) );
      return result;
    }


    inline int min ( const DofVectorPointer< int > &dofVector )
    {
      assert( !dofVector == false );
      int *array = (int *)dofVector;
      int result = std::numeric_limits< int >::max();
      FOR_ALL_DOFS( dofVector.dofSpace()->admin,
                    result = std::min( result, array[ dof ] ) );
      return result;
    }


    inline int maxAbs ( const DofVectorPointer< int > &dofVector )
    {
      assert( !dofVector == false );
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
