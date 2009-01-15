// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_REFINEMENT_HH
#define DUNE_ALBERTA_REFINEMENT_HH

/** \file
 *  \author Martin Nolte
 *  \brief  provides a wrapper for ALBERTA's refinement patches
 *          and the corners for geometryInFather
 */

#include <cassert>

#include <dune/grid/albertagrid/misc.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {


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



    // GeometryInFather
    // ----------------

    template< int dim >
    struct GeometryInFather;

    template<>
    struct GeometryInFather< 1 >
    {
      static const int dim = 1;

      typedef Real LocalVector[ dim ];

      static const LocalVector &
      coordinate ( int child, int orientation, int i )
      {
        static const Real coords[ 2 ][ dim+1 ][ dim ]
          = { { {0.0}, {0.5} }, { {0.5}, {1.0} } };
        assert( (i >= 0) && (i <= dim) );
        return coords[ child ][ i ];
      }
    };

    template<>
    struct GeometryInFather< 2 >
    {
      static const int dim = 2;

      typedef Real LocalVector[ dim ];

      static const LocalVector &
      coordinate ( int child, int orientation, int i )
      {
        static const Real coords[ 2 ][ dim+1 ][ dim ]
          = { { {0.0, 1.0}, {0.0, 0.0}, {0.5, 0.0} },
              { {1.0, 0.0}, {0.0, 1.0}, {0.5, 0.0} } };
        assert( (i >= 0) && (i <= dim) );
        return coords[ child ][ i ];
      }
    };

    template<>
    struct GeometryInFather< 3 >
    {
      static const int dim = 3;

      typedef Real LocalVector[ dim ];

      static const LocalVector &
      coordinate ( int child, int orientation, int i )
      {
        static const Real coords[ 2 ][ dim+1 ][ dim ]
          = { { {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.5, 0.0, 0.0} },
              { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.5, 0.0, 0.0} } };
        static const int flip[ 2 ][ 2 ][ dim+1 ]
          = { { {0, 1, 2, 3}, {0, 1, 2, 3} }, { {0, 2, 1, 3}, {0, 1, 2, 3} } };
        assert( (i >= 0) && (i <= dim) );
        i = flip[ child ][ orientation ][ i ];
        return coords[ child ][ i ];
      }
    };

  }

}

#endif // #if HAVE_ALBERTA

#endif
