// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ELEMENTINFO_HH
#define DUNE_ALBERTA_ELEMENTINFO_HH

/** \file
 *  \author Martin Nolte
 *  \brief  provides a wrapper for ALBERTA's el_info structure
 */

#include <cassert>

#include <dune/grid/albertagrid/misc.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    // External Forward Declarations
    // -----------------------------

    template< int dim >
    class MeshPointer;



    // ElementInfo
    // -----------

    template< int dim >
    class ElementInfo
    {
      class Instance;
      class Stack;

      template< int >
      struct Library;

      typedef Instance *InstancePtr;

    public:
      static const int dimension = dim;

      static const int numVertices = NumSubEntities< dimension, dimension >::value;
      static const int numFaces = NumSubEntities< dimension, 1 >::value;

      typedef Alberta::MeshPointer< dimension > MeshPointer;
      typedef Alberta::FillFlags< dimension > FillFlags;

      static const int maxNeighbors = N_NEIGH_MAX;

      typedef ElementInfo< dimension > LevelNeighborSet[ 1 << dimension ];

    private:
      InstancePtr instance_;

      explicit ElementInfo ( const InstancePtr &instance );

    public:
      ElementInfo ();
      ElementInfo ( const MeshPointer &mesh, const MacroElement &macroElement,
                    typename FillFlags::Flags fillFlags = FillFlags::standard );
      ElementInfo ( const ElementInfo &other );

      ~ElementInfo ();

      ElementInfo &operator= ( const ElementInfo &other );

      bool operator! () const;
      bool operator== ( const ElementInfo &other ) const;
      bool operator!= ( const ElementInfo &other ) const;

      ElementInfo father () const;
      int indexInFather () const;
      ElementInfo child ( int i ) const;
      bool isLeaf () const;

      MeshPointer mesh () const;

      bool mightVanish () const;

      int level () const;
      // see ALBERTA documentation for definition of element type
      // values are 0, 1, 2
      int type () const;

      int getMark () const;
      void setMark ( int refCount ) const;

      ElementInfo leafNeighbor ( int face ) const;
      template< int codim >
      int twist ( int subEntity ) const;
      int twistInNeighbor ( int face ) const;
      bool isBoundary ( int face ) const;
      int boundaryId ( int face ) const;
      AffineTransformation *transformation ( int face ) const;

      bool hasCoordinates () const;
      const GlobalVector &coordinate ( int vertex ) const;

      template< class Functor >
      void hierarchicTraverse ( Functor &functor ) const;

      template< class Functor >
      void leafTraverse ( Functor &functor ) const;

      const Element *element () const;
      const Element *neighbor ( int face ) const;
      Element *el () const;
      ALBERTA EL_INFO &elInfo () const;

      static ElementInfo
      createFake ( const MeshPointer &mesh,
                   const Element *element, int level, int type = 0 );
      static ElementInfo createFake ( const ALBERTA EL_INFO &elInfo );

    private:
      static bool isLeaf ( Element *element );
      static bool mightVanish ( Element *element, int depth );

      //int macroNeighbor ( int face, ElementInfo &neighbor ) const;
      //int leafNeighbor ( const int face, ElementInfo &neighbor ) const;

      void addReference () const;
      void removeReference () const;

      static InstancePtr null ();
      static Stack &stack ();
    };



    // ElementInfo::Instance
    // ---------------------

    template< int dim >
    struct ElementInfo< dim >::Instance
    {
      ALBERTA EL_INFO elInfo;
      unsigned int refCount;

      InstancePtr &parent ()
      {
        return parent_;
      }

    private:
      InstancePtr parent_;
    };



    // ElementInfo::Stack
    // ------------------

    template< int dim >
    class ElementInfo< dim >::Stack
    {
      InstancePtr top_;
      Instance null_;

    public:
      Stack ();
      ~Stack ();

      InstancePtr allocate ();
      void release ( InstancePtr &p );
      InstancePtr null ();
    };



    // ElementInfo::Library
    // --------------------

    template< int dim >
    template< int >
    struct ElementInfo< dim >::Library
    {
      typedef Alberta::ElementInfo< dim > ElementInfo;

      static int
      leafNeighbor ( const ElementInfo &element, const int face, ElementInfo &neighbor );

    private:
      static int
      macroNeighbor ( const ElementInfo &element, const int face, ElementInfo &neighbor );
    };



    // Implementation of ElementInfo
    // -----------------------------

    template< int dim >
    inline ElementInfo< dim >::ElementInfo ( const InstancePtr &instance )
      : instance_( instance )
    {
      addReference();
    }


    template< int dim >
    inline ElementInfo< dim >::ElementInfo ()
      : instance_( null() )
    {
      addReference();
    }


    template< int dim >
    inline ElementInfo< dim >
    ::ElementInfo ( const MeshPointer &mesh, const MacroElement &macroElement,
                    typename FillFlags::Flags fillFlags )
    {
      instance_ = stack().allocate();
      instance_->parent() = null();
      ++(instance_->parent()->refCount);

      addReference();

      elInfo().fill_flag = fillFlags;

      // Alberta fills opp_vertex only if there is a neighbor
      for( int k = 0; k < maxNeighbors; ++k )
        elInfo().opp_vertex[ k ] = -1;

      fill_macro_info( mesh, &macroElement, &elInfo() );
    }


    template< int dim >
    inline ElementInfo< dim >::ElementInfo ( const ElementInfo &other )
      : instance_( other.instance_ )
    {
      addReference();
    }


    template< int dim >
    inline ElementInfo< dim >::~ElementInfo ()
    {
      removeReference();
    }


    template< int dim >
    inline ElementInfo< dim > &
    ElementInfo< dim >::operator= ( const ElementInfo< dim > &other )
    {
      other.addReference();
      removeReference();
      instance_ = other.instance_;
      return *this;
    }


    template< int dim >
    inline bool ElementInfo< dim >::operator! () const
    {
      return (instance_ == null());
    }


    template< int dim >
    inline bool
    ElementInfo< dim >::operator== ( const ElementInfo< dim > &other ) const
    {
      return (instance_->elInfo.el == other.instance_->elInfo.el);
    }


    template< int dim >
    inline bool
    ElementInfo< dim >::operator!= ( const ElementInfo< dim > &other ) const
    {
      return (instance_->elInfo.el != other.instance_->elInfo.el);
    }


    template< int dim >
    inline ElementInfo< dim > ElementInfo< dim >::father () const
    {
      assert( !!(*this) );
      return ElementInfo< dim >( instance_->parent() );
    }


    template< int dim >
    inline int ElementInfo< dim >::indexInFather () const
    {
      const Element *element = elInfo().el;
#if DUNE_ALBERTA_VERSION >= 0x300
      const Element *father = elInfo().parent->el;
#else
      const Element *father = elInfo().parent;
#endif
      assert( father != NULL );

      const int index = (father->child[ 0 ] == element ? 0 : 1);
      assert( father->child[ index ] == element );
      return index;
    }


    template< int dim >
    inline ElementInfo< dim > ElementInfo< dim >::child ( int i ) const
    {
      assert( !isLeaf() );

      InstancePtr child = stack().allocate();
      child->parent() = instance_;
      addReference();

      // Alberta fills opp_vertex only if there is a neighbor
      for( int k = 0; k < maxNeighbors; ++k )
        child->elInfo.opp_vertex[ k ] = -2;

#if DUNE_ALBERTA_VERSION >= 0x300
      ALBERTA fill_elinfo( i, FILL_ANY, &elInfo(), &(child->elInfo) );
#else
      ALBERTA fill_elinfo( i, &elInfo(), &(child->elInfo) );
#endif

      return ElementInfo< dim >( child );
    }


    template< int dim >
    inline bool ElementInfo< dim >::isLeaf () const
    {
      assert( !(*this) == false );
      return isLeaf( el() );
    }


    template< int dim >
    inline typename ElementInfo< dim >::MeshPointer ElementInfo< dim >::mesh () const
    {
      return MeshPointer( elInfo().mesh );
    }


    template< int dim >
    inline bool ElementInfo< dim >::mightVanish () const
    {
      return mightVanish( el(), 0 );
    }


    template< int dim >
    inline int ElementInfo< dim >::level () const
    {
      return instance_->elInfo.level;
    }


    template< int dim >
    inline int ElementInfo< dim >::type () const
    {
      return 0;
    }


    template<>
    inline int ElementInfo< 3 >::type () const
    {
      return instance_->elInfo.el_type;
    }


    template< int dim >
    inline int ElementInfo< dim >::getMark () const
    {
      return el()->mark;
    }


    template< int dim >
    inline void ElementInfo< dim >::setMark ( int refCount ) const
    {
      assert( isLeaf() );
      assert( (refCount >= -128) && (refCount < 127) );
      el()->mark = refCount;
    }


    template< int dim >
    inline ElementInfo< dim > ElementInfo< dim >::leafNeighbor ( int face ) const
    {
      assert( (face >= 0) && (face < numFaces) );
      ElementInfo neighbor;
      Library< dimWorld >::leafNeighbor( *this, face, neighbor );
      return neighbor;
    }


    template< int dim >
    template< int codim >
    inline int ElementInfo< dim >::twist ( int subEntity ) const
    {
      return Twist< dim, dim-codim >::twist( element(), subEntity );
    }


    template< int dim >
    inline int ElementInfo< dim >::twistInNeighbor ( const int face ) const
    {
      assert( neighbor( face ) != NULL );
      return Twist< dim, dim-1 >::twist( neighbor( face ), elInfo().opp_vertex[ face ] );
    }


#if DUNE_ALBERTA_VERSION >= 0x300
    template< int dim >
    inline bool ElementInfo< dim >::isBoundary ( int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < maxNeighbors) );

      const int macroFace = elInfo().macro_wall[ face ];
      if( macroFace >= 0 )
      {
        const int id = elInfo().macro_el->wall_bound[ macroFace ];
        return (id != 0);
      }
      else
        return false;
    }
#endif // DUNE_ALBERTA_VERSION >= 0x300

#if DUNE_ALBERTA_VERSION <= 0x200
    template< int dim >
    inline bool ElementInfo< dim >::isBoundary ( int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < maxNeighbors) );
      return (elInfo().neigh[ face ] == 0);
    }
#endif // DUNE_ALBERTA_VERSION <= 0x200


#if DUNE_ALBERTA_VERSION >= 0x300
    template< int dim >
    inline int ElementInfo< dim >::boundaryId ( int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < N_WALLS_MAX) );

      const int macroFace = elInfo().macro_wall[ face ];
      const int id = elInfo().macro_el->wall_bound[ macroFace ];
      // this assertion is only allowed, if FILL_BOUND is set
      // assert( id == elInfo().wall_bound[ face ] );
      return id;
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x300


#if DUNE_ALBERTA_VERSION == 0x200
    template<>
    inline int ElementInfo< 1 >::boundaryId ( int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < N_VERTICES_MAX) );
      return elInfo().vertex_bound[ face ];
    }

    template<>
    inline int ElementInfo< 2 >::boundaryId ( int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < N_EDGES_MAX) );
      return elInfo().edge_bound[ face ];
    }

    template<>
    inline int ElementInfo< 3 >::boundaryId ( int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < N_FACES_MAX) );
      return elInfo().face_bound[ face ];
    }
#endif // #if DUNE_ALBERTA_VERSION == 0x200


#if DUNE_ALBERTA_VERSION >= 0x300
    template< int dim >
    inline AffineTransformation *
    ElementInfo< dim >::transformation ( int face ) const
    {
      assert( !!(*this) );
      assert( (face >= 0) && (face < N_WALLS_MAX) );

      const int macroFace = elInfo().macro_wall[ face ];
      return (macroFace < 0 ? NULL : elInfo().macro_el->wall_trafo[ face ]);
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x300

#if DUNE_ALBERTA_VERSION <= 0x200
    template< int dim >
    inline AffineTransformation *
    ElementInfo< dim >::transformation ( int face ) const
    {
      return NULL;
    }
#endif // #if DUNE_ALBERTA_VERSION <= 0x200


    template< int dim >
    inline bool ElementInfo< dim >::hasCoordinates () const
    {
      return ((elInfo().fill_flag & FillFlags::coords) != 0);
    }

    template< int dim >
    inline const GlobalVector &ElementInfo< dim >::coordinate ( int vertex ) const
    {
      assert( hasCoordinates() );
      assert( (vertex >= 0) && (vertex < numVertices) );
      return elInfo().coord[ vertex ];
    }


    template< int dim >
    template< class Functor >
    inline void ElementInfo< dim >::hierarchicTraverse ( Functor &functor ) const
    {
      functor( *this );
      if( !isLeaf() )
      {
        child( 0 ).hierarchicTraverse( functor );
        child( 1 ).hierarchicTraverse( functor );
      }
    }


    template< int dim >
    template< class Functor >
    inline void ElementInfo< dim >::leafTraverse ( Functor &functor ) const
    {
      if( !isLeaf() )
      {
        child( 0 ).leafTraverse( functor );
        child( 1 ).leafTraverse( functor );
      }
      else
        functor( *this );
    }


    template< int dim >
    inline const Element *ElementInfo< dim >::element () const
    {
      return elInfo().el;
    }


    template< int dim >
    inline const Element *ElementInfo< dim >::neighbor ( int face ) const
    {
      assert( (face >= 0) && (face < numFaces) );
      assert( (elInfo().fill_flag & FillFlags::neighbor) != 0 );
      return elInfo().neigh[ face ];
    }


    template< int dim >
    inline Element *ElementInfo< dim >::el () const
    {
      return elInfo().el;
    }


    template< int dim >
    inline ALBERTA EL_INFO &ElementInfo< dim >::elInfo () const
    {
      return (instance_->elInfo);
    }


    template< int dim >
    inline ElementInfo< dim >
    ElementInfo< dim >::createFake ( const MeshPointer &mesh,
                                     const Element *element, int level, int type )
    {
      InstancePtr instance = stack().allocate();
      instance->parent() = null();
      ++(instance->parent()->refCount);

      instance->elInfo.mesh = mesh;
      instance->elInfo.macro_el = NULL;
      instance->elInfo.el = const_cast< Element * >( element );
      instance->elInfo.parent = NULL;
      instance->elInfo.fill_flag = FillFlags::nothing;
      instance->elInfo.level = level;
      instance->elInfo.el_type = type;

      return ElementInfo< dim >( instance );
    }


    template< int dim >
    inline ElementInfo< dim >
    ElementInfo< dim >::createFake ( const ALBERTA EL_INFO &elInfo )
    {
      InstancePtr instance = stack().allocate();
      instance->parent() = null();
      ++(instance->parent()->refCount);

      instance->elInfo = elInfo;
      return ElementInfo< dim >( instance );
    }


    template< int dim >
    inline bool ElementInfo< dim >::isLeaf ( Element *element )
    {
      return IS_LEAF_EL( element );
    }


    template< int dim >
    inline bool ElementInfo< dim >::mightVanish ( Alberta::Element *element, int depth )
    {
      if( isLeaf( element ) )
        return (element->mark < depth);
      else
        return (mightVanish( element->child[ 0 ], depth-1 ) && mightVanish( element->child[ 1 ], depth-1 ));
    }


#if 0
    template< int dim >
    inline int ElementInfo< dim >::macroNeighbor ( int face, ElementInfo &neighbor ) const
    {
      assert( (face >= 0) && (face < numFaces) );
      const MacroElement *const macroElement = elInfo().macro_el;
      const MacroElement *const macroNeighbor = macroElement->neigh[ face ];
      if( macroNeighbor != NULL )
      {
        neighbor = ElementInfo( mesh(), *macroNeighbor, elInfo().fill_flag );
        return macroElement->opp_vertex[ face ];
      }
      else
        return -1;
    }


    template<>
    inline int ElementInfo< 1 >::leafNeighbor ( const int face, ElementInfo &neighbor ) const
    {
      static const int neighborInFather[ 2 ][ numFaces ] = { {-1, 1}, {0, -1} };

      assert( !!(*this) );

      int faceInNeighbor;
      if( level() > 0 )
      {
        assert( (face >= 0) && (face < numFaces) );

        const int myIndex = indexInFather();
        const int nbInFather = neighborInFather[ myIndex ][ face ];
        if( nbInFather >= 0 )
          return father().leafNeighbor( nbInFather, neighbor );
        else
        {
          neighbor = father().child( 1-myIndex );
          faceInNeighbor = 1-myIndex;
        }
      }
      else
        faceInNeighbor = macroNeighbor( face, neighbor );

      if( faceInNeighbor >= 0 )
      {
        // refine until we are on the leaf level (faceInNeighbor < 2 is always true)
        while( !neighbor.isLeaf() )
          neighbor = neighbor.child( 1-faceInNeighbor );
        assert( neighbor.el() == elInfo().neigh[ face ] );
      }
      return faceInNeighbor;
    }


    template<>
    inline int ElementInfo< 2 >::leafNeighbor ( const int face, ElementInfo &neighbor ) const
    {
      static const int neighborInFather[ 2 ][ numFaces ] = { {2, -1, 1}, {-1, 2, 0} };

      assert( !!(*this) );

      int faceInNeighbor;
      if( level() > 0 )
      {
        assert( (face >= 0) && (face < numFaces) );

        const int myIndex = indexInFather();
        const int nbInFather = neighborInFather[ myIndex ][ face ];
        if( nbInFather >= 0 )
        {
          faceInNeighbor = father().leafNeighbor( nbInFather, neighbor );

          // handle a common face of in refinement patch
          if( (faceInNeighbor >= 0) && (nbInFather >= 2) )
          {
            assert( faceInNeighbor >= 2 );

            int childIndex = myIndex;
            if( father().el()->dof[ 0 ][ 0 ] != neighbor.el()->dof[ 0 ][ 0 ] )
            {
              assert( father().el()->dof[ 0 ][ 0 ] == neighbor.el()->dof[ 1 ][ 0 ] );
              childIndex = 1-myIndex;
            }
            neighbor = neighbor.child( childIndex );
            faceInNeighbor = childIndex;
          }
        }
        else
        {
          neighbor = father().child( 1-myIndex );
          faceInNeighbor = myIndex;
        }
      }
      else
        faceInNeighbor = macroNeighbor( face, neighbor );

      if( faceInNeighbor >= 0 )
      {
        // refine until we share a refinement face of the neighbor
        if( !neighbor.isLeaf() && (faceInNeighbor < 2) )
        {
          neighbor = neighbor.child( 1-faceInNeighbor );
          faceInNeighbor = dimension;
        }
        assert( neighbor.el() == elInfo().neigh[ face ] );
      }
      return faceInNeighbor;
    }


    template<>
    inline int ElementInfo< 3 >::leafNeighbor ( const int face, ElementInfo &neighbor ) const
    {
      // father.neigh[ neighborInFather[ child[ i ].el_type ][ i ][ j ] == child[ i ].neigh[ j ]
      static const int neighborInFather[ 3 ][ 2 ][ numFaces ]
        = { { {-1, 2, 3, 1}, {-1, 2, 3, 0} },
            { {-1, 2, 3, 1}, {-1, 3, 2, 0} },
            { {-1, 2, 3, 1}, {-1, 2, 3, 0} } };

      assert( !!(*this) );

      int faceInNeighbor;
      if( level() > 0 )
      {
        assert( (face >= 0) && (face < numFaces) );

        const int myIndex = indexInFather();
        const int nbInFather = neighborInFather[ type() ][ myIndex ][ face ];
        if( nbInFather >= 0 )
        {
          faceInNeighbor = father().leafNeighbor( nbInFather, neighbor );

          // handle a common face of in refinement patch
          if( (faceInNeighbor >= 0) && (nbInFather >= 2) )
          {
            assert( faceInNeighbor >= 2 );

            int childIndex = myIndex;
            if( father().el()->dof[ 0 ][ 0 ] != neighbor.el()->dof[ 0 ][ 0 ] )
            {
              assert( father().el()->dof[ 0 ][ 0 ] == neighbor.el()->dof[ 1 ][ 0 ] );
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
          neighbor = father().child( 1-myIndex );
          faceInNeighbor = 0;
        }
      }
      else
        faceInNeighbor = macroNeighbor( face, neighbor );

      if( faceInNeighbor >= 0 )
      {
        // refine until we share a refinement face of the neighbor
        if( !neighbor.isLeaf() && (faceInNeighbor < 2) )
        {
          neighbor = neighbor.child( 1-faceInNeighbor );
          faceInNeighbor = dimension;
        }
        assert( neighbor.el() == elInfo().neigh[ face ] );
      }
      return faceInNeighbor;
    }
#endif


    template< int dim >
    inline void ElementInfo< dim >::addReference () const
    {
      ++(instance_->refCount);
    }


    template< int dim >
    inline void ElementInfo< dim >::removeReference () const
    {
      // this loop breaks when instance becomes null()
      for( InstancePtr instance = instance_; --(instance->refCount) == 0; )
      {
        const InstancePtr parent = instance->parent();
        stack().release( instance );
        instance = parent;
      }
    }


    template< int dim >
    inline typename ElementInfo< dim >::InstancePtr
    ElementInfo< dim >::null ()
    {
      return stack().null();
    }


    template< int dim >
    inline typename ElementInfo< dim >::Stack &
    ElementInfo< dim >::stack ()
    {
      static Stack s;
      return s;
    }



    // Implementation of ElementInfo::Stack
    // ------------------------------------

    template< int dim >
    inline ElementInfo< dim >::Stack::Stack ()
      : top_( 0 )
    {
      null_.elInfo.el = NULL;
      null_.refCount = 1;
      null_.parent() = 0;
    }


    template< int dim >
    inline ElementInfo< dim >::Stack::~Stack ()
    {
      while( top_ != 0 )
      {
        InstancePtr p = top_;
        top_ = p->parent();
        delete p;
      }
    }


    template< int dim >
    inline typename ElementInfo< dim >::InstancePtr
    ElementInfo< dim >::Stack::allocate ()
    {
      InstancePtr p = top_;
      if( p != 0 )
        top_ = p->parent();
      else
        p = new Instance;
      p->refCount = 0;
      return p;
    }


    template< int dim >
    inline void ElementInfo< dim >::Stack::release ( InstancePtr &p )
    {
      assert( (p != null()) && (p->refCount == 0) );
      p->parent() = top_;
      top_ = p;
    }


    template< int dim >
    inline typename ElementInfo< dim >::InstancePtr
    ElementInfo< dim >::Stack::null ()
    {
      return &null_;
    }

  }

}

#endif // #if HAVE_ALBERTA

#endif
