// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ELEMENTINFO_HH
#define DUNE_ALBERTA_ELEMENTINFO_HH

namespace Dune
{

  namespace Alberta
  {

    typedef ALBERTA MESH Mesh;
    typedef ALBERTA MACRO_EL MacroElement;



    // NumSubEntities
    // --------------

    template< int dim, int codim >
    struct NumSubEntities;

    template< int dim >
    struct NumSubEntities< dim, 0 >
    {
      static const int value = 1;
    };

    template< int dim >
    struct NumSubEntities< dim, dim >
    {
      static const int value = dim+1;
    };

    template<>
    struct NumSubEntities< 2, 1 >
    {
      static const int value = 3;
    };

    template<>
    struct NumSubEntities< 3, 1 >
    {
      static const int value = 4;
    };

    template<>
    struct NumSubEntities< 3, 2 >
    {
      static const int value = 6;
    };



    // ElementInfo
    // -----------

    template< int dim >
    class ElementInfo
    {
      class Instance;
      class Stack;

      typedef Instance *InstancePtr;

      InstancePtr instance_;

      explicit ElementInfo ( const InstancePtr &instance );

    public:
      static const int dimension = dim;

#if DUNE_ALBERTA_VERSION >= 0x200
      static const int maxNeighbors = N_NEIGH_MAX;
#else
      static const int maxNeighbors = N_NEIGH;
#endif

      ElementInfo ();
      ElementInfo ( Mesh &mesh, MacroElement &macroElement );
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

      int level () const;
      // see ALBERTA documentation for definition of element type
      // values are 0, 1, 2
      int type () const;

      bool isBoundary ( int face ) const;
      int boundaryId ( int face ) const;

      ALBERTA EL *el () const;
      ALBERTA EL_INFO &elInfo () const;

      static ElementInfo createFake ();

    private:
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
        //return reinterpret_cast< InstancePtr & >( elInfo.parent );
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
    inline ElementInfo< dim >::ElementInfo ( Mesh &mesh, MacroElement &macroElement )
    {
      instance_ = stack().allocate();
      instance_->parent() = null();
      ++(instance_->parent()->refCount);

      addReference();

#if DUNE_ALBERTA_VERSION >= 0x201
      elInfo().fill_flag = FILL_COORDS | FILL_NEIGH | FILL_OPP_COORDS
                           | FILL_ORIENTATION | FILL_MACRO_WALLS
                           | FILL_NON_PERIODIC;
#elif DUNE_ALBERTA_VERSION == 0x200
      elInfo().fill_flag = FILL_ANY(&mesh);
#else
      elInfo().fill_flag = FILL_ANY;
#endif

      // Alberta fills opp_vertex only if there is a neighbor
      for( int k = 0; k < maxNeighbors; ++k )
        elInfo().opp_vertex[ k ] = -1;

      fill_macro_info( &mesh, &macroElement, &elInfo() );
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
      assert( !(*this) == false );
      return ElementInfo< dim >( instance_->parent() );
    }


    template< int dim >
    inline int ElementInfo< dim >::indexInFather () const
    {
      const ALBERTA EL *element = elInfo().el;
#if DUNE_ALBERTA_VERSION >= 0x201
      const ALBERTA EL *father = elInfo().parent->el;
#else
      const ALBERTA EL *father = elInfo().parent;
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

#if DUNE_ALBERTA_VERSION >= 0x201
      ALBERTA fill_elinfo( i, ALBERTA FILL_ANY, &elInfo(), &(child->elInfo) );
#else
      ALBERTA fill_elinfo( i, &elInfo(), &(child->elInfo) );
#endif

      return ElementInfo< dim >( child );
    }


    template< int dim >
    inline bool ElementInfo< dim >::isLeaf () const
    {
      assert( !(*this) == false );
      return IS_LEAF_EL( el() );
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


#if (DUNE_ABLERTA_VERSION >= 0x200) || (DIM == 3)
    template<>
    inline int ElementInfo< 3 >::type () const
    {
      return instance_->elInfo.type;
    }
#endif


    template< int dim >
    inline bool ElementInfo< dim >::isBoundary ( int face ) const
    {
      assert( !(*this) == false );
      assert( (face >= 0) && (face < maxNeighbors) );
      return (elInfo().neigh[ face ] == 0);
    }


#if DUNE_ALBERTA_VERSION >= 0x201
    template< int dim >
    inline int ElementInfo< dim >::boundaryId ( int face ) const
    {
      assert( !(*this) == false );
      assert( (face >= 0) && (face < N_WALLS_MAX) );

      const int macroFace = elInfo().macro_wall[ face ];
      const int id = elInfo().macro_el->wall_bound[ macroFace ];
      // this assertion is only allowed, if FILL_BOUND is set
      // assert( id == elInfo().wall_bound[ face ] );
      return id;
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x201


#if DUNE_ALBERTA_VERSION == 0x200
    template<>
    inline int ElementInfo< 1 >::boundaryId ( int face ) const
    {
      assert( !(*this) == false );
      assert( (face >= 0) && (face < N_VERTICES_MAX) );
      return elInfo().vertex_bound[ face ];
    }

    template<>
    inline int ElementInfo< 2 >::boundaryId ( int face ) const
    {
      assert( !(*this) == false );
      assert( (face >= 0) && (face < N_EDGES_MAX) );
      return elInfo().edge_bound[ face ];
    }

    template<>
    inline int ElementInfo< 3 >::boundaryId ( int face ) const
    {
      assert( !(*this) == false );
      assert( (face >= 0) && (face < N_FACES_MAX) );
      return elInfo().face_bound[ face ];
    }
#endif // #if DUNE_ALBERTA_VERSION == 0x200


#if DUNE_ALBERTA_VERSION < 0x200
#if DIM == 1
    template<>
    inline int ElementInfo< 1 >::boundaryId ( int face ) const
    {
      assert( !(*this) == false );
      assert( (face >= 0) && (face < N_VERTICES) );
      return elInfo().bound[ face ];
    }
#endif // #if DIM == 1

#if DIM == 2
    template<>
    inline int ElementInfo< 2 >::boundaryId ( int face ) const
    {
      assert( !(*this) == false );
      assert( (face >= 0) && (face < N_EDGES) );
      return elInfo().boundary[ face ]->bound;
    }
#endif // #if DIM == 2

#if DIM == 3
    template<>
    inline int ElementInfo< 3 >::boundaryId ( int face ) const
    {
      assert( !(*this) == false );
      assert( (face >= 0) && (face < N_FACES) );
      return elInfo().boundary[ face ]->bound;
    }
#endif // #if DIM == 3
#endif // #if DUNE_ALBERTA_VERSION < 0x200


    template< int dim >
    inline ALBERTA EL *ElementInfo< dim >::el () const
    {
      return elInfo().el;
    }


    template< int dim >
    inline ALBERTA EL_INFO &ElementInfo< dim >::elInfo () const
    {
      return (instance_->elInfo);
    }


    template< int dim >
    inline ElementInfo< dim > ElementInfo< dim >::createFake ()
    {
      InstancePtr instance = stack().allocate();
      instance->parent() = null();
      ++(instance->parent()->refCount);
      return ElementInfo< dim >( instance );
    }


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

#endif
