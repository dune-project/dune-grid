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

    class ElementInfo
    {
      class Instance;
      class Stack;

      typedef Instance *InstancePtr;

      InstancePtr instance_;

      explicit ElementInfo ( const InstancePtr &instance );

    public:
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

      ElementInfo father () const;
      int indexInFather () const;
      ElementInfo child ( int i ) const;
      bool isLeaf () const;

      int level () const;

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

    struct ElementInfo::Instance
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

    class ElementInfo::Stack
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

    inline ElementInfo::ElementInfo ( const InstancePtr &instance )
      : instance_( instance )
    {
      addReference();
    }


    inline ElementInfo::ElementInfo ()
      : instance_( null() )
    {
      addReference();
    }


    inline ElementInfo::ElementInfo ( Mesh &mesh, MacroElement &macroElement )
    {
      instance_ = stack().allocate();
      instance_->parent() = null();
      ++(instance_->parent()->refCount);

      addReference();

#if DUNE_ALBERTA_VERSION == 0x200
      elInfo().fill_flag = FILL_ANY(&mesh);
#else
      elInfo().fill_flag = FILL_ANY;
#endif

      // Alberta fills opp_vertex only if there is a neighbor
      for( int k = 0; k < maxNeighbors; ++k )
        elInfo().opp_vertex[ k ] = -1;

      fill_macro_info( &mesh, &macroElement, &elInfo() );
    }


    inline ElementInfo::ElementInfo ( const ElementInfo &other )
      : instance_( other.instance_ )
    {
      addReference();
    }


    inline ElementInfo::~ElementInfo ()
    {
      removeReference();
    }


    inline ElementInfo &ElementInfo::operator= ( const ElementInfo &other )
    {
      other.addReference();
      removeReference();
      instance_ = other.instance_;
      return *this;
    }


    inline bool ElementInfo::operator! () const
    {
      return (instance_ == null());
    }


    inline ElementInfo ElementInfo::father () const
    {
      assert( !(*this) == false );
      return ElementInfo( instance_->parent() );
    }


    inline int ElementInfo::indexInFather () const
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


    inline ElementInfo ElementInfo::child ( int i ) const
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

      return ElementInfo( child );
    }


    inline bool ElementInfo::isLeaf () const
    {
      assert( !(*this) == false );
      return IS_LEAF_EL( el() );
    }


    inline int ElementInfo::level () const
    {
      return instance_->elInfo.level;
    }


    inline ALBERTA EL *ElementInfo::el () const
    {
      return elInfo().el;
    }


    inline ALBERTA EL_INFO &ElementInfo::elInfo () const
    {
      return (instance_->elInfo);
    }


    inline ElementInfo ElementInfo::createFake ()
    {
      InstancePtr instance = stack().allocate();
      instance->parent() = null();
      ++(instance->parent()->refCount);
      return ElementInfo( instance );
    }


    inline void ElementInfo::addReference () const
    {
      ++(instance_->refCount);
    }


    inline void ElementInfo::removeReference () const
    {
      // std::cerr << "Destructing " << (instance_ == null() ? "null" : "element")
      //           << " info (references = " << instance_->refCount << ")..."
      //           << std::flush;

      // this loop breaks when instance becomes null()
      for( InstancePtr instance = instance_; --(instance->refCount) == 0; )
      {
        const InstancePtr parent = instance->parent();
        stack().release( instance );
        instance = parent;
      }

      // std::cerr << "[done]" << std::endl;
    }


    inline ElementInfo::InstancePtr ElementInfo::null ()
    {
      return stack().null();
    }


    inline ElementInfo::Stack &ElementInfo::stack ()
    {
      static Stack s;
      return s;
    }



    // Implementation of ElementInfo::Stack
    // ------------------------------------

    inline ElementInfo::Stack::Stack ()
      : top_( 0 )
    {
      null_.elInfo.el = NULL;
      null_.refCount = 1;
      null_.parent() = 0;
    }


    inline ElementInfo::Stack::~Stack ()
    {
      while( top_ != 0 )
      {
        InstancePtr p = top_;
        top_ = p->parent();
        delete p;
      }
    }


    inline ElementInfo::InstancePtr ElementInfo::Stack::allocate ()
    {
      InstancePtr p = top_;
      if( p != 0 )
        top_ = p->parent();
      else
        p = new Instance;
      p->refCount = 0;
      return p;
    }


    inline void ElementInfo::Stack::release ( InstancePtr &p )
    {
      assert( (p != null()) && (p->refCount == 0) );
      p->parent() = top_;
      top_ = p;
    }


    inline ElementInfo::InstancePtr ElementInfo::Stack::null ()
    {
      return &null_;
    }

  }

}

#endif
