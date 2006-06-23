// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDITERATOR_HH
#define DUNE_ALU3DGRIDITERATOR_HH

// System includes

// Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/intersectioniteratorwrapper.hh>

// Local includes
#include "alu3dinclude.hh"
#include "topology.hh"
#include "faceutility.hh"
#include "myautoptr.hh"

namespace ALUGridSpace {
  //*************************************************************
  //  definition of original LeafIterators of ALUGrid
  //
  // default is element (codim = 0)
  template <int codim>
  struct BSMacroIterator
  {
    typedef AccessIterator<GitterType::helement_STI>::Handle IteratorType;
  };

  //******************************************************************
  //  LevelIterators
  //******************************************************************
  template <int codim> struct ALUHElementType
  {
    typedef GitterType :: helement_STI ElementType;
  };

  template <> struct ALUHElementType<0> {
    typedef GitterType :: helement_STI ElementType;
  };
  template <> struct ALUHElementType<1> {
    typedef GitterType :: hface_STI ElementType;
  };
  template <> struct ALUHElementType<2> {
    typedef GitterType :: hedge_STI ElementType;
  };
  template <> struct ALUHElementType<3> {
    typedef GitterType :: vertex_STI ElementType;
  };

  template <int codim> struct BSIMPLElementType
  {
    typedef GitterImplType::Objects::tetra_IMPL ElementType; // impl Element
  };

  template <> struct BSIMPLElementType<0> {
    typedef GitterImplType::Objects::tetra_IMPL ElementType; // impl Element
  };
  template <> struct BSIMPLElementType<1> {
    typedef GitterImplType::Objects::hface3_IMPL ElementType; // impl Element
  };
  template <> struct BSIMPLElementType<2> {
    typedef GitterImplType::Objects::hedge1_IMPL ElementType; // impl Element
  };

  template <> struct BSIMPLElementType<3> {
    typedef GitterType::Geometric::VertexGeo ElementType;
  };

  typedef Dune :: ALU3dGridVertexList VertexListType;

  //*********************************************************
  //  LevelIterator Wrapper
  //*********************************************************
  template <class val_t>
  class IteratorWrapperInterface : public IteratorSTI< val_t >
  {
  public:
    virtual ~IteratorWrapperInterface () {}

    virtual int size  () = 0;
    virtual void next () = 0;
    virtual void first() = 0;
    virtual int done  () const = 0;
    virtual val_t & item () const = 0;
    virtual IteratorSTI< val_t > * clone () const { assert(false); return 0; }

  };

  typedef Dune::PartitionIteratorType PartitionIteratorType;

  // defines the pair of element and boundary
  template <int codim>
  struct IteratorElType
  {
    typedef typename ALUHElementType<codim>::ElementType ElType;
    typedef pair < ElType * , HBndSegType * > val_t;
  };

  template <int codim, PartitionIteratorType pitype>
  class ALU3dGridLevelIteratorWrapper;

  // the element level iterator
  template <PartitionIteratorType pitype>
  class ALU3dGridLevelIteratorWrapper<0,pitype>
    : public IteratorWrapperInterface< typename IteratorElType<0>::val_t >
  {
    typedef ALUHElementType<0>::ElementType ElType;
    typedef ALU3DSPACE LevelIterator < ElType > IteratorType;

    // the iterator
    IteratorType it_;
  public:
    typedef typename IteratorElType<0>::val_t val_t;
    mutable val_t elem_;

    // constructor creating iterator
    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid,
                                   const VertexListType & , int level , const int nlinks )
      : it_(grid.myGrid() , level )
        , elem_(0,0)
    {}

    // copy constructor
    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : it_( org.it_ ), elem_(org.elem_)
    {}

    int size  ()    { return it_->size(); }
    void next ()    { it_->next();  }
    void first()    { it_->first(); }
    int done () const { return it_->done(); }
    val_t & item () const
    {
      assert( ! done () );
      elem_.first  = & it_->item();
      return elem_;
    }
  };

  // the face level iterator
  template <PartitionIteratorType pitype>
  class ALU3dGridLevelIteratorWrapper<1,pitype>
    : public IteratorWrapperInterface< typename IteratorElType<1>::val_t >
  {
    typedef ALUHElementType<1>::ElementType ElType;
    typedef ALU3DSPACE LevelIterator < ElType > IteratorType;

    // the iterator
    IteratorType it_;
  public:
    typedef typename IteratorElType<1>::val_t val_t;
    mutable val_t elem_;

    // constructor creating iterator
    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid,
                                   const VertexListType & , int level , const int nlinks )
      : it_(grid.myGrid() , level )
        , elem_(0,0)
    {}

    // copy constructor
    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : it_( org.it_ ), elem_(org.elem_)
    {}

    int size  ()    { return it_->size(); }
    void next ()    { it_->next();  }
    void first()    { it_->first(); }
    int done () const { return it_->done(); }
    val_t & item () const
    {
      assert( ! done () );
      elem_.first  = & it_->item();
      return elem_;
    }
  };
  // the vertex level iterator, little bit different to the others
  // this implementation uses the vertex leaf iterator and runs over all
  // vertices with level <= the given iteration level
  template <PartitionIteratorType pitype>
  class ALU3dGridLevelIteratorWrapper<3,pitype>
    : public IteratorWrapperInterface< typename IteratorElType<3>::val_t >
  {
    typedef VertexListType :: IteratorType IteratorType;

    mutable VertexListType & vxList_;

    mutable int count_;
    const int size_;

  public:
    typedef IteratorElType<3>::val_t val_t;
    mutable val_t elem_;

    // constructor creating iterator
    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid,
                                   VertexListType & vxList , int level , const int nlinks )
      : vxList_ (vxList)
        , count_(0)
        , size_(vxList.size())
        , elem_(0,0)
    {
      assert( vxList_.up2Date() );
    }

    // copy constructor
    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : vxList_(org.vxList_) , count_(org.count_) , size_(org.size_) , elem_(org.elem_)
    {}

    // returns size of leaf iterator, wrong here, return leaf size
    int size  ()  { return size_; }

    //! if level of item is larger then walk level, go next
    void next ()
    {
      if( done () ) return ;
      ++count_;
      return ;
    }

    void first()  { count_ = 0; }
    int done () const { return (count_ >= size_) ? 1 : 0; }
    val_t & item () const
    {
      assert( ! done () );
      elem_.first = vxList_.getItemList()[count_];
      assert( elem_.first );
      return elem_;
    }
  };

  template <int codim, PartitionIteratorType pitype> class ALU3dGridLeafIteratorWrapper;
  typedef pair <ALUHElementType<0>::ElementType * , HBndSegType * > LeafValType;
  typedef IteratorWrapperInterface<LeafValType> IteratorWrapperInterfaceType;

  //**********************************************************
  //  LeafIterator Wrapper
  //**********************************************************
  template <PartitionIteratorType pitype>
  class ALU3dGridLeafIteratorWrapper<0,pitype>
    : public IteratorWrapperInterface< typename IteratorElType<0>::val_t >
  {
    // type is helement_STI
    typedef IteratorElType<0>::ElType ElType;
    typedef LeafIterator < ElType > IteratorType;

    // the ALU3dGrid Iterator
    IteratorType it_;

  public:
    typedef typename IteratorElType<0>::val_t val_t;
  private:
    mutable val_t elem_;
  public:
    // constructor creating Iterator
    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level, const int links )
      : it_(grid.myGrid()), elem_(0,0) {}

    // constructor copying iterator
    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper  & org )
      : it_( org.it_ ), elem_(org.elem_)
    {}

    int size  ()    { return it_->size(); }
    void next ()    { it_->next(); }
    void first()    { it_->first(); }
    int done () const { return it_->done(); }
    val_t & item () const
    {
      assert( ! done () );
      elem_.first  = & it_->item();
      return elem_;
    }
  };

  template <class ElType, PartitionIteratorType pitype>
  struct LeafStopRule
  {
    typedef is_leaf_entity< ElType > StopRule_t;
  };

  // only in parallel we need only the interior items, in serial all items
  // are interior, to make the check fasterm this is only in parallel
  // implemented
#if ALU3DGRID_PARALLEL
  template <class ElType>
  struct LeafStopRule<ElType, Dune :: Interior_Partition>
  {
    typedef is_interior_leaf_entity< ElType > StopRule_t;
  };
#endif

  template <PartitionIteratorType pitype>
  class ALU3dGridLeafIteratorWrapper<1,pitype>
    : public IteratorWrapperInterface < typename IteratorElType<1>::val_t >
  {
    // type is hface_STI
    typedef IteratorElType<1>::ElType ElType;
    typedef typename LeafStopRule< ElType, pitype > :: StopRule_t StopRule_t;
    typedef GridIterator < ElType , StopRule_t > IteratorType;

    // the face iterator
    IteratorType it_;

  public:
    typedef IteratorElType<1>::val_t val_t;
  private:
    mutable val_t elem_;
  public:
    // constructor creating Iterator
    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level, const int links )
      : it_(grid.myGrid(), StopRule_t() ) , elem_(0,0) {}

    // constructor copying iterator
    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper  & org )
      : it_( org.it_ ), elem_(org.elem_) {}

    ~ALU3dGridLeafIteratorWrapper ()
    {}

    int size  ()    { return it_->size(); }
    void next ()    { it_->next(); }
    void first()    { it_->first(); }
    int done () const { return it_->done(); }
    val_t & item () const
    {
      assert( ! done () );
      elem_.first  = & it_->item();
      return elem_;
    }
  };

  template <PartitionIteratorType pitype>
  class ALU3dGridLeafIteratorWrapper<2,pitype>
    : public IteratorWrapperInterface < typename IteratorElType<2>::val_t >
  {
    // type of hedge_STI
    typedef IteratorElType<2>::ElType ElType;
    typedef typename LeafStopRule< ElType, pitype > :: StopRule_t StopRule_t;
    typedef GridIterator < ElType , StopRule_t > IteratorType;

  public:
    typedef IteratorElType<2>::val_t val_t;
  private:
    // the edge iterator
    IteratorType it_;

    mutable val_t elem_;
  public:
    // constructor creating Iterator
    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level, const int links )
      : it_(grid.myGrid(), StopRule_t() ), elem_(0,0) {}

    // constructor copying iterator
    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper  & org )
      : it_( org.it_ ), elem_(org.elem_) {}

    int size  ()      { return it_->size(); }
    void next ()      { it_->next(); }
    void first()      { it_->first(); }
    int done () const { return it_->done(); }
    val_t & item () const
    {
      assert( ! done () );
      elem_.first  = & it_->item();
      return elem_;
    }
  };

  template <PartitionIteratorType pitype>
  class ALU3dGridLeafIteratorWrapper<3,pitype>
    : public IteratorWrapperInterface < typename IteratorElType<3>::val_t >
  {
    typedef IteratorElType<3>::ElType ElType;
    typedef typename LeafStopRule< ElType, pitype > :: StopRule_t StopRule_t;
    // ElType is vertex_STI
    typedef LeafIterator < ElType > IteratorType;

    // the vertex iterator
    IteratorType it_;

  public:
    typedef IteratorElType<3>::val_t val_t;
  private:
    mutable val_t elem_;
    const StopRule_t rule_;
  public:
    // constructor creating Iterator
    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level, const int links )
      : it_(grid.myGrid()), elem_(0,0) , rule_ () {}

    // constructor copying iterator
    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper  & org )
      : it_( org.it_ ), elem_(org.elem_) , rule_() {}

    int size  ()    { return it_->size(); }

    void next ()
    {
      it_->next();
      if (!it_->done())
      {
        // take standard walk rule to cehck vertices again, see walk.h
        if(! rule_(it_->item()) ) next();
      }
    }

    void first()     { it_->first(); }
    int done () const { return it_->done(); }
    val_t & item () const
    {
      assert( ! done () );
      elem_.first  = & it_->item();
      return elem_;
    }
  };

#if ALU3DGRID_PARALLEL
  template <int codim>
  class LeafLevelIteratorTTProxy
  {
    // type is hface_STI or hedge_STI
    typedef typename ALUHElementType<codim>::ElementType ElType;

    typedef IteratorSTI < ElType > IteratorType;
    IteratorType * inner_;
    IteratorType * outer_;
  public:

    // constructor creating leafBorderIteratorTT
    LeafLevelIteratorTTProxy( GitterImplType & gitter , int link )
    {
      pair < IteratorSTI< ElType > * , IteratorSTI< ElType > * >
      p = gitter.leafBorderIteratorTT( (ElType *) 0 , link );

      inner_ = p.first;
      outer_ = p.second;
    }

    // constructor creating levelBorderIteratorTT
    LeafLevelIteratorTTProxy( GitterImplType & gitter , int link , int level )
    {
      pair < IteratorSTI< ElType > * , IteratorSTI< ElType > * >
      p = gitter.levelBorderIteratorTT( (ElType *) 0 , link , level );

      inner_ = p.first;
      outer_ = p.second;
    }

    LeafLevelIteratorTTProxy( const LeafLevelIteratorTTProxy & org )
      : inner_(org.inner_->clone())
        , outer_(org.outer_->clone())
    {}

    ~LeafLevelIteratorTTProxy()
    {
      delete inner_;
      delete outer_;
    }

    IteratorType & inner () { assert(inner_); return *inner_; }
    IteratorType & outer () { assert(outer_); return *outer_; }
  };

  //****************************
  //
  //  --GhostIterator
  //
  //****************************
  class ALU3dGridGhostIterator
    : public IteratorWrapperInterface< LeafValType >
  {
  protected:
    GitterImplType & gitter_;

    // this tpye is hface_STI
    typedef ALUHElementType<1>::ElementType ElType;

    typedef LeafLevelIteratorTTProxy<1> IteratorType;

    IteratorType * iterTT_;

    typedef IteratorSTI < ElType > InnerIteratorType;
    InnerIteratorType * it_;

    // number of links
    const int nl_;

    // current link
    int link_;

    bool usingInner_;
  public:
    typedef LeafValType val_t;
  private:
    // the pair of elementand boundary face
    mutable val_t elem_;
  public:
    typedef ElementPllXIF_t ItemType;

    template <class GridImp>
    ALU3dGridGhostIterator (const GridImp & grid, int level , const int nlinks )
      : gitter_(grid.myGrid())
        , iterTT_(0) , it_(0)
        , nl_(nlinks)
        , link_(nlinks) // makes default status == done
        , elem_(0,0)
    {}

    ALU3dGridGhostIterator (const ALU3dGridGhostIterator & org)
      : gitter_(org.gitter_)
        , iterTT_(0) , it_(0)
        , nl_(org.nl_)
        , link_(org.link_)
        , usingInner_(false)
        , elem_(org.elem_)
    {
      if( org.iterTT_ )
      {
        iterTT_ = new IteratorType ( *org.iterTT_ );
        usingInner_ = org.usingInner_;
        if( org.it_ )
        {
          assert( ! org.it_->done() );
          it_ = (usingInner_) ? &( iterTT_->inner() ) : &( iterTT_->outer() );
        }
      }
    }

    ~ALU3dGridGhostIterator ()
    {
      removeIterators();
    }

  protected:
    virtual IteratorType * newIterator()
    {
      return new IteratorType ( gitter_, link_ );
    }

    void removeIterators()
    {
      if(iterTT_) delete iterTT_;
      iterTT_ = 0;
      it_ = 0;
      usingInner_ = false;
    }

    void createIterator()
    {
      if (usingInner_) checkInnerOuter();

      if (!usingInner_)
      {
        ++link_;

        removeIterators();
        if(link_ < nl_)
        {
          iterTT_ = newIterator();
          assert(iterTT_);
          checkInnerOuter();
          if (!it_) createIterator();
        }
      }
    }

    void checkInnerOuter()
    {
      it_ = 0;
      if (!usingInner_)
      {
        assert(iterTT_);
        it_ = &( iterTT_->inner() );
        InnerIteratorType & it = iterTT_->inner();
        it.first();
        if(!it.done())
        {
          usingInner_ = true;
          pair < ElementPllXIF_t *, int > p = it.item ().accessPllX ().accessOuterPllX () ;
          pair < HElementType * , HBndSegType * > elems(0,0);
          p.first->getAttachedElement(elems);

          assert( elems.first || elems.second );

          if(elems.second)
          {
            return;
          }
        }
      }

      usingInner_ = false;
      InnerIteratorType & out = iterTT_->outer();
      out.first();
      if(!out.done())
      {
        pair < ElementPllXIF_t *, int > p = out.item ().accessPllX ().accessOuterPllX () ;
        pair < HElementType * , HBndSegType * > elems(0,0);
        p.first->getAttachedElement(elems);

        assert( elems.second );
        it_ = &out;
        return ;
      }

      it_ = 0;
    }

    virtual void checkLeafEntity ()
    {
      if(it_)
      {
        if(!it_->done())
        {
          val_t & el = item();
          HBndSegType * pll = el.second;
          assert( pll );

          // this occurs if internal element is leaf but the corresponding
          // ghost is not leaf, we have to go next
          if ( ! pll->isLeafEntity() ) next();
        }
      }
    }

  public:
    int size  ()    // ???? gives size only of small part of ghost cells ????
    {
      // if no iterator then size is zero
      // which can happen in the case of parallel grid with 1 processor
      if(!it_)
      {
        return 0;
      }
      return it_->size();
    }

    // go next ghost
    void next ()
    {
      if(it_)
      {
        // if not done increment
        if( !it_->done() ) it_->next();

        // if now done, create new iterator
        if( it_->done() ) createIterator();

        checkLeafEntity();
      }
    }

    void first()
    {
      link_ = -1;
      usingInner_ = false;
      // create iterator calls also first of iterators
      createIterator();
      checkLeafEntity();
      if( it_ ) assert( !it_->done());
    }

    int done () const
    {
      assert( (link_ >= nl_) ? (it_ == 0) : 1 );
      return ((link_ >= nl_ || !it_ ) ? 1 : 0);
    }

    val_t & item () const
    {
      assert(it_);
      pair < ElementPllXIF_t *, int > p = it_->item ().accessPllX ().accessOuterPllX () ;
      pair < HElementType  * , HBndSegType * > p2;
      p.first->getAttachedElement(p2);
      assert(p2.second);
      elem_.second = p2.second;
      return elem_;
    }

  }; // end ALU3dGridGhostIterator


  // the leaf ghost partition iterator
  template <>
  class ALU3dGridLeafIteratorWrapper<0,Dune::Ghost_Partition>
    : public ALU3dGridGhostIterator
  {
  protected:
    typedef LeafLevelIteratorTTProxy<1> IteratorType;
    IteratorType * newIterator()
    {
      return new IteratorType ( this->gitter_, this->link_ );
    }

    void checkLeafEntity ()
    {
      if(this->it_)
      {
        if(! this->it_->done())
        {
          val_t & el = this->item();
          HBndSegType * pll = el.second;
          assert( pll );

          // this occurs if internal element is leaf but the corresponding
          // ghost is not leaf, we have to go next
          if ( ! pll->isLeafEntity() ) this->next();
        }
      }
    }

  public:
    template <class GridImp>
    ALU3dGridLeafIteratorWrapper(const GridImp & grid, int level , const int nlinks )
      : ALU3dGridGhostIterator(grid,level,nlinks) {}

    ALU3dGridLeafIteratorWrapper(const ALU3dGridLeafIteratorWrapper & org)
      : ALU3dGridGhostIterator(org) {}
  };

  // the level ghost partition iterator
  template <>
  class ALU3dGridLevelIteratorWrapper<0,Dune::Ghost_Partition>
    : public ALU3dGridGhostIterator
  {
    const int level_;
    const int mxl_;
  protected:
    typedef LeafLevelIteratorTTProxy<1> IteratorType;
    IteratorType * newIterator()
    {
      // create new level Iterator Proxy
      return new IteratorType ( this->gitter_, this->link_ , level_ );
    }

    // for level iterators don't check leaf entity
    void checkLeafEntity ()
    {
      if(this->it_)
      {
        if(! this->it_->done())
        {
          val_t & el = this->item();

          assert( el.second );
          HBndSegType & pll = *(el.second);

          // this occurs if internal element is leaf but the corresponding
          // ghost is not leaf, we have to go next if level of ghost is not
          // our level
          if ( ! pll.down() )
          {
            if( pll.ghostLevel() != level_ ) this->next();
          }
        }
      }
    }

  public:
    template <class GridImp>
    ALU3dGridLevelIteratorWrapper(const GridImp & grid,
                                  const VertexListType & , int level , const int nlinks )
      : ALU3dGridGhostIterator(grid,level,nlinks) , level_(level) ,
        mxl_(grid.maxLevel()){}

    template <class GridImp>
    ALU3dGridLevelIteratorWrapper(const GridImp & grid, int level , const int nlinks )
      : ALU3dGridGhostIterator(grid,level,nlinks) , level_(level),
        mxl_(grid.maxLevel()) {}

    ALU3dGridLevelIteratorWrapper(const ALU3dGridLevelIteratorWrapper & org)
      : ALU3dGridGhostIterator(org) , level_(org.level_) , mxl_(org.mxl_){}
  };

  template <int codim>
  class ALU3dGridGhostIteratorHigherCodim
    : public IteratorWrapperInterface < typename IteratorElType<codim>::val_t >
  {
    typedef typename IteratorElType<codim>::ElType ElType;
    typedef typename IteratorElType<codim>::val_t val_t;

    template<Dune :: ALU3dGridElementType elType, int cd>
    struct SelectVector
    {};

    template<Dune :: ALU3dGridElementType elType>
    struct SelectVector<elType,1>
    {
      typedef typename Dune :: ALU3dImplTraits<elType>::GEOElementType GEOElementType;
      static const vector<int> & getNotOnItemVector(int face)
      {
        return GEOElementType :: facesNotOnFace( face );
      }
    };

    template<Dune :: ALU3dGridElementType elType>
    struct SelectVector<elType,2>
    {
      typedef typename Dune :: ALU3dImplTraits<elType>::GEOElementType GEOElementType;
      static const vector<int> & getNotOnItemVector(int face)
      {
        return GEOElementType :: edgesNotOnFace( face );
      }
    };

    template<Dune :: ALU3dGridElementType elType>
    struct SelectVector<elType,3>
    {
      typedef typename Dune :: ALU3dImplTraits<elType>::GEOElementType GEOElementType;
      static const vector<int> & getNotOnItemVector(int face)
      {
        return GEOElementType :: verticesNotOnFace( face );
      }
    };

    template<class GridImp, int cd>
    struct GetItem;

    template<class GridImp>
    struct GetItem<GridImp,1>
    {
      static const Dune::ALU3dGridElementType elType = GridImp::elementType ;

      typedef typename Dune::ALU3dImplTraits<elType>::GEOElementType GEOElementType;
      typedef typename IteratorElType<1>::ElType ItemType;

      static ItemType * getItemFromEl(GEOTetraElementType & el, int i)
      {
        return el.myhface3(i);
      }

      static ItemType * getItemFromEl(GEOHexaElementType & el, int i)
      {
        return el.myhface4(i);
      }

      static ItemType * getItem(HElementType & el, int i)
      {
        return getItemFromEl(static_cast<GEOElementType &> (el),i);
      }

      static int numItems()
      {
        return Dune :: EntityCount<elType>::numFaces;
      }
    };

    template<class GridImp>
    struct GetItem<GridImp,2>
    {
      static const Dune::ALU3dGridElementType elType = GridImp::elementType ;
      typedef typename Dune::ALU3dImplTraits<elType>::GEOElementType GEOElementType;
      typedef typename IteratorElType<2>::ElType ItemType;
      static ItemType * getItem(HElementType & el, int i)
      {
        return static_cast<GEOElementType &> (el).myhedge1(i);
      }

      static int numItems()
      {
        return Dune :: EntityCount<elType>::numEdges;
      }
    };

    template<class GridImp>
    struct GetItem<GridImp,3>
    {
      static const Dune::ALU3dGridElementType elType = GridImp::elementType ;
      typedef typename Dune::ALU3dImplTraits<elType>::GEOElementType GEOElementType;
      typedef typename IteratorElType<3>::ElType ItemType;
      static ItemType * getItem(HElementType & el, int i)
      {
        return static_cast<GEOElementType &> (el).myvertex(i);
      }

      static int numItems()
      {
        return Dune :: EntityCount<elType>::numVertices;
      }
    };

    typedef ElType * getItemFunc_t (HElementType & el, int i);

  private:
    typedef Dune :: ALU3dGridItemListType GhostItemListType;
    mutable GhostItemListType & ghList_;
    typedef typename GhostItemListType :: IteratorType IteratorType;
    IteratorType curr_;
    IteratorType end_;
    mutable val_t elem_;
    mutable int count_;
  public:
    template <class GhostElementIteratorImp, class GridImp>
    ALU3dGridGhostIteratorHigherCodim(GhostElementIteratorImp *, const GridImp & grid,
                                      int level , const int nlinks, GhostItemListType & ghList)
      : ghList_( ghList )
        , elem_(0,0)
        , count_(0)
    {
      if( ! ghList_.up2Date() )
      {
        GhostElementIteratorImp ghostIter(grid,level,nlinks);
        updateGhostList(grid,ghostIter,ghList_);
      }
    }

    ALU3dGridGhostIteratorHigherCodim(const ALU3dGridGhostIteratorHigherCodim & org)
      : ghList_( org.ghList_ )
        , elem_(org.elem_)
        , count_(org.count_)
    {}

    int size  () { return ghList_.getItemList().size(); }
    void first() { count_ = 0; }
    void next () { ++count_; }
    int done () const { return (count_ >= (int) ghList_.getItemList().size() ? 1 : 0); }
    val_t & item () const
    {
      assert( ! done() );
      void * item = ghList_.getItemList()[count_];
      elem_.first = ((ElType * ) item);
      assert( elem_.first );
      return elem_;
    }

  protected:
    template <class GridImp, class GhostElementIteratorImp>
    void updateGhostList(const GridImp & grid, GhostElementIteratorImp & ghostIter, GhostItemListType & ghList)
    {
      int count = 0;
      for( ghostIter.first(); !ghostIter.done(); ghostIter.next() )
      {
        ++count;
      }

      const int numItems = SelectVector<GridImp::elementType,codim>::getNotOnItemVector(0).size();
      const int maxSize = numItems * count;

      ghList.getItemList().reserve(maxSize);
      ghList.getItemList().resize(0);
      map< int , int > visited;

      for( ghostIter.first(); !ghostIter.done(); ghostIter.next() )
      {
        GhostPairType ghPair = ghostIter.item().second->getGhost();
        const vector<int> & notOnFace = SelectVector<GridImp::elementType,codim>::
                                        getNotOnItemVector(ghPair.second);
        for(int i=0; i<numItems; ++i)
        {
          ElType * item = GetItem<GridImp,codim>::getItem( *(ghPair.first) , notOnFace[i] );
          int idx = item->getIndex();
          if( visited.find(idx) == visited.end() )
          {
            ghList.getItemList().push_back( (void *) item );
            visited[idx] = 1;
          }
        }
      }
      ghList.markAsUp2Date();
    }
  };

  // the leaf ghost partition iterator
  template <>
  class ALU3dGridLeafIteratorWrapper<1,Dune::Ghost_Partition>
    : public ALU3dGridGhostIteratorHigherCodim<1>
  {
    enum { codim = 1 };
    typedef ALU3dGridLeafIteratorWrapper<0,Dune::Ghost_Partition> GhostElementIteratorType;
  public:
    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLeafList(codim)) {}

    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the leaf ghost partition iterator
  template <>
  class ALU3dGridLeafIteratorWrapper<2,Dune::Ghost_Partition>
    : public ALU3dGridGhostIteratorHigherCodim<2>
  {
    enum { codim = 2 };
    typedef ALU3dGridLeafIteratorWrapper<0,Dune::Ghost_Partition> GhostElementIteratorType;
  public:
    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLeafList(codim)) {}

    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the leaf ghost partition iterator
  template <>
  class ALU3dGridLeafIteratorWrapper<3,Dune::Ghost_Partition>
    : public ALU3dGridGhostIteratorHigherCodim<3>
  {
    enum { codim = 3 };
    typedef ALU3dGridLeafIteratorWrapper<0,Dune::Ghost_Partition> GhostElementIteratorType;
  public:
    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLeafList(codim)) {}

    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the level ghost partition iterator
  template <>
  class ALU3dGridLevelIteratorWrapper<1,Dune::Ghost_Partition>
    : public ALU3dGridGhostIteratorHigherCodim<1>
  {
    enum { codim = 1 };
    typedef ALU3dGridLevelIteratorWrapper<0,Dune::Ghost_Partition> GhostElementIteratorType;
  public:
    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, const VertexListType & , int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLevelList(codim,level)) {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the level ghost partition iterator
  template <>
  class ALU3dGridLevelIteratorWrapper<2,Dune::Ghost_Partition>
    : public ALU3dGridGhostIteratorHigherCodim<2>
  {
    enum { codim = 2 };
    typedef ALU3dGridLevelIteratorWrapper<0,Dune::Ghost_Partition> GhostElementIteratorType;
  public:
    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, const VertexListType & , int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLevelList(codim,level)) {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the level ghost partition iterator
  template <>
  class ALU3dGridLevelIteratorWrapper<3,Dune::Ghost_Partition>
    : public ALU3dGridGhostIteratorHigherCodim<3>
  {
    enum { codim = 3 };
    typedef ALU3dGridLevelIteratorWrapper<0,Dune::Ghost_Partition> GhostElementIteratorType;
  public:
    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, const VertexListType & , int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLevelList(codim,level)) {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the all partition iterator
  template <>
  class ALU3dGridLeafIteratorWrapper<0,Dune::All_Partition>
    : public IteratorWrapperInterface < IteratorElType<0>::val_t >
  {
    enum { codim = 0 };
    typedef ALU3dGridLeafIteratorWrapper<codim,Dune::InteriorBorder_Partition> InteriorIteratorType;
    typedef ALU3dGridLeafIteratorWrapper<codim,Dune::Ghost_Partition> GhostIteratorType;

  public:
    typedef IteratorElType<codim>::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType , val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, level, nlinks ) ,
                GhostIteratorType    ( grid, level, nlinks ) )
    {}

    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper & org )
      : iter_ (org.iter_) {}

    int size  () { return iter_.size(); }
    void next () { iter_.next(); }
    void first() { iter_.first(); }
    int done () const {return iter_.done(); }
    val_t & item () const { assert( ! done() ); return iter_.item(); }
  };

  // the all partition iterator
  template <>
  class ALU3dGridLeafIteratorWrapper<1,Dune::All_Partition>
    : public IteratorWrapperInterface < IteratorElType<1>::val_t >
  {
    enum { codim = 1 };
    typedef ALU3dGridLeafIteratorWrapper<codim,Dune::InteriorBorder_Partition> InteriorIteratorType;
    typedef ALU3dGridLeafIteratorWrapper<codim,Dune::Ghost_Partition> GhostIteratorType;

  public:
    typedef IteratorElType<codim>::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType , val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, level, nlinks ) ,
                GhostIteratorType    ( grid, level, nlinks ) )
    {}

    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper & org )
      : iter_ (org.iter_) {}

    int size  () { return iter_.size(); }
    void next () { iter_.next(); }
    void first() { iter_.first(); }
    int done () const {return iter_.done(); }
    val_t & item () const { assert( ! done() ); return iter_.item(); }
  };

  // the all partition iterator
  template <>
  class ALU3dGridLeafIteratorWrapper<2,Dune::All_Partition>
    : public IteratorWrapperInterface < IteratorElType<2>::val_t >
  {
    enum { codim = 2 };
    typedef ALU3dGridLeafIteratorWrapper<codim,Dune::InteriorBorder_Partition> InteriorIteratorType;
    typedef ALU3dGridLeafIteratorWrapper<codim,Dune::Ghost_Partition> GhostIteratorType;

  public:
    typedef IteratorElType<codim>::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType , val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, level, nlinks ) ,
                GhostIteratorType    ( grid, level, nlinks ) )
    {}

    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper & org )
      : iter_ (org.iter_) {}

    int size  () { return iter_.size(); }
    void next () { iter_.next(); }
    void first() { iter_.first(); }
    int done () const {return iter_.done(); }
    val_t & item () const { assert( ! done() ); return iter_.item(); }
  };

  // the all partition iterator
  template <>
  class ALU3dGridLeafIteratorWrapper<3,Dune::All_Partition>
    : public IteratorWrapperInterface < IteratorElType<3>::val_t >
  {
    enum { codim = 3 };
    typedef ALU3dGridLeafIteratorWrapper<codim,Dune::InteriorBorder_Partition> InteriorIteratorType;
    typedef ALU3dGridLeafIteratorWrapper<codim,Dune::Ghost_Partition> GhostIteratorType;

  public:
    typedef IteratorElType<codim>::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType , val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, level, nlinks ) ,
                GhostIteratorType    ( grid, level, nlinks ) )
    {}

    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper & org )
      : iter_ (org.iter_) {}

    int size  () { return iter_.size(); }
    void next () { iter_.next(); }
    void first() { iter_.first(); }
    int done () const {return iter_.done(); }
    val_t & item () const { assert( ! done() ); return iter_.item(); }
  };

  // the all partition iterator
  template <>
  class ALU3dGridLevelIteratorWrapper<0,Dune::All_Partition>
    : public IteratorWrapperInterface< LeafValType >
  {
    typedef ALU3dGridLevelIteratorWrapper<0,Dune::InteriorBorder_Partition> InteriorIteratorType;
    typedef ALU3dGridLevelIteratorWrapper<0,Dune::Ghost_Partition> GhostIteratorType;

  public:
    typedef LeafValType val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType , val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, const VertexListType & vxList, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, vxList, level, nlinks ) ,
                GhostIteratorType    ( grid, vxList, level, nlinks ) )
    {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org)
      : iter_(org.iter_) {}

    int size  () { return iter_.size(); }
    void next () { iter_.next(); }
    void first() { iter_.first(); }
    int done () const {return iter_.done(); }
    val_t & item () const { assert( ! done() ); return iter_.item(); }
  };
#endif // end ALU3DGRID_PARALLEL

  // placed here because we need ALU3dGridLevelIteratorWrapper<0,Dune::All_Partition> here
  // the edge level iterator
  template <PartitionIteratorType pitype>
  class ALU3dGridLevelIteratorWrapper<2,pitype>
    : public IteratorWrapperInterface< typename IteratorElType<2>::val_t >
  {
  public:
    typedef ALUHElementType<2>::ElementType ElType;
    typedef ALU3DSPACE GEOEdgeT GEOEdgeType;

    typedef typename IteratorElType<2>::val_t val_t;
  private:
    mutable val_t elem_;
    const int level_;

    typedef Dune :: ALU3dGridItemListType ItemListType;
    mutable ItemListType & edgeList_;

    size_t count_ ;
    bool maxLevel_;

  public:
    // constructor creating iterator
    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid,
                                   const VertexListType & vxList , int level , const int nlinks )
      : elem_(0,0)
        , level_(level)
        , edgeList_( grid.getEdgeList(level) )
        , count_(0)
    {
      if( ! edgeList_.up2Date() )
        updateEdgeList(grid,vxList,level,nlinks);
    }

    // copy constructor
    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : elem_(org.elem_)
        , level_(org.level_)
        , edgeList_( org.edgeList_ )
        , count_(org.count_)
    {}

    int size  () { return edgeList_.getItemList().size(); }
    void next ()
    {
      ++count_;
    }

    void first()
    {
      count_ = 0;
    }

    int done () const { return ((count_ >= edgeList_.size()) ? 1 : 0); }

    val_t & item () const
    {
      assert( ! done () );
      elem_.first = ( (ElType *) edgeList_.getItemList()[count_]);

      assert( elem_.first );
      return elem_;
    }
  private:
    template <class GridImp>
    void updateEdgeList(const GridImp & grid, const VertexListType & vxList, int level, int nlinks)
    {
      typedef ALU3dGridLevelIteratorWrapper<0,Dune::All_Partition> ElementLevelIterator;
      typedef typename ElementLevelIterator :: val_t el_val_t;
      ElementLevelIterator iter(grid,vxList,level,nlinks);

      edgeList_.getItemList().resize(0);
      map < int , int > visited;

      for( iter.first(); ! iter.done(); iter.next() )
      {
        typedef typename Dune :: ALU3dImplTraits<GridImp::elementType>::GEOElementType GEOElementType;
        enum { numEdges = Dune :: EntityCount<GridImp::elementType>::numEdges };

        GEOElementType * elem = 0;
        el_val_t & item = iter.item();

        if( item.first )
          elem = static_cast<GEOElementType *> (item.first);
        else if( item.second )
          elem = static_cast<GEOElementType *> (item.second->getGhost().first);

        assert( elem );
        for(int e=0; e<numEdges; ++e)
        {
          ElType * edge = elem->myhedge1(e);
          if( edge->isGhost() ) continue;

          int idx = edge->getIndex();
          if( visited.find(idx) == visited.end() )
          {
            edgeList_.getItemList().push_back( (void *) edge );
            visited[idx] = 1;
          }
        }
      }
      edgeList_.markAsUp2Date();
    }
  };

#if ALU3DGRID_PARALLEL
  // the all partition iterator
  template <>
  class ALU3dGridLevelIteratorWrapper<1,Dune::All_Partition>
    : public IteratorWrapperInterface < IteratorElType<1>::val_t >
  {
    enum { codim = 1 };
    typedef ALU3dGridLevelIteratorWrapper<codim,Dune::InteriorBorder_Partition> InteriorIteratorType;
    typedef ALU3dGridLevelIteratorWrapper<codim,Dune::Ghost_Partition> GhostIteratorType;

  public:
    typedef IteratorElType<codim>::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType , val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, const VertexListType & vxList, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, vxList, level, nlinks ) ,
                GhostIteratorType    ( grid, vxList, level, nlinks ) )
    {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : iter_ (org.iter_) {}

    int size  () { return iter_.size(); }
    void next () { iter_.next(); }
    void first() { iter_.first(); }
    int done () const {return iter_.done(); }
    val_t & item () const { assert( ! done() ); return iter_.item(); }
  };

  // the all partition iterator
  template <>
  class ALU3dGridLevelIteratorWrapper<2,Dune::All_Partition>
    : public IteratorWrapperInterface < IteratorElType<2>::val_t >
  {
    enum { codim = 2 };
    typedef ALU3dGridLevelIteratorWrapper<codim,Dune::InteriorBorder_Partition> InteriorIteratorType;
    typedef ALU3dGridLevelIteratorWrapper<codim,Dune::Ghost_Partition> GhostIteratorType;

  public:
    typedef IteratorElType<codim>::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType , val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, VertexListType & vxList, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, vxList, level, nlinks ) ,
                GhostIteratorType    ( grid, vxList, level, nlinks ) )
    {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : iter_ (org.iter_) {}

    int size  () { return iter_.size(); }
    void next () { iter_.next(); }
    void first() { iter_.first(); }
    int done () const {return iter_.done(); }
    val_t & item () const { assert( ! done() ); return iter_.item(); }
  };

  // the all partition iterator
  template <>
  class ALU3dGridLevelIteratorWrapper<3,Dune::All_Partition>
    : public IteratorWrapperInterface < IteratorElType<3>::val_t >
  {
    enum { codim = 3 };
    typedef ALU3dGridLevelIteratorWrapper<codim,Dune::InteriorBorder_Partition> InteriorIteratorType;
    typedef ALU3dGridLevelIteratorWrapper<codim,Dune::Ghost_Partition> GhostIteratorType;

  public:
    typedef IteratorElType<codim>::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType , val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, VertexListType & vxList, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, vxList, level, nlinks ) ,
                GhostIteratorType    ( grid, vxList, level, nlinks ) )
    {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : iter_ (org.iter_) {}

    int size  () { return iter_.size(); }
    void next () { iter_.next(); }
    void first() { iter_.first(); }
    int done () const {return iter_.done(); }
    val_t & item () const { assert( ! done() ); return iter_.item(); }
  };
#endif // end if ALU3DGRID_PARALLEL

  typedef PureElementLeafIterator < GitterType::helement_STI > BSLeafIteratorMaxLevel;

} //end namespace ALU3dGrid


namespace Dune {
  // Forward declarations
  template<int cd, int dim, class GridImp>
  class ALU3dGridEntity;
  template<int cd, PartitionIteratorType pitype, class GridImp >
  class ALU3dGridLevelIterator;
  template<int cd, class GridImp >
  class ALU3dGridEntityPointer;
  template<int mydim, int coorddim, class GridImp>
  class ALU3dGridGeometry;
  template<class GridImp>
  class ALU3dGridHierarchicIterator;
  template<class GridImp>
  class ALU3dGridIntersectionIterator;
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class ALU3dGridLeafIterator;
  template<int dim, int dimworld, ALU3dGridElementType elType>
  class ALU3dGrid;
  template <ALU3dGridElementType type>
  class ALU3dGridFaceInfo;
  template <ALU3dGridElementType elType>
  class ALU3dGridGeometricFaceInfo;

  //**********************************************************************
  //
  // --ALU3dGridIntersectionIterator
  // --IntersectionIterator
  /*!
     Mesh entities of codimension 0 ("elements") allow to visit all neighbors, wh
     a neighbor is an entity of codimension 0 which has a common entity of codimens
     These neighbors are accessed via a IntersectionIterator. This allows the implement
     non-matching meshes. The number of neigbors may be different from the number o
     of an element!
   */
  template<class GridImp>
  class ALU3dGridIntersectionIterator :
    public IntersectionIteratorDefaultImplementation <GridImp,ALU3dGridIntersectionIterator>
  {
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    typedef ALU3dImplTraits<GridImp::elementType> ImplTraits;
    typedef typename ImplTraits::GEOElementType GEOElementType;
    typedef typename ImplTraits::IMPLElementType IMPLElementType;
    typedef typename ImplTraits::GEOFaceType GEOFaceType;
    typedef typename ImplTraits::NeighbourPairType NeighbourPairType;
    typedef typename ImplTraits::PLLBndFaceType PLLBndFaceType;
    typedef typename ImplTraits::BNDFaceType BNDFaceType;

    typedef ALU3dGridFaceInfo<GridImp::elementType> FaceInfoType;
    typedef typename std::auto_ptr<FaceInfoType> FaceInfoPointer;

    typedef typename SelectType<
        SameType<Int2Type<tetra>, Int2Type<GridImp::elementType> >::value,
        ALU3dGridGeometricFaceInfoTetra,
        ALU3dGridGeometricFaceInfoHexa
        >::Type GeometryInfoType;

    typedef ElementTopologyMapping<GridImp::elementType> ElementTopo;
    typedef FaceTopologyMapping<GridImp::elementType> FaceTopo;

    enum { numFaces = EntityCount<GridImp::elementType>::numFaces };
    enum { numVerticesPerFace =
             EntityCount<GridImp::elementType>::numVerticesPerFace };
    enum { numVertices = EntityCount<GridImp::elementType>::numVertices };

    friend class ALU3dGridEntity<0,dim,GridImp>;
    friend class IntersectionIteratorWrapper<GridImp>;

  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    typedef ALU3dGridGeometry<dim-1,dimworld,GridImp> GeometryImp;
    typedef MakeableInterfaceObject<Geometry> GeometryObject;

    typedef FieldVector<alu3d_ctype, dimworld> NormalType;
    typedef ALU3dGridEntityPointer<0,GridImp> EntityPointer;

    //! The default Constructor , level tells on which level we want
    //! neighbours
    ALU3dGridIntersectionIterator(const GridImp & grid,
                                  ALU3DSPACE HElementType *el,
                                  int wLevel,bool end=false);

    ALU3dGridIntersectionIterator(const GridImp & grid,int wLevel);

    //! The copy constructor
    ALU3dGridIntersectionIterator(const ALU3dGridIntersectionIterator<GridImp> & org);

    //! assignment of iterators
    void assign(const ALU3dGridIntersectionIterator<GridImp> & org);

    //! The copy constructor
    bool equals (const ALU3dGridIntersectionIterator<GridImp> & i) const;

    //! increment iterator
    void increment ();

    //! access neighbor
    EntityPointer outside() const;

    //! access entity where iteration started
    EntityPointer inside() const;

    //! return true if intersection is with boundary.
    bool boundary () const;

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const;

    //! return true if across the edge an neighbor on this level exists
    bool levelNeighbor () const;

    //! return true if across the edge an neighbor on leaf level exists
    bool leafNeighbor () const;

    //! return information about the Boundary
    int boundaryId () const;

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    const LocalGeometry & intersectionSelfLocal () const;

    //! intersection of codimension 1 of this neighbor with element where
    //!  iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where
    //! iteration started.
    const Geometry & intersectionGlobal () const;

    //! local number of codim 1 entity in self where intersection is contained
    //!  in
    int numberInSelf () const;

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    const LocalGeometry & intersectionNeighborLocal () const;

    //! local number of codim 1 entity in neighbor where intersection is
    //! contained
    int numberInNeighbor () const;

    int twistInSelf() const;
    int twistInNeighbor() const;

    //! return unit outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    NormalType & unitOuterNormal (const FieldVector<alu3d_ctype, dim-1>& local) const ;

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    NormalType & outerNormal (const FieldVector<alu3d_ctype, dim-1>& local) const;

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    NormalType & integrationOuterNormal (const FieldVector<alu3d_ctype, dim-1>& local) const;

    //! return level of iterator
    int level () const;

  private:
    // set interator to end iterator
    void done () ;

    void outputElementInfo() const;

    void outputFaceInfo() const;

    template <typename T>
    void printToScreen(int duneIdx, int aluIdx,
                       const T& info) const;

    void printToScreen(int duneIdx, int aluIdx) const;

    // used in printToScreen
    NormalType convert2FV(const alu3d_ctype (&p)[3]) const;

    // reset IntersectionIterator to first neighbour
    void setFirstItem(const ALU3DSPACE HElementType & elem, int wLevel);

    // reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    // set new face
    void setNewFace(const GEOFaceType& newFace);

    // is there a refined element at the outer side of the face which needs to be considered when incrementing the iterator?
    bool canGoDown(const GEOFaceType& nextFace) const;

    void buildLocalGeometries() const;

    void buildGlobalGeometry() const;

    // get the face corresponding to the index
    const typename ALU3dImplTraits<tetra>::GEOFaceType*
    getFace(const ALU3DSPACE GEOTetraElementType& elem, int index) const;

    const typename ALU3dImplTraits<hexa>::GEOFaceType*
    getFace(const ALU3DSPACE GEOHexaElementType& elem, int index) const;

    //! structure containing the topological and geometrical information about
    //! the face which the iterator points to
    mutable FaceInfoType connector_;
    mutable GeometryInfoType geoProvider_; // need to initialise

    // reference to grid
    const GridImp & grid_;

    //! current element from which we started the intersection iterator
    const IMPLElementType* item_;

    mutable int nFaces_;
    mutable int walkLevel_;
    mutable int index_;

    mutable bool generatedGlobalGeometry_;
    mutable bool generatedLocalGeometries_;

    mutable GeometryObject intersectionGlobal_;
    mutable GeometryImp &  intersectionGlobalImp_;
    mutable GeometryObject intersectionSelfLocal_;
    mutable GeometryImp &  intersectionSelfLocalImp_;
    mutable GeometryObject intersectionNeighborLocal_;
    mutable GeometryImp &  intersectionNeighborLocalImp_;

    // unit outer normal
    mutable NormalType unitOuterNormal_;

    // true if end iterator
    bool done_;

    bool levelNeighbor_;
    bool leafNeighbor_;
    bool goneDown_;
    bool isLeafItem_;
  };

  //////////////////////////////////////////////////////////////////////////////
  //
  //  --IterationImpl
  //
  //////////////////////////////////////////////////////////////////////////////
  template <class InternalIteratorType >
  struct ALU3dGridTreeIterator
  {
    typedef typename InternalIteratorType :: val_t val_t;
  protected:
    // set iterator to first item
    template <class IteratorImp>
    void firstItem(IteratorImp & it)
    {
      InternalIteratorType & iter = it.internalIterator();
      iter.first();
      if( ! iter.done() )
      {
        assert( iter.size() > 0 );
        setItem(it,iter);
      }
      else
      {
        it.removeIter();
      }
    }

    // set the iterators entity to actual item
    template <class IteratorImp>
    void setItem (IteratorImp & it, InternalIteratorType & iter)
    {
      val_t & item = iter.item();
      assert( item.first || item.second );
      if( item.first )
        it.updateEntityPointer( item.first );
      else
        it.updateGhostPointer( *item.second );
    }

    // increment iterator
    template <class IteratorImp>
    void incrementIterator(IteratorImp & it)
    {
      // if iter_ is zero, then end iterator
      InternalIteratorType & iter = it.internalIterator();

      iter.next();

      if(iter.done())
      {
        it.removeIter();
        return ;
      }

      setItem(it,iter);
      return ;
    }
  };

  //**********************************************************************
  //
  // --ALU3dGridLevelIterator
  // --LevelIterator
  /*!
     Enables iteration over all entities of a given codimension and level of a grid.
   */
  template<int cd, PartitionIteratorType pitype, class GridImp>
  class ALU3dGridLevelIterator :
    public ALU3dGridEntityPointer <cd,GridImp> ,
    public LevelIteratorDefaultImplementation <cd,pitype,GridImp,ALU3dGridLevelIterator> ,
    public ALU3dGridTreeIterator< ALU3DSPACE ALU3dGridLevelIteratorWrapper<cd,pitype> >
  {
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    friend class ALU3dGridEntity<3,dim,GridImp>;
    friend class ALU3dGridEntity<2,dim,GridImp>;
    friend class ALU3dGridEntity<1,dim,GridImp>;
    friend class ALU3dGridEntity<0,dim,GridImp>;
    friend class ALU3dGrid < dim , dimworld, GridImp::elementType >;

    friend class ALU3dGridTreeIterator< ALU3DSPACE ALU3dGridLevelIteratorWrapper<cd,pitype> > ;
  public:
    typedef typename GridImp::template Codim<cd>::Entity Entity;
    typedef ALU3DSPACE VertexListType VertexListType;

    //! typedef of my type
    typedef ALU3dGridLevelIterator<cd,pitype,GridImp> ThisType;
    // the wrapper for the original iterator of the ALU3dGrid
    typedef typename ALU3DSPACE ALU3dGridLevelIteratorWrapper<cd,pitype> IteratorType;
    typedef IteratorType InternalIteratorType;
    typedef typename ALU3DSPACE IteratorElType<cd>::val_t val_t;

    //! Constructor for begin iterator
    ALU3dGridLevelIterator(const GridImp & grid, VertexListType & vxList, int level);

    //! Constructor for end iterator
    ALU3dGridLevelIterator(const GridImp & grid, int level);

    //! Constructor
    ALU3dGridLevelIterator(const ThisType & org);

    // destructor
    ~ALU3dGridLevelIterator();

    //! prefix increment
    void increment ();

    //! dereference Entity, faster then the entity pointersmethod
    Entity & dereference () const;

    //! assignment of iterators
    ThisType & operator = (const ThisType & org);
  private:

    // actual level
    int level_;

    // the wrapper for the original iterator of the ALU3dGrid
    typedef typename ALU3DSPACE ALU3dGridLevelIteratorWrapper<cd,pitype> IteratorType;

    // the internal iterator
    IteratorType * iter_ ;

    // deletes iter_
    void removeIter ();
    IteratorType & internalIterator ()
    {
      assert( iter_ );
      return *iter_;
    }
  };

  //********************************************************************
  //
  //  --ALU3dGridLeafIterator
  //  --LeafIterator
  //
  //********************************************************************
  //! Leaf iterator
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  class ALU3dGridLeafIterator :
    public LeafIteratorDefaultImplementation<cdim, pitype, GridImp, ALU3dGridLeafIterator>,
    public ALU3dGridEntityPointer<cdim,GridImp> ,
    public ALU3dGridTreeIterator< ALU3DSPACE ALU3dGridLeafIteratorWrapper<cdim, pitype> >
  {
    enum { dim = GridImp :: dimension };

    friend class ALU3dGridEntity<cdim,dim,GridImp>;
    enum { codim = cdim };

    friend class ALU3dGridTreeIterator< ALU3DSPACE ALU3dGridLeafIteratorWrapper<cdim, pitype> > ;

  public:
    typedef typename GridImp::template Codim<cdim>::Entity Entity;

    // the wrapper for the original iterator of the ALU3dGrid
    typedef typename ALU3DSPACE ALU3dGridLeafIteratorWrapper<cdim, pitype> IteratorType;
    typedef IteratorType InternalIteratorType;
    typedef typename ALU3DSPACE IteratorElType<cdim>::val_t val_t;

    typedef ALU3dGridLeafIterator<cdim, pitype, GridImp> ThisType;

    //! Constructor for end iterators
    ALU3dGridLeafIterator(const GridImp & grid, int level);

    //! Constructor for begin Iterators
    ALU3dGridLeafIterator(const GridImp & grid, int level , bool isBegin);

    //! copy Constructor
    ALU3dGridLeafIterator(const ThisType & org);

    //! destructor deleting real iterator
    ~ALU3dGridLeafIterator();

    //! prefix increment
    void increment ();

    //! dereference Entity, faster then the entity pointersmethod
    Entity & dereference () const;

    //! assignment of iterators
    ThisType & operator = (const ThisType & org);

  private:
    // the internal iterator
    IteratorType * iter_;

    //! do assignment
    void assign (const ThisType & org, bool clone = true );

    // deletes iter_
    void removeIter () ;

    // return reference to iter_
    InternalIteratorType & internalIterator ()
    {
      assert( iter_ );
      return *iter_;
    }
  };

  // - HierarchicIteraor
  // --HierarchicIterator
  template<class GridImp>
  class ALU3dGridHierarchicIterator :
    public ALU3dGridEntityPointer<0,GridImp> ,
    public HierarchicIteratorDefaultImplementation <GridImp,ALU3dGridHierarchicIterator>
  {
    enum { dim = GridImp::dimension };
    typedef ALU3dGridHierarchicIterator<GridImp> ThisType;
  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::ctype ctype;

    //! the normal Constructor
    ALU3dGridHierarchicIterator(const GridImp &grid,
                                const ALU3DSPACE HElementType & elem, int maxlevel, bool end=false);

    //! the normal Constructor
    ALU3dGridHierarchicIterator(const ALU3dGridHierarchicIterator<GridImp> &org);

    //! increment
    void increment();

    //! dereference Entity, faster then the entity pointersmethod
    Entity & dereference () const;

    //! the assignment operator
    ThisType & operator = (const ThisType & org);

  private:
    // go to next valid element
    ALU3DSPACE HElementType * goNextElement (ALU3DSPACE HElementType * oldEl);

    //! element from where we started
    const ALU3DSPACE HElementType * elem_;

    //! maximal level to go down
    int maxlevel_;
  };


} // end namespace Dune

#include "iterator_imp.cc"

#endif
