// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DITERATORS_HH
#define DUNE_ALU3DITERATORS_HH

// System includes

// Dune includes
#include <dune/grid/common/grid.hh>

// Local includes
#include "alu3dinclude.hh"
#include "topology.hh"

namespace ALUGridSpace
{

  //*************************************************************
  //  definition of original LeafIterators of ALUGrid
  //
  // default is element (codim = 0)
  template< int codim, class Comm >
  struct BSMacroIterator
  {
    typedef typename Dune::ALU3dBasicImplTraits< Comm >::HElementType HElementType;
    typedef typename AccessIterator< HElementType >::Handle IteratorType;
  };

  //******************************************************************
  //  LevelIterators
  //******************************************************************
  template< int codim, class Comm >
  struct ALUHElementType;

  template< class Comm >
  struct ALUHElementType< 0, Comm >
  {
    typedef typename Dune::ALU3dBasicImplTraits< Comm >::HElementType ElementType;
  };

  template< class Comm >
  struct ALUHElementType< 1, Comm >
  {
    typedef typename Dune::ALU3dBasicImplTraits< Comm >::HFaceType ElementType;
  };

  template< class Comm >
  struct ALUHElementType< 2, Comm >
  {
    typedef typename Dune::ALU3dBasicImplTraits< Comm >::HEdgeType ElementType;
  };

  template< class Comm >
  struct ALUHElementType< 3, Comm >
  {
    typedef typename Dune::ALU3dBasicImplTraits< Comm >::VertexType ElementType;
  };


  //*********************************************************
  //  LevelIterator Wrapper
  //*********************************************************
  template< class val_t >
  class IteratorWrapperInterface
    : public IteratorSTI< val_t >
  {
  public:
    virtual ~IteratorWrapperInterface () {}

    virtual int size  () = 0;
    virtual void next () = 0;
    virtual void first() = 0;
    virtual int done  () const = 0;
    virtual val_t & item () const = 0;
    virtual IteratorSTI< val_t > * clone () const { assert(false); abort(); return 0; }
  };

  typedef Dune::PartitionIteratorType PartitionIteratorType;

  // defines the pair of element and boundary
  template< int codim, class Comm >
  struct IteratorElType
  {
    typedef typename ALUHElementType< codim, Comm >::ElementType ElType;
    typedef typename Dune::ALU3dBasicImplTraits< Comm >::HBndSegType HBndSegType;
    typedef pair< ElType *, HBndSegType * > val_t;
  };

  template< int codim, PartitionIteratorType pitype, class Comm >
  class ALU3dGridLevelIteratorWrapper;

  // the element level iterator
  template< PartitionIteratorType pitype, class Comm >
  class ALU3dGridLevelIteratorWrapper< 0, pitype, Comm >
    : public IteratorWrapperInterface< typename IteratorElType< 0, Comm >::val_t >
  {
    typedef typename IteratorElType< 0, Comm >::ElType ElType;
    typedef typename IteratorElType< 0, Comm >::HBndSegType HBndSegType;
    typedef ALU3DSPACE LevelIterator< ElType > IteratorType;

    // the iterator
    IteratorType it_;

  public:
    typedef typename IteratorElType< 0, Comm >::val_t val_t;
    mutable val_t elem_;

    // constructor creating iterator
    template< class GridImp >
    ALU3dGridLevelIteratorWrapper ( const GridImp &grid, int level, const int nlinks )
      : it_( grid.myGrid(), level ),
        elem_( (ElType *) 0, (HBndSegType *) 0 )
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
  template< PartitionIteratorType pitype, class Comm >
  class ALU3dGridLevelIteratorWrapper< 1, pitype, Comm >
    : public IteratorWrapperInterface< typename IteratorElType< 1, Comm >::val_t >
  {
    typedef typename IteratorElType< 1, Comm >::ElType ElType;
    typedef typename IteratorElType< 1, Comm >::HBndSegType HBndSegType;
#ifdef ALUGRID_PERIODIC_BOUNDARY
    typedef ALU3DSPACE any_has_level_periodic< ElType > StopRule_t;
    typedef GridIterator< ElType, StopRule_t > IteratorType;
#else
    typedef ALU3DSPACE LevelIterator< ElType > IteratorType;
#endif

    // the iterator
    IteratorType it_;

  public:
    typedef typename IteratorElType< 1, Comm >::val_t val_t;
    mutable val_t elem_;

    // constructor creating iterator
    template< class GridImp >
    ALU3dGridLevelIteratorWrapper ( const GridImp &grid, int level, const int nlinks )
#ifdef ALUGRID_PERIODIC_BOUNDARY
      : it_( grid.myGrid(), StopRule_t(level) ),
#else
      : it_( grid.myGrid(), level ),
#endif
        elem_( (ElType *) 0, (HBndSegType*) 0 )
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
  template< PartitionIteratorType pitype, class Comm >
  class ALU3dGridLevelIteratorWrapper< 3, pitype, Comm >
    : public IteratorWrapperInterface< typename IteratorElType< 3, Comm >::val_t >
  {
    typedef typename IteratorElType< 3, Comm >::ElType ElType;
    typedef typename IteratorElType< 3, Comm >::HBndSegType HBndSegType;
    typedef Dune::ALU3dGridVertexList< Comm > VertexListType;
    typedef typename VertexListType::IteratorType IteratorType;

    VertexListType & vxList_;

    mutable int count_;
    const int size_;

  public:
    typedef typename IteratorElType< 3, Comm >::val_t val_t;
    mutable val_t elem_;

    // constructor creating iterator
    template< class GridImp >
    ALU3dGridLevelIteratorWrapper ( const GridImp &grid, int level, const int nlinks )
      : vxList_ ( grid.getVertexList( level ) ),
        count_( 0 ),
        size_( vxList_.size() ),
        elem_( (ElType *) 0, (HBndSegType *) 0 )
    {
      assert( vxList_.up2Date() );
    }

    // copy constructor
    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : vxList_(org.vxList_) , count_(org.count_) , size_(org.size_)
        , elem_(org.elem_)
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

  template< int codim, PartitionIteratorType pitype, class Comm >
  class ALU3dGridLeafIteratorWrapper;

  //typedef pair< ALUHElementType<0>::ElementType * , HBndSegType * > LeafValType;
  //typedef IteratorWrapperInterface<LeafValType> IteratorWrapperInterfaceType;

  //**********************************************************
  //  LeafIterator Wrapper
  //**********************************************************
  template< PartitionIteratorType pitype, class Comm >
  class ALU3dGridLeafIteratorWrapper< 0, pitype, Comm >
    : public IteratorWrapperInterface< typename IteratorElType< 0, Comm >::val_t >
  {
    typedef typename IteratorElType< 0, Comm >::ElType ElType;
    typedef typename IteratorElType< 0, Comm >::HBndSegType HBndSegType;
    typedef LeafIterator< ElType > IteratorType;

    // the ALU3dGrid Iterator
    IteratorType it_;

  public:
    typedef typename IteratorElType< 0, Comm >::val_t val_t;

  private:
    mutable val_t elem_;

  public:
    // constructor creating Iterator
    template< class GridImp >
    ALU3dGridLeafIteratorWrapper ( const GridImp &grid, int level, const int links )
      : it_( grid.myGrid() ),
        elem_( (ElType *) 0, (HBndSegType *) 0 )
    {}

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

  template< class ElType, PartitionIteratorType pitype, class Comm >
  struct LeafStopRule
  {
    typedef is_leaf_entity< ElType > StopRule_t;
  };

  // only in parallel we need only the interior items, in serial all items
  // are interior, to make the check fasterm this is only in parallel
  // implemented
#if ALU3DGRID_PARALLEL
  template< class ElType >
  struct LeafStopRule< ElType, Dune::Interior_Partition, MPI_Comm >
  {
    typedef is_interior_leaf_entity< ElType > StopRule_t;
  };
#endif // #if ALU3DGRID_PARALLEL


  template< PartitionIteratorType pitype, class Comm >
  class ALU3dGridLeafIteratorWrapper< 1, pitype, Comm >
    : public IteratorWrapperInterface< typename IteratorElType< 1, Comm >::val_t >
  {
    typedef typename IteratorElType< 1, Comm >::ElType ElType;
    typedef typename IteratorElType< 1, Comm >::HBndSegType HBndSegType;
    typedef typename LeafStopRule< ElType, pitype, Comm >::StopRule_t StopRule_t;
    typedef GridIterator< ElType, StopRule_t > IteratorType;

    // the face iterator
    IteratorType it_;

  public:
    typedef typename IteratorElType< 1, Comm >::val_t val_t;
  private:
    mutable val_t elem_;
  public:
    // constructor creating Iterator
    template< class GridImp >
    ALU3dGridLeafIteratorWrapper ( const GridImp &grid, int level, const int links )
      : it_( grid.myGrid(), StopRule_t() ),
        elem_( (ElType *) 0, (HBndSegType *) 0 )
    {}

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

  template< PartitionIteratorType pitype, class Comm >
  class ALU3dGridLeafIteratorWrapper< 2, pitype, Comm >
    : public IteratorWrapperInterface< typename IteratorElType< 2, Comm >::val_t >
  {
    typedef typename IteratorElType< 2, Comm >::ElType ElType;
    typedef typename IteratorElType< 2, Comm >::HBndSegType HBndSegType;
    typedef typename LeafStopRule< ElType, pitype, Comm >::StopRule_t StopRule_t;
    typedef GridIterator< ElType, StopRule_t > IteratorType;

  public:
    typedef typename IteratorElType< 2, Comm >::val_t val_t;

  private:
    // the edge iterator
    IteratorType it_;

    mutable val_t elem_;

  public:
    // constructor creating Iterator
    template< class GridImp >
    ALU3dGridLeafIteratorWrapper ( const GridImp &grid, int level, const int links )
      : it_( grid.myGrid(), StopRule_t() ),
        elem_( (ElType *) 0, (HBndSegType *) 0 )
    {}

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


  // the vertex leaf iterator, little bit different to the others
  template< PartitionIteratorType pitype, class Comm >
  class ALU3dGridLeafIteratorWrapper< 3, pitype, Comm >
    : public IteratorWrapperInterface< typename IteratorElType< 3, Comm >::val_t >
  {
    typedef typename IteratorElType< 3, Comm >::ElType ElType;
    typedef typename IteratorElType< 3, Comm >::HBndSegType HBndSegType;
    typedef Dune::ALU3dGridLeafVertexList< Comm > LeafVertexListType;
    typedef typename LeafVertexListType::IteratorType IteratorType;
    typedef typename LeafVertexListType::ItemType VxItemType;
    typedef typename LeafStopRule< ElType, pitype, Comm >::StopRule_t StopRule_t;

    LeafVertexListType & vxList_;
    typedef typename LeafVertexListType :: IteratorType ListIteratorType;

    mutable int count_;
    const int size_;

  public:
    typedef typename IteratorElType< 3, Comm >::val_t val_t;
    mutable val_t elem_;
    const StopRule_t rule_;

    // constructor creating iterator
    template< class GridImp >
    ALU3dGridLeafIteratorWrapper ( const GridImp &grid, int level, const int nlinks )
      : vxList_( grid.getLeafVertexList() ),
        count_( 0 ),
        size_( vxList_.size() ),
        elem_( (ElType *) 0, (HBndSegType *) 0 ),
        rule_()
    {
      assert( vxList_.up2Date() );
    }

    // copy constructor
    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper & org )
      : vxList_(org.vxList_)
        , count_(org.count_) , size_(org.size_)
        , elem_(org.elem_)
        , rule_()
    {}

    // returns size of leaf iterator, wrong here, return leaf size
    int size  ()  { return size_; }

    //! if level of item is larger then walk level, go next
    void next ()
    {
      ++count_;
      goNextValid();
      return ;
    }

    void first()
    {
      count_ = 0;
      goNextValid();
    }
    int done () const { return (count_ >= size_) ? 1 : 0; }
    val_t & item () const
    {
      assert( ! done () );
      assert( elem_.first );
      return elem_;
    }
  private:
    val_t & getItem () const
    {
      //elem_.first = vxList_.getItemList()[count_].first;
      assert( ! done () );
      elem_.first = vxList_.getItemList()[count_].first;
      return elem_;
    }
    void goNextValid()
    {
      if( done() ) return ;
      if( getItem().first == 0)
      {
        ++count_;
        goNextValid();
      }
      else
      {
        assert( elem_.first );
        if(! rule_( elem_.first ) )
        {
          ++count_;
          goNextValid();
        }
      }
    }
  };

  /*
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
     int done () const{ return it_->done(); }
     val_t & item () const
     {
      assert( ! done () );
      elem_.first  = & it_->item();
      return elem_;
     }
     };
   */

#if ALU3DGRID_PARALLEL
  template< int codim >
  class LeafLevelIteratorTTProxy
  {
    // type is hface_STI or hedge_STI
    typedef typename ALUHElementType< codim, MPI_Comm >::ElementType ElType;

    typedef typename Dune::ALU3dBasicImplTraits< MPI_Comm >::GitterImplType GitterImplType;

    typedef IteratorSTI< ElType > IteratorType;
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



  typedef pair< ALUHElementType< 0, MPI_Comm >::ElementType *, Dune::ALU3dBasicImplTraits< MPI_Comm >::HBndSegType * > LeafValType;

  //****************************
  //
  //  --GhostIterator
  //
  //****************************
  class ALU3dGridGhostIterator
    : public IteratorWrapperInterface< LeafValType >
  {
  public:
    typedef Dune::ALU3dBasicImplTraits< MPI_Comm >::GitterImplType GitterImplType;

    typedef Dune::ALU3dBasicImplTraits< MPI_Comm >::HElementType HElementType;
    typedef Dune::ALU3dBasicImplTraits< MPI_Comm >::HBndSegType HBndSegType;

  protected:
    GitterImplType & gitter_;

    // this tpye is hface_STI
    typedef ALUHElementType< 1, MPI_Comm >::ElementType ElType;

    typedef LeafLevelIteratorTTProxy< 1 > IteratorType;

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

    template< class GridImp >
    ALU3dGridGhostIterator ( const GridImp &grid, int level, const int nlinks )
      : gitter_( grid.myGrid() ),
        iterTT_( 0 ),
        it_( 0 ),
        nl_( nlinks ),
        link_( nlinks ), // makes default status == done
        elem_( (HElementType *) 0, (HBndSegType *) 0 )
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
          pair< HElementType *, HBndSegType * > elems( (HElementType *)0, (HBndSegType *)0 );
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
        pair< HElementType *, HBndSegType * > elems( (HElementType *)0, (HBndSegType *)0 );
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
  template<>
  class ALU3dGridLeafIteratorWrapper< 0, Dune::Ghost_Partition, MPI_Comm >
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
  template<>
  class ALU3dGridLevelIteratorWrapper< 0, Dune::Ghost_Partition, MPI_Comm >
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
    ALU3dGridLevelIteratorWrapper(const GridImp & grid,int level , const int nlinks )
      : ALU3dGridGhostIterator(grid,level,nlinks)
        , level_(level) , mxl_(grid.maxLevel()){}

    ALU3dGridLevelIteratorWrapper(const ALU3dGridLevelIteratorWrapper & org)
      : ALU3dGridGhostIterator(org) , level_(org.level_) , mxl_(org.mxl_){}
  };

  ///////////////////////////////////////////
  //
  //  Helper class to get item from Helement
  //
  //////////////////////////////////////////
  template< class GridImp, int cd >
  struct GetItem;

  template< class GridImp >
  struct GetItem< GridImp, 1 >
  {
    enum { cd = 1 };
    enum { elType = GridImp::elementType };

    typedef typename GridImp::MPICommunicatorType Comm;

    typedef typename Dune::ALU3dBasicImplTraits< Comm >::HElementType HElementType;
    typedef typename Dune::ALU3dImplTraits< GridImp::elementType, Comm >::GEOElementType GEOElementType;
    typedef typename IteratorElType< 1, Comm >::ElType ItemType;

    static ItemType *getItemFromEl ( typename Dune::ALU3dImplTraits< Dune::tetra, Comm >::GEOElementType &el, int i )
    {
      return el.myhface3( i );
    }

    static ItemType *getItemFromEl ( typename Dune::ALU3dImplTraits< Dune::hexa, Comm >::GEOElementType &el, int i )
    {
      return el.myhface4( i );
    }

    static ItemType *getItem ( HElementType &el, int i )
    {
      return getItemFromEl( static_cast< GEOElementType & >( el ), i );
    }

    static int numItems ()
    {
      return Dune::EntityCount< elType >::numFaces;
    }
  };

  template< class GridImp >
  struct GetItem< GridImp, 2 >
  {
    enum { cd = 2 };
    enum { elType = GridImp::elementType };

    typedef typename GridImp::MPICommunicatorType Comm;

    typedef typename Dune::ALU3dBasicImplTraits< Comm >::HElementType HElementType;
    typedef typename Dune::ALU3dImplTraits< GridImp::elementType, Comm >::GEOElementType GEOElementType;
    typedef typename IteratorElType< 2, Comm >::ElType ItemType;

    static ItemType *getItem ( HElementType &el, int i )
    {
      return static_cast< GEOElementType & >( el ).myhedge1( i );
    }

    static int numItems ()
    {
      return Dune::EntityCount<elType>::numEdges;
    }
  };

  template< class GridImp >
  struct GetItem< GridImp, 3 >
  {
    enum { cd = 3 };
    enum { elType = GridImp::elementType };

    typedef typename GridImp::MPICommunicatorType Comm;

    typedef typename Dune::ALU3dBasicImplTraits< Comm >::HElementType HElementType;
    typedef typename Dune::ALU3dImplTraits< GridImp::elementType, Comm >::GEOElementType GEOElementType;
    typedef typename IteratorElType< 3, Comm >::ElType ItemType;

    static ItemType *getItem ( HElementType &el, int i )
    {
      return static_cast< GEOElementType & >( el ).myvertex( i );
    }

    static int numItems ()
    {
      return Dune::EntityCount< elType >::numVertices;
    }
  };


  //! Ghost Iterator
  template< int codim >
  class ALU3dGridGhostIteratorHigherCodim
    : public IteratorWrapperInterface< typename IteratorElType< codim, MPI_Comm >::val_t >
  {
  public:
    typedef typename Dune::ALU3dBasicImplTraits< MPI_Comm >::HElementType HElementType;
    typedef typename Dune::ALU3dBasicImplTraits< MPI_Comm >::HBndSegType HBndSegType;
    typedef typename Dune::ALU3dBasicImplTraits< MPI_Comm >::GhostPairType GhostPairType;
    typedef typename IteratorElType< codim, MPI_Comm >::ElType ElType;
    typedef typename IteratorElType< codim, MPI_Comm >::val_t val_t;

  private:
    template< Dune::ALU3dGridElementType elType, int cd >
    struct SelectVector;

    template< Dune::ALU3dGridElementType elType >
    struct SelectVector< elType, 1 >
    {
      typedef typename Dune::ALU3dImplTraits< elType, MPI_Comm >::GEOElementType GEOElementType;

      static const vector< int > &getNotOnItemVector ( int face )
      {
        return GEOElementType::facesNotOnFace( face );
      }
    };

    template< Dune::ALU3dGridElementType elType >
    struct SelectVector< elType, 2 >
    {
      typedef typename Dune::ALU3dImplTraits< elType, MPI_Comm >::GEOElementType GEOElementType;
      static const vector< int > &getNotOnItemVector( int face )
      {
        return GEOElementType::edgesNotOnFace( face );
      }
    };

    template< Dune::ALU3dGridElementType elType >
    struct SelectVector< elType, 3 >
    {
      typedef typename Dune::ALU3dImplTraits< elType, MPI_Comm >::GEOElementType GEOElementType;
      static const vector< int > &getNotOnItemVector ( int face )
      {
        return GEOElementType::verticesNotOnFace( face );
      }
    };

    typedef ElType *getItemFunc_t ( HElementType &el, int i );

  private:
    typedef Dune :: ALU3dGridItemListType GhostItemListType;
    GhostItemListType &ghList_;
    typedef typename GhostItemListType :: IteratorType IteratorType;
    IteratorType curr_;
    IteratorType end_;
    mutable val_t elem_;
    mutable int count_;

  public:
    template< class GhostElementIteratorImp, class GridImp >
    ALU3dGridGhostIteratorHigherCodim ( GhostElementIteratorImp *, const GridImp &grid,
                                        int level, const int nlinks, GhostItemListType &ghList )
      : ghList_( ghList ),
        elem_( (ElType *) 0, (HBndSegType *) 0 ),
        count_( 0 )
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
  template<>
  class ALU3dGridLeafIteratorWrapper< 1, Dune::Ghost_Partition, MPI_Comm >
    : public ALU3dGridGhostIteratorHigherCodim< 1 >
  {
    enum { codim = 1 };
    typedef ALU3dGridLeafIteratorWrapper< 0, Dune::Ghost_Partition, MPI_Comm > GhostElementIteratorType;

  public:
    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLeafList(codim)) {}

    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the leaf ghost partition iterator
  template<>
  class ALU3dGridLeafIteratorWrapper< 2, Dune::Ghost_Partition, MPI_Comm >
    : public ALU3dGridGhostIteratorHigherCodim< 2 >
  {
    enum { codim = 2 };
    typedef ALU3dGridLeafIteratorWrapper< 0, Dune::Ghost_Partition, MPI_Comm > GhostElementIteratorType;

  public:
    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLeafList(codim)) {}

    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the leaf ghost partition iterator
  template<>
  class ALU3dGridLeafIteratorWrapper< 3, Dune::Ghost_Partition, MPI_Comm >
    : public ALU3dGridGhostIteratorHigherCodim< 3 >
  {
    enum { codim = 3 };
    typedef ALU3dGridLeafIteratorWrapper< 0, Dune::Ghost_Partition, MPI_Comm > GhostElementIteratorType;

  public:
    template <class GridImp>
    ALU3dGridLeafIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLeafList(codim)) {}

    ALU3dGridLeafIteratorWrapper (const ALU3dGridLeafIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the level ghost partition iterator
  template<>
  class ALU3dGridLevelIteratorWrapper< 1, Dune::Ghost_Partition, MPI_Comm >
    : public ALU3dGridGhostIteratorHigherCodim< 1 >
  {
    enum { codim = 1 };
    typedef ALU3dGridLevelIteratorWrapper< 0, Dune::Ghost_Partition, MPI_Comm > GhostElementIteratorType;

  public:
    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLevelList(codim,level)) {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the level ghost partition iterator
  template<>
  class ALU3dGridLevelIteratorWrapper< 2, Dune::Ghost_Partition, MPI_Comm >
    : public ALU3dGridGhostIteratorHigherCodim< 2 >
  {
    enum { codim = 2 };
    typedef ALU3dGridLevelIteratorWrapper< 0, Dune::Ghost_Partition, MPI_Comm > GhostElementIteratorType;

  public:
    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLevelList(codim,level)) {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the level ghost partition iterator
  template<>
  class ALU3dGridLevelIteratorWrapper< 3, Dune::Ghost_Partition, MPI_Comm >
    : public ALU3dGridGhostIteratorHigherCodim< 3 >
  {
    enum { codim = 3 };
    typedef ALU3dGridLevelIteratorWrapper< 0, Dune::Ghost_Partition, MPI_Comm > GhostElementIteratorType;

  public:
    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : ALU3dGridGhostIteratorHigherCodim<codim>((GhostElementIteratorType *) 0,grid,level,nlinks,grid.getGhostLevelList(codim,level)) {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : ALU3dGridGhostIteratorHigherCodim<codim>(org) {}
  };

  // the all partition iterator
  template<>
  class ALU3dGridLeafIteratorWrapper< 0, Dune::All_Partition, MPI_Comm >
    : public IteratorWrapperInterface< IteratorElType< 0, MPI_Comm >::val_t >
  {
    enum { codim = 0 };
    typedef ALU3dGridLeafIteratorWrapper< codim, Dune::InteriorBorder_Partition, MPI_Comm > InteriorIteratorType;
    typedef ALU3dGridLeafIteratorWrapper< codim, Dune::Ghost_Partition, MPI_Comm > GhostIteratorType;

  public:
    typedef IteratorElType< codim, MPI_Comm >::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType, val_t > IteratorType;
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
  template<>
  class ALU3dGridLeafIteratorWrapper< 1, Dune::All_Partition, MPI_Comm >
    : public IteratorWrapperInterface< IteratorElType< 1, MPI_Comm >::val_t >
  {
    enum { codim = 1 };
    typedef ALU3dGridLeafIteratorWrapper< codim, Dune::InteriorBorder_Partition, MPI_Comm > InteriorIteratorType;
    typedef ALU3dGridLeafIteratorWrapper< codim, Dune::Ghost_Partition, MPI_Comm > GhostIteratorType;

  public:
    typedef IteratorElType< codim, MPI_Comm >::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType, val_t > IteratorType;
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
  template<>
  class ALU3dGridLeafIteratorWrapper< 2, Dune::All_Partition, MPI_Comm >
    : public IteratorWrapperInterface< IteratorElType< 2, MPI_Comm >::val_t >
  {
    enum { codim = 2 };
    typedef ALU3dGridLeafIteratorWrapper< codim, Dune::InteriorBorder_Partition, MPI_Comm > InteriorIteratorType;
    typedef ALU3dGridLeafIteratorWrapper< codim, Dune::Ghost_Partition, MPI_Comm > GhostIteratorType;

  public:
    typedef IteratorElType< codim, MPI_Comm >::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType, val_t > IteratorType;
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
  template<>
  class ALU3dGridLeafIteratorWrapper< 3, Dune::All_Partition, MPI_Comm >
    : public IteratorWrapperInterface< IteratorElType< 3, MPI_Comm >::val_t >
  {
    enum { codim = 3 };
    typedef ALU3dGridLeafIteratorWrapper< codim, Dune::InteriorBorder_Partition, MPI_Comm > InteriorIteratorType;
    typedef ALU3dGridLeafIteratorWrapper< codim, Dune::Ghost_Partition, MPI_Comm > GhostIteratorType;

  public:
    typedef IteratorElType< codim, MPI_Comm >::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType, val_t > IteratorType;
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
  template<>
  class ALU3dGridLevelIteratorWrapper< 0, Dune::All_Partition, MPI_Comm >
    : public IteratorWrapperInterface< LeafValType >
  {
    typedef ALU3dGridLevelIteratorWrapper< 0, Dune::InteriorBorder_Partition, MPI_Comm > InteriorIteratorType;
    typedef ALU3dGridLevelIteratorWrapper< 0, Dune::Ghost_Partition, MPI_Comm > GhostIteratorType;

  public:
    typedef LeafValType val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType, val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, level, nlinks ) ,
                GhostIteratorType    ( grid, level, nlinks ) )
    {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org)
      : iter_(org.iter_) {}

    int size  () { return iter_.size(); }
    void next () { iter_.next(); }
    void first() { iter_.first(); }
    int done () const {return iter_.done(); }
    val_t & item () const { assert( ! done() ); return iter_.item(); }
  };
#endif // #if ALU3DGRID_PARALLEL

  // placed here because we need ALU3dGridLevelIteratorWrapper<0,Dune::All_Partition> here
  // the edge level iterator
  template< PartitionIteratorType pitype, class Comm >
  class ALU3dGridLevelIteratorWrapper< 2, pitype, Comm >
    : public IteratorWrapperInterface< typename IteratorElType< 2, Comm >::val_t >
  {
  public:
    typedef typename ALUHElementType< 2, Comm >::ElementType ElType;
    typedef typename Dune::ALU3dBasicImplTraits< Comm >::HBndSegType HBndSegType;
    typedef typename Dune::ALU3dBasicImplTraits< Comm >::GEOEdgeType GEOEdgeType;

    typedef typename IteratorElType< 2, Comm >::val_t val_t;

  private:
    mutable val_t elem_;
    const int level_;

    typedef Dune :: ALU3dGridItemListType ItemListType;
    ItemListType & edgeList_;

    size_t count_ ;
    bool maxLevel_;

  public:
    // constructor creating iterator
    template< class GridImp >
    ALU3dGridLevelIteratorWrapper ( const GridImp &grid, int level, const int nlinks )
      : elem_( (ElType *) 0, (HBndSegType *) 0 ),
        level_( level ),
        edgeList_( grid.getEdgeList( level ) ),
        count_( 0 )
    {
      if( ! edgeList_.up2Date() )
        updateEdgeList(grid,level,nlinks);
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
    void updateEdgeList(const GridImp & grid, int level, int nlinks)
    {
      typedef ALU3dGridLevelIteratorWrapper< 0, Dune::All_Partition, Comm > ElementLevelIterator;
      typedef typename ElementLevelIterator :: val_t el_val_t;
      ElementLevelIterator iter(grid,level,nlinks);

      edgeList_.getItemList().resize(0);
      map < int , int > visited;

      for( iter.first(); ! iter.done(); iter.next() )
      {
        typedef typename Dune::ALU3dImplTraits< GridImp::elementType, Comm >::GEOElementType GEOElementType;
        enum { numEdges = Dune::EntityCount< GridImp::elementType >::numEdges };

        GEOElementType *elem = 0;
        el_val_t & item = iter.item();

        if( item.first )
          elem = static_cast< GEOElementType * > (item.first);
        else if( item.second )
          elem = static_cast< GEOElementType * > (item.second->getGhost().first);

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
  template<>
  class ALU3dGridLevelIteratorWrapper< 1, Dune::All_Partition, MPI_Comm >
    : public IteratorWrapperInterface< IteratorElType< 1, MPI_Comm >::val_t >
  {
    enum { codim = 1 };
    typedef ALU3dGridLevelIteratorWrapper< codim, Dune::InteriorBorder_Partition, MPI_Comm > InteriorIteratorType;
    typedef ALU3dGridLevelIteratorWrapper< codim, Dune::Ghost_Partition, MPI_Comm > GhostIteratorType;

  public:
    typedef IteratorElType< codim, MPI_Comm >::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType, val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, level, nlinks ) ,
                GhostIteratorType    ( grid, level, nlinks ) )
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
  template<>
  class ALU3dGridLevelIteratorWrapper< 2, Dune::All_Partition, MPI_Comm >
    : public IteratorWrapperInterface< IteratorElType< 2, MPI_Comm >::val_t >
  {
    enum { codim = 2 };
    typedef ALU3dGridLevelIteratorWrapper< codim, Dune::InteriorBorder_Partition, MPI_Comm > InteriorIteratorType;
    typedef ALU3dGridLevelIteratorWrapper< codim, Dune::Ghost_Partition, MPI_Comm > GhostIteratorType;

  public:
    typedef IteratorElType< codim, MPI_Comm >::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType, val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, level, nlinks ) ,
                GhostIteratorType    ( grid, level, nlinks ) )
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
  template<>
  class ALU3dGridLevelIteratorWrapper< 3, Dune::All_Partition, MPI_Comm >
    : public IteratorWrapperInterface < IteratorElType< 3, MPI_Comm >::val_t >
  {
    enum { codim = 3 };
    typedef ALU3dGridLevelIteratorWrapper< codim, Dune::InteriorBorder_Partition, MPI_Comm > InteriorIteratorType;
    typedef ALU3dGridLevelIteratorWrapper< codim, Dune::Ghost_Partition, MPI_Comm > GhostIteratorType;

  public:
    typedef IteratorElType< codim, MPI_Comm >::val_t val_t;
    // use ALUGrids AlignIterator to combine Interior and Ghost Iterator
    typedef AlignIterator< InteriorIteratorType, GhostIteratorType, val_t > IteratorType;
  private:
    IteratorType iter_;
  public:

    template <class GridImp>
    ALU3dGridLevelIteratorWrapper (const GridImp & grid, int level , const int nlinks )
      : iter_ ( InteriorIteratorType ( grid, level, nlinks ) ,
                GhostIteratorType    ( grid, level, nlinks ) )
    {}

    ALU3dGridLevelIteratorWrapper (const ALU3dGridLevelIteratorWrapper & org )
      : iter_ (org.iter_) {}

    int size  () { return iter_.size(); }
    void next () { iter_.next(); }
    void first() { iter_.first(); }
    int done () const {return iter_.done(); }
    val_t & item () const { assert( ! done() ); return iter_.item(); }
  };
#endif // #if ALU3DGRID_PARALLEL

} // end namespace ALUGridSpace

#endif // #ifndef DUNE_ALU3DITERATORS_HH
