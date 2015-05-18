// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/common/exceptions.hh>

#include "geometry.hh"
#include "grid.hh"

namespace Dune {

  template<int cd, int dim, class GridImp>
  inline void ALU3dGridEntity<cd,dim,GridImp> ::
  reset( int l )
  {
    item_  = 0;
    level_ = l;
    twist_ = 0;
    face_  = -1;
  }

  template<int cd, int dim, class GridImp>
  inline void ALU3dGridEntity<cd,dim,GridImp> ::
  removeElement()
  {
    item_ = 0;
    geo_.invalidate();
  }

  template<int cd, int dim, class GridImp>
  inline bool ALU3dGridEntity<cd,dim,GridImp> ::
  equals(const ALU3dGridEntity<cd,dim,GridImp> & org) const
  {
    return (item_ == org.item_);
  }

  template<int cd, int dim, class GridImp>
  inline int ALU3dGridEntity<cd,dim,GridImp> :: getIndex () const
  {
    return gIndex_;
  }

  template<int cd, int dim, class GridImp>
  inline int ALU3dGridEntity<cd,dim,GridImp> :: level () const
  {
    return level_;
  }

  template<int cd, int dim, class GridImp>
  inline PartitionType ALU3dGridEntity<cd,dim,GridImp> ::
  partitionType () const
  {
    return partitionType_;
  }

  template<int cd, int dim, class GridImp>
  inline GeometryType
  ALU3dGridEntity<cd,dim,GridImp>:: type () const
  {
    return geo_.type();
  }

  /////////////////////////////////////////////////////////////////
  //
  //  --Entity0
  //  --Codim0Entity
  //
  ////////////////////////////////////////////////////////////////
  template<int dim, class GridImp>
  inline void ALU3dGridEntity<0,dim,GridImp> ::
  removeElement ()
  {
    item_  = 0;
    ghost_ = 0;
    geo_.invalidate();
  }

  template<int dim, class GridImp>
  inline void ALU3dGridEntity<0,dim,GridImp> ::
  reset (int walkLevel )
  {
    item_       = 0;
    ghost_      = 0;
    level_      = -1;
    isLeaf_     = false;

    // reset geometry information
    geo_.invalidate();
  }

  // works like assignment
  template<int dim, class GridImp>
  inline void
  ALU3dGridEntity<0,dim,GridImp> :: setEntity(const ALU3dGridEntity<0,dim,GridImp> & org)
  {
    item_          = org.item_;
    ghost_         = org.ghost_;
    level_         = org.level_;
    isLeaf_        = org.isLeaf_;

    // reset geometry information
    geo_.invalidate();
  }

  template<int dim, class GridImp>
  inline void
  ALU3dGridEntity<0,dim,GridImp>::
  setElement(const EntitySeed& key )
  {
    if( ! key.isGhost() )
      setElement( *key.interior() );
    else
      setGhost( *key.ghost() );
  }

  template<int dim, class GridImp>
  inline void
  ALU3dGridEntity<0,dim,GridImp>::
  setElement(HElementType & element)
  {
    item_ = static_cast<IMPLElementType *> (&element);
    assert( item_ );
    // make sure this method is not called for ghosts
    assert( ! item_->isGhost() );
    ghost_   = 0;
    level_   = (*item_).level();
    isLeaf_  = ((*item_).down() == 0);

    // reset geometry information
    geo_.invalidate();
  }

  template<int dim, class GridImp>
  inline void
  ALU3dGridEntity<0,dim,GridImp> :: setGhost(HBndSegType & ghost)
  {
    // use element as ghost
    item_  = static_cast<IMPLElementType *> ( ghost.getGhost().first );

    // method getGhost can return 0, but then is something wrong
    assert(item_);
    assert(item_->isGhost());

    level_   = item_->level();
    // remember pointer to ghost face
    ghost_ = static_cast<BNDFaceType *> (&ghost);
    assert( ghost_ );

    BNDFaceType * dwn = static_cast<BNDFaceType *> (ghost.down());
    if ( ! dwn ) isLeaf_ = true;
    else
    {
      assert( ghost.level() == level_ );
      if(dwn->ghostLevel() == level_)
        isLeaf_ = true;
      else
        isLeaf_ = false;
    }
    // check wether ghost is leaf or not, ghost leaf means
    // that this is the ghost that we want in the leaf iterator
    // not necessarily is real leaf element
    // see Intersection Iterator, same story

    // reset geometry information
    geo_.invalidate();
  }

  template<int dim, class GridImp>
  inline int
  ALU3dGridEntity<0,dim,GridImp> :: level() const
  {
    return level_;
  }

  template<int dim, class GridImp>
  inline bool ALU3dGridEntity<0,dim,GridImp> ::
  equals (const ALU3dGridEntity<0,dim,GridImp> &org ) const
  {
    return (item_ == org.item_);
  }

  template<int dim, class GridImp>
  inline GeometryType
  ALU3dGridEntity<0,dim,GridImp> :: type () const
  {
    return geo_.type();
  }

  template<int dim, class GridImp>
  inline int ALU3dGridEntity<0,dim,GridImp> :: getIndex() const
  {
    assert( item_ );
    return (*item_).getIndex();
  }

  template<int dim, class GridImp>
  template<int cc>
  inline int ALU3dGridEntity<0,dim,GridImp> :: count () const
  {
    return grid().referenceElement().size(cc);
  }

  template<int dim, class GridImp>
  inline unsigned int ALU3dGridEntity<0,dim,GridImp> :: subEntities (unsigned int codim) const
  {
    return grid().referenceElement().size(codim);
  }

  template<int dim, class GridImp>
  inline PartitionType ALU3dGridEntity<0,dim,GridImp> ::
  partitionType () const
  {
    assert( item_ );
    // make sure we really got a ghost
    assert( (isGhost()) ? item_->isGhost() : true );
    return (isGhost() ?  GhostEntity : InteriorEntity);
  }

  template<int dim, class GridImp>
  inline bool ALU3dGridEntity<0,dim,GridImp> :: isLeaf() const
  {
    return isLeaf_;
  }

  template<int dim, class GridImp>
  inline ALU3dGridHierarchicIterator<GridImp>
  ALU3dGridEntity<0,dim,GridImp> :: hbegin (int maxlevel) const
  {
    assert(item_ != 0);
    // if isGhost is true the end iterator will be returned
    if( isGhost() )
    {
      return ALU3dGridHierarchicIterator<GridImp>(factory_,*ghost_,maxlevel, isLeaf() );
    }
    return ALU3dGridHierarchicIterator<GridImp>(factory_,*item_,maxlevel, isLeaf() );
  }

  template<int dim, class GridImp>
  inline ALU3dGridHierarchicIterator<GridImp> ALU3dGridEntity<0,dim,GridImp> :: hend (int maxlevel) const
  {
    assert(item_ != 0);
    return ALU3dGridHierarchicIterator<GridImp> (factory_, *item_, maxlevel, true);
  }

  template<int dim, class GridImp>
  inline typename ALU3dGridEntity<0,dim,GridImp> :: ALU3dGridLeafIntersectionIteratorType
  ALU3dGridEntity<0,dim,GridImp> :: ileafbegin () const
  {
    assert(item_ != 0);
    return ALU3dGridIntersectionIteratorType (*this, this->level(), false);
  }

  template<int dim, class GridImp>
  inline typename ALU3dGridEntity<0,dim,GridImp> :: ALU3dGridLeafIntersectionIteratorType
  ALU3dGridEntity<0,dim,GridImp> :: ileafend () const
  {
    assert(item_ != 0);
    return ALU3dGridLeafIntersectionIteratorType (*this, this->level(), true);
  }

  template<int dim, class GridImp>
  inline typename ALU3dGridEntity<0,dim,GridImp> :: ALU3dGridLevelIntersectionIteratorType
  ALU3dGridEntity<0,dim,GridImp> :: ilevelbegin () const
  {
    assert(item_ != 0);
    // disable level intersection iterator for conforming refinement
    return ALU3dGridLevelIntersectionIteratorType (*this, this->level(), grid().conformingRefinement() );
  }

  template<int dim, class GridImp>
  inline typename ALU3dGridEntity<0,dim,GridImp> :: ALU3dGridLevelIntersectionIteratorType
  ALU3dGridEntity<0,dim,GridImp> :: ilevelend () const
  {
    assert(item_ != 0);
    return ALU3dGridLevelIntersectionIteratorType (*this, this->level(), true);
  }

  // Adaptation methods
  template<int dim, class GridImp>
  inline bool ALU3dGridEntity<0,dim,GridImp> :: isNew () const
  {
    assert( item_ );
    return item_->hasBeenRefined();
  }

  template<int dim, class GridImp>
  inline bool ALU3dGridEntity<0,dim,GridImp> :: mightVanish () const
  {
    assert( item_ );
    return ((*item_).requestrule() == coarse_element_t);
  }

  //*******************************************************************
  //
  //  --EntityPointer
  //  --EnPointer
  //
  //*******************************************************************
  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const FactoryType& factory,
                             const HElementType &item)
    : factory_(factory)
      , seed_( item )
      , entity_( 0 )
  {}

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const FactoryType& factory,
                             const HBndSegType & ghostFace )
    : factory_(factory)
      , seed_( ghostFace )
      , entity_ ( factory_.template getNewEntity<codim> ( ghostFace.level() ))
  {
    // sets entity and item pointer
    updateGhostPointer( const_cast<HBndSegType &> (ghostFace) );
  }

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const FactoryType& factory,
                             const ALU3dGridEntitySeedType& key )
    : factory_(factory)
      , seed_( key )
      , entity_ ( 0 )
  {}

  // constructor Level,Leaf and HierarchicIterator
  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const FactoryType& factory, int level )
    : factory_(factory)
      , seed_()
      , entity_ ( factory_.template getNewEntity<codim> ( level ) )
  {
    // this needs to be called
    // have to investigate why
    entityImp().reset(level);
  }


  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const FactoryType& factory,
                             const HElementType &item,
                             const int level,
                             const int twist,
                             const int duneFace )
    : factory_(factory)
      , seed_( item, level, twist, duneFace )
      , entity_( 0 )
  {}

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const ALU3dGridEntityPointerType & org)
    : factory_(org.factory_)
      , seed_( org.seed_ )
      , entity_( 0 )
  {
    // if entity exists then copy entity
    getEntity( org );
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointerBase<codim,GridImp> ::
  getEntity(const ALU3dGridEntityPointerType & org)
  {
    // if entity existed for original pointer then copy
    if( org.entity_ )
    {
      assert( entity_ == 0 );
      entity_ = factory_.template getNewEntity<codim> ();
      // set entity right away
      entityImp().setEntity( org.entityImp() );
    }
  }

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> &
  ALU3dGridEntityPointerBase<codim,GridImp> ::
  operator = (const ALU3dGridEntityPointerType & org)
  {
    clone( org );
    return *this;
  }

  template<int codim, class GridImp >
  inline void
  ALU3dGridEntityPointerBase<codim,GridImp> ::
  clone (const ALU3dGridEntityPointerType & org)
  {
    assert( &factory_ == &org.factory_ );

    // set item
    seed_ = org.seed_;

    HElementType* item = seed_.item();

    if( item )
    {
      // if no entity check org entity
      // if no org entity then nothing is done
      if( !entity_ )
      {
        getEntity(org);
      }
      else
      {
        // in case of ghost element use different set method
        if( item->isGhost() )
        {
          // on ghosts entity pointers entity always exists
          assert( org.entity_ );
          entityImp().setEntity( org.entityImp() );
        }
        else
        {
          // otherwise item is set
          entityImp().setElement( seed_ );
        }
      }
    }
    else
    {
      this->done();
    }
    return ;
  }

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ~ALU3dGridEntityPointerBase()
  {
    this->done();
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointerBase<codim,GridImp>::done ()
  {
    seed_.clear();
    // free entity
    freeEntity();
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointerBase<codim,GridImp>::freeEntity ()
  {
    // sets entity pointer in the status of an empty entity
    if( entity_ )
    {
      entityImp().removeElement();
      factory_.template freeEntity<codim> ( (EntityObject *) entity_ );
      entity_ = 0;
    }
  }

  template<int codim, class GridImp >
  inline bool ALU3dGridEntityPointerBase<codim,GridImp>::
  equals (const ALU3dGridEntityPointerBase<codim,GridImp>& i) const
  {
    // check equality of underlying items
    return (seed_.equals( i.seed_ ));
  }

  template<int codim, class GridImp >
  inline typename ALU3dGridEntityPointerBase<codim,GridImp>::Entity &
  ALU3dGridEntityPointerBase<codim,GridImp>::dereference () const
  {
    // don't dereference empty entity pointer
    assert( seed_.item() );
    if( ! entity_ )
    {
      entity_ = factory_.template getNewEntity<codim> ();
      entityImp().setElement( seed_ );
    }
    assert( seed_.item() == & entityImp().getItem() );
    return (*entity_);
  }

  template<int codim, class GridImp >
  inline int ALU3dGridEntityPointerBase<codim,GridImp>::level () const
  {
    assert( seed_.item() );
    return seed_.item()->level();
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointerBase<codim,GridImp>::
  updateGhostPointer( HBndSegType & ghostFace )
  {
    assert( entity_ );
    entityImp().setGhost( ghostFace );
    // inside the method setGhost the method getGhost of the ghostFace is
    // called and set as item
    seed_.set( ghostFace );
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointerBase<codim,GridImp>::
  updateEntityPointer( HElementType * item , int )
  {
    seed_.set( *item );
    if( item && entity_ )
    {
      entityImp().setElement( seed_ );
    }
  }

  ///////////////////////////////////////////////////////////////////
  //
  //  specialisation for higher codims
  //
  ///////////////////////////////////////////////////////////////////

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointer<codim,GridImp> ::
  ALU3dGridEntityPointer(const FactoryType& factory,
                         const int level,
                         const HElementType &item,
                         const int twist,
                         const int duneFace )
    : ALU3dGridEntityPointerBase<codim,GridImp> (factory,item,level,twist,duneFace)
  {}

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointer<codim,GridImp> ::
  ALU3dGridEntityPointer(const ALU3dGridEntityPointerType & org)
    : ALU3dGridEntityPointerBase<codim,GridImp>(org)
  {}

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointer<codim,GridImp> &
  ALU3dGridEntityPointer<codim,GridImp>::
  operator = (const ALU3dGridEntityPointerType & org)
  {
    // clone pointer
    clone(org);
    return *this;
  }

  template<int codim, class GridImp >
  inline void
  ALU3dGridEntityPointer<codim,GridImp> ::
  clone (const ALU3dGridEntityPointerType & org)
  {
    // copy key
    seed_ = org.seed_;

    assert( &factory_ == &org.factory_ );

    // if entity exists, just remove item pointer
    if( seed_.item() )
    {
      if( ! entity_ )
        getEntity(org);
      else
        entityImp().setElement( seed_ );
    }
    else
      this->done();
    return ;
  }

  template<int codim, class GridImp >
  inline typename ALU3dGridEntityPointer<codim,GridImp>::Entity &
  ALU3dGridEntityPointer<codim,GridImp>::dereference () const
  {
    // don't dereference empty entity pointer
    assert( seed_.item() );
    if( ! entity_ )
    {
      entity_ = factory_.template getNewEntity<codim> ();
      entityImp().setElement( seed_ );
    }
    assert( seed_.item() == & entityImp().getItem() );
    return (*entity_);
  }

  template<int codim, class GridImp >
  inline int ALU3dGridEntityPointer<codim,GridImp>::level () const
  {
    return seed_.level();
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointer<codim,GridImp>::
  updateEntityPointer( HElementType * item, int level)
  {
    seed_.set( *item, level );
    if( item && entity_ )
    {
      entityImp().setElement( seed_ );
    }
  }


} // end namespace Dune
