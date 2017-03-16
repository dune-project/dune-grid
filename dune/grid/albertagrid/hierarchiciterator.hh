// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_HIERARCHICITERATOR_HH
#define DUNE_ALBERTA_HIERARCHICITERATOR_HH

#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/common/entityiterator.hh>

#if HAVE_ALBERTA

namespace Dune
{

  // AlbertaGridHierarchicIterator
  // -----------------------------

  /*!
     Mesh entities of codimension 0 ("elements") allow to visit all entities of
     codimension 0 obtained through nested, hierarchic refinement of the entity.
     Iteration over this set of entities is provided by the HIerarchicIterator,
     starting from a given entity.
     This is redundant but important for memory efficient implementations of unstru
     hierarchically refined meshes.
   */
  template< class GridImp >
  class AlbertaGridHierarchicIterator
  {
    typedef AlbertaGridHierarchicIterator< GridImp > This;

  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::ctype ctype;

    typedef MakeableInterfaceObject< Entity > EntityObject;
    typedef typename EntityObject::ImplementationType EntityImp;

    typedef typename EntityImp::ElementInfo ElementInfo;

    AlbertaGridHierarchicIterator ()
    {}

    //! the normal Constructor
    AlbertaGridHierarchicIterator ( const GridImp &grid,
                                    const ElementInfo &elementInfo,
                                    int maxLevel );

    //! the default Constructor
    AlbertaGridHierarchicIterator ( const GridImp &grid, int actLevel, int maxLevel );

    //! copy onstructor
    AlbertaGridHierarchicIterator ( const This &other );

    //! assignment operator
    This &operator= ( const This &other );

    //! increment
    void increment();

    //! equality
    bool equals ( const This &other ) const
    {
      return entityImp().equals( other.entityImp() );
    }

    //! dereferencing
    Entity &dereference () const
    {
      return entity_;
    }

    //! ask for level of entities
    int level () const
    {
      return entityImp().level();
    }

  protected:
    //! obtain reference to internal entity implementation
    EntityImp &entityImp ()
    {
      return GridImp::getRealImplementation( entity_ );
    }

    //! obtain const reference to internal entity implementation
    const EntityImp &entityImp () const
    {
      return GridImp::getRealImplementation( entity_ );
    }

    //! obtain a reference to the grid
    const GridImp &grid () const
    {
      return entityImp().grid();
    }

  private:
    void increment ( ElementInfo elementInfo );

    mutable  Entity entity_;

    // level on which the iterator was started
    int startLevel_;

    // maximal level to go down to
    int maxlevel_;
  };


  template< class GridImp >
  inline AlbertaGridHierarchicIterator< GridImp >
  ::AlbertaGridHierarchicIterator( const GridImp &grid, int actLevel, int maxLevel )
    : entity_( EntityImp( grid ) ),
      startLevel_( actLevel ),
      maxlevel_( maxLevel )
  {}


  template< class GridImp >
  inline AlbertaGridHierarchicIterator< GridImp >
  ::AlbertaGridHierarchicIterator ( const GridImp &grid,
                                    const ElementInfo &elementInfo,
                                    int maxLevel )
    : entity_( EntityImp( grid ) ),
      startLevel_( elementInfo.level() ),
      maxlevel_( maxLevel )
  {
    increment( elementInfo );
  }


  template< class GridImp >
  inline AlbertaGridHierarchicIterator< GridImp >
  ::AlbertaGridHierarchicIterator( const This &other )
    : entity_( other.entity_ ),
      startLevel_( other.startLevel_ ),
      maxlevel_( other.maxlevel_ )
  {}


  template< class GridImp >
  inline typename AlbertaGridHierarchicIterator< GridImp >::This &
  AlbertaGridHierarchicIterator< GridImp >::operator= ( const This &other )
  {
    entity_ = other.entity_;
    startLevel_ = other.startLevel_;
    maxlevel_ = other.maxlevel_;
    return *this;
  }


  template< class GridImp >
  inline void AlbertaGridHierarchicIterator< GridImp >::increment ()
  {
    increment( entityImp().elementInfo() );
  }

  template< class GridImp >
  inline void AlbertaGridHierarchicIterator< GridImp >
  ::increment ( ElementInfo elementInfo )
  {
    assert( !elementInfo == false );
    if( (elementInfo.level() >= maxlevel_) || elementInfo.isLeaf() )
    {
      while( (elementInfo.level() > startLevel_) && (elementInfo.indexInFather() == 1) )
        elementInfo = elementInfo.father();
      if( elementInfo.level() > startLevel_ )
        entityImp().setElement( elementInfo.father().child( 1 ), 0 );
      else
        entityImp().clearElement();
    }
    else
      entityImp().setElement( elementInfo.child( 0 ), 0 );
  }

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_HIERARCHICITERATOR_HH
