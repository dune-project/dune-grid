// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
      return entity_.impl().equals( other.entity_.impl() );
    }

    //! dereferencing
    Entity &dereference () const
    {
      return entity_;
    }

    //! ask for level of entities
    int level () const
    {
      return entity_.impl().level();
    }

  protected:
    //! obtain a reference to the grid
    const GridImp &grid () const
    {
      return entity_.impl().grid();
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
    increment( entity_.impl().elementInfo() );
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
        entity_.impl().setElement( elementInfo.father().child( 1 ), 0 );
      else
        entity_.impl().clearElement();
    }
    else
      entity_.impl().setElement( elementInfo.child( 0 ), 0 );
  }

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_HIERARCHICITERATOR_HH
