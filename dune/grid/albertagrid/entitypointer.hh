// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ENTITYPOINTER_HH
#define DUNE_ALBERTA_ENTITYPOINTER_HH

#include <dune/grid/albertagrid/elementinfo.hh>

#if HAVE_ALBERTA

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int dim, int dimworld >
  class AlbertaGrid;



  /** \class   AlbertaGridEntityPointer
   *  \ingroup AlbertaGrid
   *  \brief   EntityPointer implementation for AlbertaGrid
   */
  template< int codim, class GridImp >
  class AlbertaGridEntityPointer
  {
    typedef AlbertaGridEntityPointer< codim, GridImp > This;

    friend class AlbertaGrid< GridImp::dimension, GridImp::dimensionworld >;

  public:
    static const int dimension = GridImp::dimension;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;
    static const int dimensionworld = GridImp::dimensionworld;

    typedef typename GridImp::template Codim< codimension >::Entity Entity;

  protected:
    typedef MakeableInterfaceObject< Entity > EntityObject;
    typedef typename EntityObject::ImplementationType EntityImp;

  public:
    typedef AlbertaGridEntityPointer< codimension, GridImp > EntityPointerImp;

    typedef typename EntityImp::ElementInfo ElementInfo;

    AlbertaGridEntityPointer ();

    //! make an EntityPointer that points to an element
    AlbertaGridEntityPointer ( const GridImp &grid,
                               const ElementInfo &elementInfo,
                               int subEntity );

    //! constructor for invalid EntityPointer
    AlbertaGridEntityPointer ( const GridImp &grid );

    //! make entity pointer from entity
    AlbertaGridEntityPointer ( const EntityImp &entity );

#if 0
    //! Destructor
    ~AlbertaGridEntityPointer();
#endif

    //! equality
    bool equals ( const This &other ) const;

    //! dereferencing
    Entity &dereference () const;

    //! ask for level of entities
    int level () const;

  protected:
    //! obtain reference to internal entity implementation
    EntityImp &entityImp ();

    //! obtain const reference to internal entity implementation
    const EntityImp &entityImp () const;

    //! obtain a reference to the grid
    const GridImp &grid () const;

  private:
    mutable Entity entity_;
  };



  template< int codim, class GridImp >
  inline AlbertaGridEntityPointer< codim, GridImp >
  ::AlbertaGridEntityPointer ()
  {}


  template< int codim, class GridImp >
  inline AlbertaGridEntityPointer< codim, GridImp >
  ::AlbertaGridEntityPointer ( const GridImp &grid,
                               const ElementInfo &elementInfo,
                               int subEntity )
    : entity_( EntityImp( grid, elementInfo, subEntity ) )
  {}


  template<int codim, class GridImp >
  inline AlbertaGridEntityPointer< codim, GridImp >
  ::AlbertaGridEntityPointer ( const GridImp &grid )
    : entity_( EntityImp( grid ) )
  {}


  template< int codim, class GridImp >
  inline AlbertaGridEntityPointer< codim, GridImp >
  ::AlbertaGridEntityPointer ( const EntityImp &entity )
    : entity_( entity )
  {}


#if 0
  template<int codim, class GridImp >
  inline AlbertaGridEntityPointer< codim, GridImp >::~AlbertaGridEntityPointer ()
  {}
#endif


  template<int codim, class GridImp >
  inline bool
  AlbertaGridEntityPointer< codim, GridImp >::equals ( const This &other ) const
  {
    return entityImp().equals( other.entityImp() );
  }


  template<int codim, class GridImp >
  inline typename AlbertaGridEntityPointer< codim, GridImp >::Entity &
  AlbertaGridEntityPointer< codim, GridImp >::dereference () const
  {
    return entity_;
  }


  template< int codim, class GridImp >
  inline int AlbertaGridEntityPointer< codim, GridImp >::level () const
  {
    return entityImp().level();
  }


  template< int codim, class GridImp >
  inline typename AlbertaGridEntityPointer< codim, GridImp >::EntityImp &
  AlbertaGridEntityPointer< codim, GridImp >::entityImp ()
  {
    return GridImp::getRealImplementation( entity_ );
  }


  template< int codim, class GridImp >
  inline const typename AlbertaGridEntityPointer< codim, GridImp >::EntityImp &
  AlbertaGridEntityPointer< codim, GridImp >::entityImp () const
  {
    return GridImp::getRealImplementation( entity_ );
  }


  template< int codim, class GridImp >
  inline const GridImp &AlbertaGridEntityPointer< codim, GridImp >::grid () const
  {
    return entityImp().grid();
  }

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_ENTITYPOINTER_HH
