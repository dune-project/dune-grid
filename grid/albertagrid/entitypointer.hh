// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ENTITYPOINTER_HH
#define DUNE_ALBERTA_ENTITYPOINTER_HH

#include <dune/grid/common/entitypointer.hh>

#include <dune/grid/albertagrid/elementinfo.hh>

namespace Dune
{

  /*!
     Enables iteration over all entities of a given codimension and level of a grid.
   */
  template< int codim, class GridImp >
  class AlbertaGridEntityPointer
  {
    typedef AlbertaGridEntityPointer< codim, GridImp > This;

    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    friend class AlbertaGridEntity< codim, dim, GridImp >;
    friend class AlbertaGridEntity< 0, dim, GridImp >;
    friend class AlbertaGrid< dim, dimworld >;

  public:
    static const int dimension = GridImp::dimension;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;
    static const int dimensionworld = GridImp::dimensionworld;

    typedef typename GridImp::template Codim< codimension >::Entity Entity;
    typedef MakeableInterfaceObject< Entity > EntityObject;
    typedef typename EntityObject::ImplementationType EntityImp;

    typedef AlbertaGridEntityPointer< codimension, GridImp > EntityPointerImp;

    typedef Alberta::ElementInfo< dim > ElementInfo;

    //! Constructor for EntityPointer that points to an element
    AlbertaGridEntityPointer ( const GridImp &grid,
                               int level,
                               const ElementInfo &elementInfo,
                               int subEntity );

    //! Constructor for EntityPointer init of Level- and LeafIterator
    AlbertaGridEntityPointer( const GridImp &grid, int level, bool done );

    //! make entity pointer from entity
    AlbertaGridEntityPointer ( const EntityImp &entity );

    //! make empty entity pointer (to be revised)
    AlbertaGridEntityPointer ( const This &other );

    //! make empty entity pointer (to be revised)
    AlbertaGridEntityPointer(const GridImp & , const EntityImp & en);

    //! Destructor
    ~AlbertaGridEntityPointer();

    //! assignment operator
    This &operator= ( const This &other );

    //! equality
    bool equals ( const This &other ) const;

    //! dereferencing
    Entity &dereference () const;

    //! ask for level of entities
    int level () const;

    //! has to be called when iterator is finished
    void done ();

    //! reduce memory
    void compactify ()
    {}

  protected:
    //! return reference to internal entity imp
    EntityImp &entityImp ();

    //! return const reference to internal entity imp
    const EntityImp &entityImp () const;

    //! obtain a reference to the grid
    const GridImp &grid () const;

  private:
    // entity that this EntityPointer points to
    EntityObject *entity_;
  };

}

#endif
