// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDENTITYPOINTER_HH
#define DUNE_GRID_YASPGRIDENTITYPOINTER_HH

/** \file
 * \brief The YaspEntityPointer class
 */

namespace Dune {

  /** \brief A pointer to a YaspGrid::Entity
   */
  template<int codim, class GridImp>
  class YaspEntityPointer
  {
    //! know your own dimension
    enum { dim=GridImp::dimension };
    //! know your own dimension of world
    typedef typename GridImp::ctype ctype;

  public:
    typedef typename GridImp::template Codim<codim>::Entity Entity;
    typedef typename GridImp::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;
    typedef YaspEntityPointer<codim,GridImp> EntityPointerImp;
  protected:
    typedef YaspEntity<codim, dim, GridImp> YaspEntityImp;

  public:
    //! codimension of entity pointer
    enum { codimension = codim };

    //! constructor
    YaspEntityPointer (const GridImp * yg, const YGLI & g, const TSI & it)
      : _g(g), _it(it),
        _entity(MakeableInterfaceObject<Entity>(YaspEntity<codim,dim,GridImp>(yg, _g,_it)))
    {
      if (codim>0 && codim<dim)
      {
        DUNE_THROW(GridError, "YaspEntityPointer: codim not implemented");
      }
    }

    //! copy constructor
    YaspEntityPointer (const YaspEntityImp& entity)
      : _g(entity.gridlevel()),
        _it(entity.transformingsubiterator()),
        _entity(MakeableInterfaceObject<Entity>(YaspEntity<codim,dim,GridImp>(entity.yaspgrid(), _g,_it)))
    {
      if (codim>0 && codim<dim)
      {
        DUNE_THROW(GridError, "YaspEntityPointer: codim not implemented");
      }
    }

    //! copy constructor
    YaspEntityPointer (const YaspEntityPointer& rhs)
      : _g(rhs._g), _it(rhs._it), _entity(MakeableInterfaceObject<Entity>(YaspEntity<codim,dim,GridImp>(GridImp::getRealImplementation(rhs._entity).yaspgrid(),_g,_it)))
    {
      if (codim>0 && codim<dim)
      {
        DUNE_THROW(GridError, "YaspEntityPointer: codim not implemented");
      }
    }

    //! equality
    bool equals (const YaspEntityPointer& rhs) const
    {
      return (_it==rhs._it && _g == rhs._g);
    }

    //! dereferencing
    Entity& dereference() const
    {
      return _entity;
    }

    //! ask for level of entity
    int level () const {return _g->level();}

    const YaspEntityPointer&
    operator = (const YaspEntityPointer& rhs)
    {
      _g = rhs._g;
      _it = rhs._it;
      /* _entity = i._entity
       * is done implicitely, as the entity is completely
       * defined via the iterator it belongs to
       */
      return *this;
    }

    const TSI& transformingsubiterator () const
    {
      return _it;
    }

    const YGLI& gridlevel () const
    {
      return _g;
    }

    TSI& transformingsubiterator ()
    {
      return _it;
    }

    YGLI& gridlevel ()
    {
      return _g;
    }

  protected:
    YGLI _g;             // access to grid level
    TSI _it;             // position in the grid level
    mutable MakeableInterfaceObject<Entity> _entity; //!< virtual entity
  };

}   // namespace Dune

#endif   // DUNE_GRID_YASPGRIDENTITYPOINTER_HH
