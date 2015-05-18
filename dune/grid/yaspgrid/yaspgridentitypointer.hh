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
    typedef typename GridImp::YGrid::Iterator I;
    typedef YaspEntityPointer<codim,GridImp> EntityPointerImp;
  protected:
    typedef YaspEntity<codim, dim, GridImp> YaspEntityImp;

  public:
    //! codimension of entity pointer
    enum { codimension = codim };

    //! default constructor
    YaspEntityPointer () :
      _entity(YaspEntityImp())
    {}

    //! constructor
    YaspEntityPointer (const YGLI & g, const I& it)
      : _entity(YaspEntityImp(g,it))
    {}


// skip this constructor for GCC 4.4, which has a number of nasty bugs in its rvalue reference support
// As this behavior is hard to trigger in small configuration tests and because we'll probably drop GCC 4.4
// after the next release anyway, I hacked in this hardcoded check for the compiler version
#if not (defined(__GNUC__) && (__GNUC__ < 5) && (__GNUC_MINOR__ < 5))

    YaspEntityPointer (YGLI&& g, I&& it)
      : _entity(YaspEntityImp(std::move(g),std::move(it)))
    {}

#endif

    //! copying and moving
    YaspEntityPointer (const YaspEntityImp& entity)
      : _entity(entity)
    {}

    YaspEntityPointer (YaspEntityImp&& entity)
      : _entity(std::move(entity))
    {}

    //! copying and moving -- use default implementations

    //! equality
    bool equals (const YaspEntityPointer& rhs) const
    {
      return (_entity == rhs._entity);
    }

    //! dereferencing
    const Entity& dereference() const
    {
      return _entity;
    }

    //! ask for level of entity
    int level () const {return _entity.level();}

    //! use default assignment operator

  protected:
    Entity _entity; //!< entity
  };

}   // namespace Dune

#endif   // DUNE_GRID_YASPGRIDENTITYPOINTER_HH
