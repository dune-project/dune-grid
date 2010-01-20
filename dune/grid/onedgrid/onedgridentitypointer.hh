// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONEDGRID_ENTITY_POINTER_HH
#define DUNE_ONEDGRID_ENTITY_POINTER_HH

#include "onedgridentity.hh"

namespace Dune {

  /*! Acts as a pointer to an  entities of a given codimension.
   */
  template<int codim, class GridImp>
  class OneDGridEntityPointer
  {
    enum { dim = GridImp::dimension };
    template <class GridImp_>
    friend class OneDGridLevelIntersection;
    template <class GridImp_>
    friend class OneDGridLeafIntersection;
    friend class OneDGridEntity<0,dim,GridImp>;

  public:
    typedef typename GridImp::template Codim<codim>::Entity Entity;
    typedef OneDGridEntityPointer<codim,GridImp> Base;

    /** \brief The type of the class itself
        Do we really need this?
     */
    typedef OneDGridEntityPointer<codim,GridImp> EntityPointerImp;

    //! codimension of entity pointer
    enum { codimension = codim };

    //! equality
    bool equals(const OneDGridEntityPointer<codim,GridImp>& other) const {
      return GridImp::getRealImplementation(other.virtualEntity_).target_
             == GridImp::getRealImplementation(virtualEntity_).target_;
    }

    //! dereferencing
    Entity& dereference() const {return virtualEntity_;}

    //! ask for level of entity
    int level () const {return virtualEntity_.level();}

    OneDGridEntityPointer()
      : virtualEntity_(OneDGridEntity<codim, dim, GridImp>())
    {}

    /** \brief Constructor from a given entity  */
    OneDGridEntityPointer(const OneDGridEntity<codim, dim, GridImp> & entity)
      : virtualEntity_(entity)
    {}

    //! empty method since internal entity is not a pointer
    void compactify () {}

  protected:

    /** \brief Constructor from a given iterator */
    OneDGridEntityPointer(OneDEntityImp<dim-codim>* it)
      : virtualEntity_(OneDGridEntity<codim, dim, GridImp>())
    {
      GridImp::getRealImplementation(virtualEntity_).setToTarget(it);
    };

  protected:

    mutable MakeableInterfaceObject<Entity> virtualEntity_;

  };


} // end namespace Dune

#endif
