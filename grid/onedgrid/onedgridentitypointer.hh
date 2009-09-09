// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONEDGRID_ENTITY_POINTER_HH
#define DUNE_ONEDGRID_ENTITY_POINTER_HH


namespace Dune {

  /*! Acts as a pointer to an  entities of a given codimension.
   */
  template<int codim, class GridImp>
  class OneDGridEntityPointer
  {
    enum { dim = GridImp::dimension };
    template <class GridImp_>
    friend class OneDGridLevelIntersectionIterator;
    template <class GridImp_>
    friend class OneDGridLeafIntersectionIterator;
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
      return other.virtualEntity_.target() == virtualEntity_.target();
    }

    //! dereferencing
    Entity& dereference() const {return virtualEntity_;}

    //! ask for level of entity
    int level () const {return virtualEntity_.level();}

    OneDGridEntityPointer() {}

    /** \brief Constructor from a given entity  */
    OneDGridEntityPointer(const OneDGridEntity<codim, dim, GridImp> & entity)
    {
      virtualEntity_.setToTarget(entity.target_);
    }

    //! empty method since internal entity is not a pointer
    void compactify () {}

  protected:

    /** \brief Constructor from a given iterator */
    OneDGridEntityPointer(OneDEntityImp<dim-codim>* it) {
      virtualEntity_.setToTarget(it);
    };

  protected:

    mutable OneDEntityWrapper<codim,GridImp::dimension,GridImp> virtualEntity_;

  };


} // end namespace Dune

#endif
