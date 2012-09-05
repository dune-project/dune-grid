// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_ENTITY_POINTER_HH
#define DUNE_UGGRID_ENTITY_POINTER_HH


namespace Dune {

  /*! Acts as a pointer to an  entities of a given codimension.
   */
  template<int codim, class GridImp>
  class UGGridEntityPointer
  {
    enum { dim = GridImp::dimension };
  public:
    //! export the type of the EntityPointer Implementation.
    //! Necessary for the type conversion between Iterators and EntityPointer
    typedef UGGridEntityPointer EntityPointerImp;

    //! codimension of entity pointer
    enum { codimension = codim };

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! constructor
    UGGridEntityPointer ()  {
      virtualEntity_.setToTarget(nullptr,nullptr);
    }

    //! constructor
    UGGridEntityPointer (typename UG_NS<dim>::template Entity<codim>::T* target, const GridImp* gridImp)
      : virtualEntity_(target, gridImp)
    {}

    //! construct entity pointer from given entity
    UGGridEntityPointer (const UGGridEntity<codim,dim,GridImp>& entity)
      : virtualEntity_(entity.target_,entity.gridImp_)
    {}

    /** \brief Construct face entity pointer from given element and side number
     *
     * Only used if dim==3 and codim==1
     */
    UGGridEntityPointer (const UGGridEntity<0,dim,GridImp>& centerEntity, unsigned int side)
      : virtualEntity_(centerEntity.target_, side)
    {}

    void setToTarget(typename UG_NS<dim>::template Entity<codim>::T* target, const GridImp* gridImp) {
      virtualEntity_.setToTarget(target,gridImp);
    }

    typename UG_NS<dim>::template Entity<codim>::T* getTarget()
    {
      return virtualEntity_.getTarget();
    }

    const typename UG_NS<dim>::template Entity<codim>::T* getTarget() const
    {
      return virtualEntity_.getTarget();
    }

    //! equality
    bool equals(const UGGridEntityPointer<codim,GridImp>& i) const {
      return getTarget() == i.getTarget();
    }

    //! dereferencing
    Entity& dereference() const {
      return virtualEntity_;
    }

    //! ask for level of entity
    int level () const {return virtualEntity_.level();}

  protected:

    mutable UGMakeableEntity<codim,dim,GridImp> virtualEntity_; //!< virtual entity

  };


} // end namespace Dune

#endif
