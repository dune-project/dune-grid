// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRID_ENTITY_POINTER_HH
#define DUNE_IDENTITYGRID_ENTITY_POINTER_HH

#include "identitygridentity.hh"

/** \file
 * \brief The IdentityGridEntityPointer class
 */

namespace Dune {


  /** Acts as a pointer to an  entities of a given codimension.
   */
  template<int codim, class GridImp>
  class IdentityGridEntityPointer
  {
  private:

    enum { dim = GridImp::dimension };


  public:

    //! export the type of the EntityPointer Implementation.
    //! Necessary for the typeconversion between Iterators and EntityPointer
    typedef IdentityGridEntityPointer EntityPointerImp;

    /** \brief Codimension of entity pointed to */
    enum { codimension = codim };

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    typedef IdentityGridEntityPointer<codim,GridImp> Base;

    // The codimension of this entitypointer wrt the host grid
    enum {CodimInHostGrid = GridImp::HostGridType::dimension - GridImp::dimension + codim};

    // EntityPointer to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::EntityPointer HostGridEntityPointer;


    //! constructor
    template< class HostGridEntityPointer >
    IdentityGridEntityPointer (const GridImp* identityGrid, const HostGridEntityPointer& hostEntity_) :
      identityGrid_(identityGrid),
      virtualEntity_(identityGrid, hostEntity_)
    {}

    //! Constructor from an IdentityGrid entity
    IdentityGridEntityPointer (const IdentityGridEntity<codim,dim,GridImp>& entity)
      : identityGrid_(entity.identityGrid_),
        virtualEntity_(entity.identityGrid_, entity.hostEntity_)
    {}

    //! equality
    bool equals(const IdentityGridEntityPointer<codim,GridImp>& i) const {
      return virtualEntity_.getTarget() == i.virtualEntity_.getTarget();
    }


    //! dereferencing
    Entity& dereference() const {
      return virtualEntity_;
    }

    //! Make this pointer as small as possible
    void compactify () {
      //virtualEntity_.getTarget().compactify();
    }

    //! ask for level of entity
    int level () const {
      return virtualEntity_.level();
    }


  protected:

    const GridImp* identityGrid_;

    //! virtual entity
    mutable IdentityGridMakeableEntity<codim,dim,GridImp> virtualEntity_;


  };


} // end namespace Dune

#endif
