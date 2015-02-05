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
  template<int codim, class GridImp, class HostGridEntityPointer_>
  class IdentityGridEntityPointer
  {
  private:

    enum { dim = GridImp::dimension };

    template<int, typename, typename>
    friend class IdentityGridEntityPointer;


  public:

    //! export the type of the EntityPointer Implementation.
    //! Necessary for the typeconversion between Iterators and EntityPointer
    typedef IdentityGridEntityPointer EntityPointerImp;

    /** \brief Codimension of entity pointed to */
    enum { codimension = codim };

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    // The codimension of this entitypointer wrt the host grid
    enum {CodimInHostGrid = GridImp::HostGridType::dimension - GridImp::dimension + codim};

    // EntityPointer to the equivalent entity in the host grid
    typedef HostGridEntityPointer_ HostGridEntityPointer;


    //! constructor
    IdentityGridEntityPointer (const GridImp* identityGrid, const HostGridEntityPointer& hostEntityPointer)
      : identityGrid_(identityGrid)
      , hostEntityPointer_(hostEntityPointer)
    {}

    ///! copy constructor from EntityPointer storing different host EntityPointer
    template<typename ForeignHostGridEntityPointer>
    explicit IdentityGridEntityPointer (const IdentityGridEntityPointer<codim,GridImp,ForeignHostGridEntityPointer>& entityPointer)
      : identityGrid_(entityPointer.identityGrid_)
      , hostEntityPointer_(entityPointer.hostEntityPointer_)
    {}

    ///! assignment operator from EntityPointer storing different host EntityPointer
    template<typename ForeignHostGridEntityPointer>
    IdentityGridEntityPointer& operator=(const IdentityGridEntityPointer<codim,GridImp,ForeignHostGridEntityPointer>& entityPointer)
    {
      hostEntityPointer_ = entityPointer.hostEntityPointer_;
      return *this;
    }

    //! Move constructor to avoid copying the host EntityPointer
    IdentityGridEntityPointer (const GridImp* identityGrid, HostGridEntityPointer&& hostEntityPointer)
      : identityGrid_(identityGrid)
      , hostEntityPointer_(std::move(hostEntityPointer))
    {}

    //! Constructor from an IdentityGrid entity
    IdentityGridEntityPointer (const IdentityGridEntity<codim,dim,GridImp>& entity)
      : identityGrid_(entity.identityGrid_)
      , hostEntityPointer_(entity.hostEntity_)
    {}

    //! equality
    bool equals(const IdentityGridEntityPointer& i) const {
      return hostEntityPointer_ == i.hostEntityPointer_;
    }

    //! equality with EntityPointer based on different host EntityPointer
    template<typename ForeignHostGridEntityPointer>
    bool equals(const IdentityGridEntityPointer<codim,GridImp,ForeignHostGridEntityPointer>& entityPointer) const
    {
      return  hostEntityPointer_ == entityPointer.hostEntityPointer_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity{{identityGrid_,*hostEntityPointer_}};
    }

    //! Make this pointer as small as possible
    void compactify () {
      //virtualEntity_.getTarget().compactify();
    }

    //! ask for level of entity
    int level () const {
      return hostEntityPointer_->level();
    }


  protected:

    const GridImp* identityGrid_;

    //! host EntityPointer
    HostGridEntityPointer hostEntityPointer_;


  };


} // end namespace Dune

#endif
