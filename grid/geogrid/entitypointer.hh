// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITYPOINTER_HH
#define DUNE_GEOGRID_ENTITYPOINTER_HH

namespace Dune
{

  /** Acts as a pointer to an  entities of a given codimension.
   */
  template<int codim, class GridImp>
  class GeometryGridEntityPointer :
    public EntityPointerDefaultImplementation <codim, GridImp, Dune::GeometryGridEntityPointer<codim,GridImp> >
  {
  private:

    enum { dim = GridImp::dimension };


  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    typedef GeometryGridEntityPointer<codim,GridImp> Base;

    // The codimension of this entitypointer wrt the host grid
    enum {CodimInHostGrid = GridImp::HostGridType::dimension - GridImp::dimension + codim};

    // EntityPointer to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::EntityPointer HostGridEntityPointer;


    //! constructor
    GeometryGridEntityPointer (const GridImp* identityGrid, const HostGridEntityPointer& hostEntity_) :
      identityGrid_(identityGrid),
      virtualEntity_(identityGrid, hostEntity_)
    {}


    //! equality
    bool equals(const GeometryGridEntityPointer<codim,GridImp>& i) const {
      return virtualEntity_.getTarget() == i.virtualEntity_.getTarget();
    }


    //! dereferencing
    Entity& dereference() const {
      return virtualEntity_;
    }


    //! ask for level of entity
    int level () const {
      return virtualEntity_.level();
    }


  protected:

    const GridImp* identityGrid_;

    //! virtual entity
    mutable GeometryGridMakeableEntity<codim,dim,GridImp> virtualEntity_;


  };


} // end namespace Dune

#endif
