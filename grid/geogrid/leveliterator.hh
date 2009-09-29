// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_LEVELITERATOR_HH
#define DUNE_GEOGRID_LEVELITERATOR_HH

namespace Dune
{

  //**********************************************************************
  //
  // --GeometryGridLevelIterator
  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup GeometryGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class GeometryGridLevelIterator :
    public Dune::GeometryGridEntityPointer <codim,GridImp>,
    public LevelIteratorDefaultImplementation <codim,pitype,GridImp,GeometryGridLevelIterator>
  {
  private:

    enum {dim = GridImp::dimension};


  public:

    //! Constructor
    explicit GeometryGridLevelIterator(const GridImp* identityGrid, int level)
      : GeometryGridEntityPointer<codim,GridImp>(identityGrid, identityGrid->hostgrid_->template lbegin<codim>(level)),
        hostGridLevelIterator_(identityGrid->hostgrid_->template lbegin<codim>(level)),
        hostGridLevelEndIterator_(identityGrid->hostgrid_->template lend<codim>(level))
    {
      this->virtualEntity_.setToTarget(hostGridLevelIterator_);
    }


    /** \brief Constructor which create the end iterator
        \param endDummy Here only to distinguish it from the other constructor
     */
    explicit GeometryGridLevelIterator(const GridImp* identityGrid, int level, bool endDummy)
      :
        GeometryGridEntityPointer<codim,GridImp>(identityGrid, identityGrid->hostgrid_->template lend<codim>(level)),
        hostGridLevelIterator_(identityGrid->hostgrid_->template lend<codim>(level)),
        hostGridLevelEndIterator_(identityGrid->hostgrid_->template lend<codim>(level))
    {}


    //! prefix increment
    void increment() {
      ++hostGridLevelIterator_;
      this->virtualEntity_.setToTarget(hostGridLevelIterator_);
    }


  private:

    // LevelIterator to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::Traits::template Codim<codim>::LevelIterator HostGridLevelIterator;

    //! \todo Please doc me !
    HostGridLevelIterator hostGridLevelIterator_;

    //! \todo Please doc me !
    HostGridLevelIterator hostGridLevelEndIterator_;

  };


}  // namespace Dune

#endif
