// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_LEAFITERATOR_HH
#define DUNE_GEOGRID_LEAFITERATOR_HH

namespace Dune
{

  /** \brief Iterator over all entities of a given codimension and level of a grid.
   *  \ingroup GeometryGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class GeometryGridLeafIterator :
    public Dune::GeometryGridEntityPointer <codim,GridImp>
  {
  private:

    enum {dim = GridImp::dimension};


  public:

    //! \todo Please doc me !
    explicit GeometryGridLeafIterator(const GridImp* identityGrid) :
      GeometryGridEntityPointer<codim,GridImp>(identityGrid, identityGrid->hostgrid_->template leafbegin<codim>()),
      hostGridLeafIterator_(identityGrid->hostgrid_->template leafbegin<codim>()),
      hostGridLeafEndIterator_(identityGrid->hostgrid_->template leafend<codim>())
    {
      this->virtualEntity_.setToTarget(hostGridLeafIterator_);
    }


    /** \brief Constructor which create the end iterator
     *  \param endDummy Here only to distinguish it from the other constructor
     */
    explicit GeometryGridLeafIterator(const GridImp* identityGrid, bool endDummy) :
      GeometryGridEntityPointer<codim,GridImp>(identityGrid, identityGrid->hostgrid_->template leafend<codim>()),
      hostGridLeafIterator_(identityGrid->hostgrid_->template leafbegin<codim>()),
      hostGridLeafEndIterator_(identityGrid->hostgrid_->template leafend<codim>())
    {}


    //! prefix increment
    void increment() {
      ++hostGridLeafIterator_;
      this->virtualEntity_.setToTarget(hostGridLeafIterator_);
    }


  private:

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////

    // LevelIterator to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::template Codim<codim>::LeafIterator HostGridLeafIterator;

    //! \todo Please doc me !
    HostGridLeafIterator hostGridLeafIterator_;

    //! \todo Please doc me !
    HostGridLeafIterator hostGridLeafEndIterator_;

  };


}  // namespace Dune

#endif
