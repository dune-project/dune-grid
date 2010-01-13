// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU_BNDPROJECTION_HH
#define DUNE_ALU_BNDPROJECTION_HH

namespace Dune {

  //! \brief ALUGrid boundary projection implementation
  //!  DuneBndProjection has to fulfil the DuneBoundaryProjection interface
  template <class GridImp>
  class ALUGridBoundaryProjection
#ifdef ALUGRID_VERTEX_PROJECTION
    : public GridImp :: ALUGridVertexProjectionType
#endif
  {
    typedef GridImp GridType;
    // type of double coordinate vector
    typedef double coord_t[ GridType :: dimension ];
  protected:

    //! reference to boundary projection implementation
    const GridType& grid_;
  public:
    //! type of boundary projection
    typedef typename GridType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

    //! type of coordinate vector
    typedef typename DuneBoundaryProjectionType :: CoordinateType CoordinateType;

    //! constructor storing reference to boundary projection implementation
    ALUGridBoundaryProjection(const GridType& grid)
      : grid_( grid )
    {}

    //! (old) method projection vertices defaults to segment 0
    int operator () (const coord_t &orig,
                     coord_t &prj) const
    {
      return this->operator()( orig, 0, prj);
    }

    //! projection operator
    int operator () (const coord_t &orig,
                     const int segmentIndex,
                     coord_t &prj) const
    {
#ifdef ALUGRID_VERTEX_PROJECTION
      // get boundary projection
      const DuneBoundaryProjectionType& bndPrj = grid_.boundaryProjection( segmentIndex );

      // call projection operator
      reinterpret_cast<CoordinateType &> (* (&prj[0])) =
        bndPrj( reinterpret_cast<const CoordinateType &> (* (&orig[0])) );
#endif

      // return 1 for success
      return 1;
    }
  };

} // end namespace Dune
#endif
