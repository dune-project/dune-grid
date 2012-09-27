// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_2D_LOCALGEOMETRY_HH
#define DUNE_ALUGRID_2D_LOCALGEOMETRY_HH

#include <dune/geometry/mapping/affinemapping.hh>

//#include <dune/grid/alugrid/common/matrixhelper.hh>
#include <dune/grid/alugrid/2d/geometry.hh>
#include <dune/grid/common/geometry.hh>

namespace Dune
{

  // ALU2dGridLocalGeometryTraits
  // ----------------------------

  struct ALU2dGridLocalGeometryTraits
  {
    typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< alu2d_ctype > > MatrixHelper;
    //typedef ALUMatrixHelper< alu2d_ctype > MatrixHelper;

    struct UserData {};
  };


  // ALU2dGridLocalGeometry
  // ----------------------

  template< int mydim, int cdim, class Grid >
  class ALU2dGridLocalGeometry
    : public AffineMapping< alu2d_ctype, mydim, cdim, ALU2dGridLocalGeometryTraits >
  {
    typedef AffineMapping< alu2d_ctype, mydim, cdim, ALU2dGridLocalGeometryTraits > Base;

  public:
    template< class CoordVector >
    ALU2dGridLocalGeometry ( GeometryType gt, const CoordVector &coordVector )
      : Base( gt, coordVector )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_ALUGRID_2D_LOCALGEOMETRY_HH
