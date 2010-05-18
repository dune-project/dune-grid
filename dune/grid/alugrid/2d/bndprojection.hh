// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2D_BNDPROJECTION_HH
#define DUNE_ALU2D_BNDPROJECTION_HH

#include <dune/grid/alugrid/bndprojection.hh>

#include <dune/grid/alugrid/2d/alu2dinclude.hh>

namespace Dune
{

#ifdef ALUGRID_SURFACE_2D

  template< class Grid >
  class ALU2dGridBoundaryProjection
    : public ALU2DSPACE VtxProjection ALU2DDIMWORLD(Grid::dimensionworld,Grid::elementType)
  {
    typedef ALU2DSPACE VtxProjection ALU2DDIMWORLD (Grid::dimensionworld,Grid::elementType) Base;

  public:
    enum { ncoord = Base::ncoord };

    typedef typename Base::hbndel_t hbndel_t;
    typedef typename Base::helement_t helement_t;

    typedef typename Grid::DuneBoundaryProjectionType DuneBoundaryProjectionType;

    typedef typename DuneBoundaryProjectionType::CoordinateType CoordinateType;

    explicit ALU2dGridBoundaryProjection ( const Grid &grid )
      : grid_( grid )
    {}

    int operator() ( const hbndel_t *hbndel, const double local, double (&global)[ ncoord ] ) const
    {
      return callProjection( grid_.boundaryProjection( hbndel->segmentIndex() ), global );
    }

    int operator() ( const helement_t *helement, const double (&local)[ 2 ], double (&global)[ ncoord ] ) const
    {
      return callProjection( grid_.globalProjection(), global );
    }

  private:
    static int callProjection ( const DuneBoundaryProjectionType *prj, double (&global)[ ncoord ] )
    {
      if( prj )
      {
        CoordinateType x, y;
        for( int i = 0; i < ncoord; ++i )
          x[ i ] = global[ i ];
        y = (*prj)( x );
        for( int i = 0; i < ncoord; ++i )
          global[ i ] = y[ i ];
      }
      return 1;
    }

    const Grid &grid_;
  };

#else // #ifdef ALUGRID_SURFACE_2D

  template< class Grid >
  class ALU2dGridBoundaryProjection
    : public ALUGridBoundaryProjection< Grid >
  {
    typedef ALUGridBoundaryProjection< Grid > Base;

  public:
    explicit ALU2dGridBoundaryProjection ( const Grid &grid )
      : Base( grid )
    {}
  };

#endif // #else // #ifdef ALUGRID_SURFACE_2D

} // end namespace Dune

#endif // #ifndef DUNE_ALU2D_BNDPROJECTION_HH
