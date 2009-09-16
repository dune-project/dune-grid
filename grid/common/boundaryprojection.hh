// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BOUNDARYPROJECTION_HH
#define DUNE_BOUNDARYPROJECTION_HH

//- system includes
#include <cmath>

//- Dune includes
#include <dune/common/fvector.hh>

namespace Dune {

  /** \brief Interface class for vertex projection at the boundary.
   */
  template <int dimworld>
  class DuneBoundaryProjection
  {
  public:
    //! \brief type of coordinate vector
    typedef FieldVector< double, dimworld> CoordinateType;
  protected:
    //! default constructor
    DuneBoundaryProjection() {}
  public:
    //! \brief destructor
    virtual ~DuneBoundaryProjection() {}

    //! \brief projection operator projection a global coordinate
    virtual CoordinateType operator() (const CoordinateType& global) const = 0;
  };

  //////////////////////////////////////////////////////////////////////
  //
  // Example of boundary projection projection to a circle
  //
  //////////////////////////////////////////////////////////////////////
  template <int dimworld>
  class CircleBoundaryProjection : public DuneBoundaryProjection< dimworld >
  {
  protected:
    //! radius of circ
    const double radius_;

  public:
    //! \brief type of coordinate vector
    typedef FieldVector< double, dimworld> CoordinateType;

    //! constructor taking radius of circle (default = sqrt( dimworld ) )
    CircleBoundaryProjection(const double radius = std::sqrt( (double)dimworld ))
      : radius_( radius ) {}

    //! \brief destructor
    virtual ~CircleBoundaryProjection() {}

    //! \brief projection operator projection a global coordinate
    virtual CoordinateType operator() (const CoordinateType& global) const
    {
      CoordinateType prj( global );
      // get adjustment factor
      const double factor = radius_  / global.two_norm();
      // adjust
      prj *= factor;
      return prj;
    }
  };

} // end namespace
#endif
