// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_UTILITY_MOCKGEOMETRY_HH
#define DUNE_GRID_UTILITY_MOCKGEOMETRY_HH

#include <cstddef>

#include <dune/common/fmatrix.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/genericgeometry/geometry.hh>
#include <dune/grid/genericgeometry/geometrytraits.hh>

namespace Dune {

  //! Grid-independent geometry
  /**
   * This geometry can be used when somewhing mostly like a Dune::Geometry is
   * required, but a full grid is a little bit too much.  It provides the full
   * interface of Dune::Geometry, except for the grid-specific member
   * constants \c dimension and \c dimensionworld.
   *
   * One further difference is that the jacobian methods return by value
   * instead of by reference.  The Jacobian depends on the local coordinate;
   * returning it by reference is asking for trouble.
   *
   * \tparam ctype    Field type for coordinates.
   * \tparam mydim    Dimension of the local coordinates.
   * \tparam coorddim Dimension of the global coordinates.
   */
  template<class ctype, std::size_t mydim, std::size_t coorddim>
  class MockGeometry :
    public GenericGeometry::BasicGeometry<
        mydim, GenericGeometry::DefaultGeometryTraits<ctype, coorddim, coorddim>
        >
  {
    typedef GenericGeometry::DefaultGeometryTraits<ctype, coorddim, coorddim>
    Traits;
    typedef GenericGeometry::BasicGeometry<mydim, Traits> Base;

    // Hide members of BasicGeometry that are not part of Dune::Geometry
    using Base::JacobianInverseTransposed;
    using Base::operator[];
    using Base::normal;

  public:
    typedef typename Base::Allocator Allocator;

    //! type of jacobian (also of jacobian inverse transposed)
    typedef FieldMatrix<ctype, coorddim, mydim> Jacobian;
    //! type of jacobian transposed
    typedef FieldMatrix<ctype, mydim, coorddim> JacobianTransposed;

    //! Default constructor.
    MockGeometry(const Allocator &allocator = Allocator()) : Base(allocator) {}
    //! constructor
    template<class CoordVector>
    MockGeometry(unsigned topologyId, const CoordVector &coords,
                 const Allocator &allocator=Allocator()) :
      Base(topologyId, coords, allocator)
    { }
    //! constructor
    template<class CoordVector>
    MockGeometry(unsigned topologyId, const CoordVector &coords, bool affine,
                 const Allocator &allocator=Allocator()) :
      Base(topologyId, coords, affine, allocator)
    { }
    //! Constructor using a GeometryType and a list of corner coordinates.
    template<class CoordVector>
    MockGeometry(const GeometryType &type, const CoordVector &coords,
                 const Allocator &allocator=Allocator()) :
      Base(type, coords, allocator)
    { }
    //! obtain a geometry for a subentity
    template<int fatherdim>
    MockGeometry(const MockGeometry<ctype, fatherdim, coorddim> &father,
                 int i) :
      Base(static_cast<const GenericGeometry::BasicGeometry<
                         fatherdim,
                         GenericGeometry::DefaultGeometryTraits<ctype, coorddim, coorddim>
                         > &>(father),
           i)
    { }

    //! Return the transposed of the Jacobian.
    JacobianTransposed jacobianTransposed
      (const typename Base::LocalCoordinate &local) const
    { return Base::jacobianTransposed(local); }
    //! Return inverse of transposed of Jacobian.
    Jacobian jacobianInverseTransposed
      (const typename Base::LocalCoordinate &local) const
    { return Base::jacobianInverseTransposed(local); }
  };

}  // namespace Dune

#endif // DUNE_GRID_UTILITY_MOCKGEOMETRY_HH
