// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDGEOMETRY_HH
#define DUNE_GRID_YASPGRIDGEOMETRY_HH

/** \file
 * \brief The YaspGeometry class and its specializations

   YaspGeometry realizes the concept of the geometric part of a mesh entity.

   We have specializations for dim == dimworld (elements) and dim == 0
   (vertices).  The general version implements all other codimensions.

   As of September 2014, the functionality of YaspGeometry is identical
   to that of AxisAlignedCubeGeometry. The latter cannot be used directly
   due to the grid interface facade construction (it isn't templated to the
   GridImp). As soon as template aliases are available, this header boils
   down to one line.
 */

namespace Dune {

  /** \brief The general version that handles all codimensions but 0 and dim.
   * \tparam mydim the codimension
   * \tparam cdim the coordinate dimension (NOT codim)
   */
  template<int mydim,int cdim, class GridImp>
  class YaspGeometry : public AxisAlignedCubeGeometry<typename GridImp::ctype,mydim,cdim>
  {
  public:
    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

    //! constructor from midpoint and extension and a bitset defining which unit vectors span the entity
    YaspGeometry (const FieldVector<ctype, cdim>& ll, const FieldVector<ctype, cdim>& ur, const std::bitset<cdim>& shift)
      : AxisAlignedCubeGeometry<ctype,mydim,cdim>(ll,ur,shift)
    {
      assert(mydim == shift.count());
    }
  };

  //! specialize for dim=dimworld, i.e. a volume element
  template<int mydim, class GridImp>
  class YaspGeometry<mydim,mydim,GridImp> : public AxisAlignedCubeGeometry<typename GridImp::ctype,mydim,mydim>
  {
  public:
    typedef typename GridImp::ctype ctype;

    //! constructor from midpoint and extension
    YaspGeometry (const FieldVector<ctype, mydim>& ll, const FieldVector<ctype, mydim>& ur)
      : AxisAlignedCubeGeometry<ctype,mydim,mydim>(ll,ur)
    {}

    //! copy constructor (skipping temporary variables)
    YaspGeometry (const YaspGeometry& other)
      : AxisAlignedCubeGeometry<ctype,mydim,mydim>(other)
    {}
  };

  //! specialization for dim=0, this is a vertex
  template<int cdim, class GridImp>
  class YaspGeometry<0,cdim,GridImp> : public AxisAlignedCubeGeometry<typename GridImp::ctype,0,cdim>
  {
  public:
    typedef typename GridImp::ctype ctype;

    //! constructor
    explicit YaspGeometry ( const FieldVector< ctype, cdim > &p )
      : AxisAlignedCubeGeometry<typename GridImp::ctype,0,cdim>( p )
    {}

    YaspGeometry ( const FieldVector< ctype, cdim > &p, const FieldVector< ctype, cdim > &, const std::bitset<cdim> &)
      : AxisAlignedCubeGeometry<typename GridImp::ctype,0,cdim>( p )
    {}
  };
}  // namespace Dune

#endif // DUNE_GRID_YASPGRIDGEOMETRY_HH
