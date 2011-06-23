// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_LOCALGEOMETRY_HH
#define DUNE_UGGRID_LOCALGEOMETRY_HH

/** \file
 * \brief The UGGridElement class and its specializations
 */

#include "uggridrenumberer.hh"
#include <dune/common/array.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/genericgeometry/geometry.hh>

namespace Dune {


  /** \brief Geometry of an entity embedded in another element, not in the world space
   * \ingroup UGGrid

     \tparam mydim Dimension of the corresponding reference element
     \tparam coorddim Dimension of the coordinate space

   */
  template<int mydim, int coorddim, class GridImp>
  class UGGridLocalGeometry :
    public GenericGeometry::BasicGeometry<mydim, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,mydim,coorddim> >
  {

    typedef typename GenericGeometry::BasicGeometry<mydim, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,mydim,coorddim> > Base;

  public:

    /** \brief Setup method with a geometry type and a set of corners
        \param coordinates The corner coordinates in DUNE numbering
     */
    void setup(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,coorddim> >& coordinates)
    {
      // set up base class
      // Yes, a strange way, but the only way, as BasicGeometry doesn't have a setup method
      static_cast< Base & >( *this ) = Base( type, coordinates );
    }

  };

}  // namespace Dune

#endif
