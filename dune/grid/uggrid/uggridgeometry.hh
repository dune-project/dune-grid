// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRIDGEOMETRY_HH
#define DUNE_UGGRIDGEOMETRY_HH

/** \file
 * \brief The UGGridGeometry class and its specializations
 */

#include "uggridrenumberer.hh"
#include <dune/common/array.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/genericgeometry/geometry.hh>

namespace Dune {


  /** \brief Defines the geometry part of a mesh entity.
   * \ingroup UGGrid

     \tparam mydim Dimension of the corresponding reference element
     \tparam coorddim Each corner is a point with coorddim coordinates.

     This version is actually used only for mydim==coorddim. The mydim<coorddim cases
     are in specializations below.
   */
  template<int mydim, int coorddim, class GridImp>
  class UGGridGeometry :
    public GeometryDefaultImplementation <mydim, coorddim, GridImp, UGGridGeometry>
  {
    typedef typename GridImp::ctype UGCtype;

    template <int codim_, int dim_, class GridImp_>
    friend class UGGridEntity;

  public:

    /** Default constructor. Puts geometry in element mode
     */
    UGGridGeometry()
    {
      mode_ = element_mode;
      jacobianIsUpToDate_ = false;
      jacobianInverseIsUpToDate_ = false;
    }

    //! put object in coord_mode
    void coordmode ()
    {
      // set the mode
      mode_ = coord_mode;

      // initialize pointers to data
      for (int i=0; i<((mydim==2) ? 4 : 8); i++)
        cornerpointers_[i] = &(coord_[i][0]);

      jacobianIsUpToDate_ = false;
      jacobianInverseIsUpToDate_ = false;
    }

    /** \brief Return the element type identifier
     *
     * UGGrid supports triangles and quadrilaterals in 2D, and
     * tetrahedra, pyramids, prisms, and hexahedra in 3D.
     */
    GeometryType type () const;

    //! returns true if type is simplex, false otherwise (impl could be improved)
    bool affine() const { return type().isSimplex(); }

    //! return the number of corners of this element.
    int corners () const {
      return UG_NS<coorddim>::Corners_Of_Elem(target_);
    }

    //! access to coordinates of corners. Index is the number of the corner
    FieldVector<UGCtype, coorddim> corner (int i) const;

    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    FieldVector<UGCtype, coorddim> global (const FieldVector<UGCtype, mydim>& local) const;

    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    FieldVector<UGCtype, mydim> local (const FieldVector<UGCtype, coorddim>& global) const;

    /**
       Integration over a general element is done by integrating over the reference element
       and using the transformation from the reference element to the global element as follows:
       \f[\int\limits_{\Omega_e} f(x) dx = \int\limits_{\Omega_{ref}} f(g(l)) A(l) dl \f] where
       \f$g\f$ is the local to global mapping and \f$A(l)\f$ is the integration element.

       For a general map \f$g(l)\f$ involves partial derivatives of the map (surface element of
       the first kind if \f$d=2,w=3\f$, determinant of the Jacobian of the transformation for
       \f$d=w\f$, \f$\|dg/dl\|\f$ for \f$d=1\f$).

       For linear elements, the derivatives of the map with respect to local coordinates
       do not depend on the local coordinates and are the same over the whole element.

       For a structured mesh where all edges are parallel to the coordinate axes, the
       computation is the length, area or volume of the element is very simple to compute.

       Each grid module implements the integration element with optimal efficiency. This
       will directly translate in substantial savings in the computation of finite element
       stiffness matrices.
     */
    UGCtype integrationElement (const FieldVector<UGCtype, mydim>& local) const;

    UGCtype volume() const {

      if (mydim==0)
        return 1;

      if (mode_==element_mode) {
        // coorddim*coorddim is an upper bound for the number of vertices
        UGCtype* cornerCoords[coorddim*coorddim];
        UG_NS<coorddim>::Corner_Coordinates(target_, cornerCoords);
        return UG_NS<coorddim>::Area_Of_Element(corners(),
                                                const_cast<const double**>(cornerCoords));
      } else {
        return UG_NS<coorddim>::Area_Of_Element(corners(),
                                                const_cast<const double**>(cornerpointers_));
      }
    }

    //! The inverse transpose of the Jacobian matrix of the mapping from the reference element to this element
    const FieldMatrix<UGCtype, coorddim,mydim>& jacobianInverseTransposed (const FieldVector<UGCtype, mydim>& local) const;
    //! The transpose of the Jacobian matrix of the mapping from the reference element to this element
    const FieldMatrix<UGCtype, mydim,coorddim>& jacobianTransposed (const FieldVector<UGCtype, mydim>& local) const;


  private:
    // mode that decides whether coordinates are taken from the element or given explicitely
    enum SourceMode {element_mode, coord_mode};

    // mode is set by constructor
    SourceMode mode_;

    /** \brief Init the element with a given UG element */
    void setToTarget(typename UG_NS<coorddim>::template Entity<coorddim-mydim>::T* target)
    {
      target_ = target;
      jacobianIsUpToDate_ = false;
      jacobianInverseIsUpToDate_ = false;
    }

    //! \brief set a corner
    void setCoords (int i, const UGCtype* pos)
    {
      if (mode_!=coord_mode)
        DUNE_THROW(GridError,"mode must be coord_mode!");

      for (int j=0; j<coorddim; j++)
        coord_[i][j] = pos[j];

      jacobianIsUpToDate_ = false;
      jacobianInverseIsUpToDate_ = false;
    }

    //! the vertex coordinates
    mutable array<FieldVector<UGCtype, coorddim>, (mydim==2) ? 4 : 8> coord_;

    //! The jacobian inverse transposed
    mutable FieldMatrix<UGCtype,coorddim,mydim> jac_inverse_;
    //! The jacobian transposed
    mutable FieldMatrix<UGCtype,mydim,coorddim> jac_;

    /** \brief For simplices the Jacobian matrix is a constant, hence it needs
        to be computed only once for each new element.  This can save some
        assembly time. */
    mutable bool jacobianInverseIsUpToDate_;
    mutable bool jacobianIsUpToDate_;

    // in element mode this points to the element we map to
    // in coord_mode this is the element whose reference element is mapped into the father's one
    typename UG_NS<coorddim>::template Entity<coorddim-mydim>::T* target_;

    // in coord_mode we explicitely store an array of coordinates
    // containing the position in the fathers reference element
    mutable UGCtype* cornerpointers_[(mydim==2) ? 4 : 8];
  };



  /****************************************************************/
  /*                                                              */
  /*       Specialization for faces in 3d                         */
  /*                                                              */
  /****************************************************************/

  template<class GridImp>
  class UGGridGeometry<2, 3, GridImp> :
    public GenericGeometry::BasicGeometry<2, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,2,3> >
  {

    template <int codim_, int dim_, class GridImp_>
    friend class UGGridEntity;

    template <class GridImp_>
    friend class UGGridIntersectionIterator;

    typedef typename GenericGeometry::BasicGeometry<2, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,2,3> > Base;

  public:

    /** \brief Setup method with a geometry type and a set of corners
        \param coordinates The corner coordinates in DUNE numbering
     */
    void setup(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,3> >& coordinates)
    {
      // set up base class
      // Yes, a strange way, but the only way, as BasicGeometry doesn't have a setup method
      static_cast< Base & >( *this ) = Base( type, coordinates );
    }

  };


  /****************************************************************/
  /*                                                              */
  /*       Specialization for faces in 2d                         */
  /*                                                              */
  /****************************************************************/

  template<class GridImp>
  class UGGridGeometry <1, 2, GridImp> :
    public GenericGeometry::BasicGeometry<1, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,1,2> >
  {

    template <int codim_, int dim_, class GridImp_>
    friend class UGGridEntity;

    template <class GridImp_>
    friend class UGGridIntersectionIterator;

    typedef typename GenericGeometry::BasicGeometry<1, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,1,2> > Base;

  public:

    /** \brief Constructor with a geometry type and a set of corners */
    void setup(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,2> >& coordinates)
    {
      // set up base class
      // Yes, a strange way, but the only way, as BasicGeometry doesn't have a setup method
      static_cast< Base & >( *this ) = Base( type, coordinates );
    }

  };

}  // namespace Dune

#endif
