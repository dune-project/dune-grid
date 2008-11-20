// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GEOMETRY
#define DUNE_ALBERTA_GEOMETRY

namespace Dune
{

  // AlbertaGridGeometry
  /*!
     Defines the geometry part of a mesh entity. Works for all dimensions, element types and dime
     of world. Provides reference element and mapping between local and global coordinates.
     The element may have different implementations because the mapping can be
     done more efficient for structured meshes than for unstructured meshes.

     dim: An element is a polygonal in a hyperplane of dimension dim. 0 <= dim <= 3 is typically
     dim=0 is a point.

     dimworld: Each corner is a point with dimworld coordinates.
   */
  template <int mydim, int cdim, class GridImp>
  class AlbertaGridGeometry
    : public GeometryDefaultImplementation<mydim,cdim,GridImp,AlbertaGridGeometry>
  {

    typedef AlbertaGridGeometry<mydim,cdim,GridImp> ThisType;

    //! know dimension of barycentric coordinates
    enum { dimbary=mydim+1};
  public:
    //! Default constructor
    AlbertaGridGeometry();

    //! constructor building geometry in father
    AlbertaGridGeometry(const int child, const int orientation );

    //! return the element type identifier
    //! line , triangle or tetrahedron, depends on dim
    const GeometryType & type () const;

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const;

    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector<albertCtype, cdim> & operator[] (int i) const;

    //! maps a local coordinate within reference element to
    //! global coordinate in element
    FieldVector<albertCtype, cdim> global (const FieldVector<albertCtype, mydim>& local) const;

    //! maps a global coordinate within the element to a
    //! local coordinate in its reference element
    FieldVector<albertCtype, mydim> local (const FieldVector<albertCtype, cdim>& global) const;

    //! returns true if the point in local coordinates is inside reference element
    bool checkInside(const FieldVector<albertCtype, mydim>& local) const;

    /*!
       Copy from sgrid.hh:

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

       Each grid module implements the integration element with optimal efficieny. This
       will directly translate in substantial savings in the computation of finite element
       stiffness matrices.
     */

    // A(l)
    albertCtype integrationElement (const FieldVector<albertCtype, mydim>& local) const;

    // volume if geometry
    albertCtype volume () const;

    //! can only be called for dim=dimworld!
    //! Note that if both methods are called on the same element, then
    //! call jacobianInverseTransposed first because integration element is calculated
    //! during calculation of the transposed of the jacobianInverse
    const FieldMatrix<albertCtype,cdim,mydim>& jacobianInverseTransposed (const FieldVector<albertCtype, mydim>& local) const;

    //***********************************************************************
    //!  Methods that not belong to the Interface, but have to be public
    //***********************************************************************
    //! generate the geometry for the ALBERTA EL_INFO
    //! no interface method
    typedef GridImp GridType;
    bool builtGeom(const GridImp & grid, ALBERTA EL_INFO *elInfo, int face, int edge, int vertex);

    //! build geometry for intersectionSelfLocal and
    //! intersectionNeighborLocal
    template <class GeometryType, class LocalGeomType >
    bool builtLocalGeom(const GeometryType & geo , const LocalGeomType & lg,
                        ALBERTA EL_INFO *elInfo, int face);

    // init geometry with zeros
    //! no interface method
    void initGeom();
    FieldVector<albertCtype, cdim>& getCoordVec (int i);

    //! print internal data
    //! no interface method
    void print (std::ostream& ss) const;

  private:
    // build geometry with local coords of child in reference element
    void buildGeomInFather(const int child, const int orientation);

    // calculate Matrix for Mapping from reference element to actual element
    void calcElMatrix () const;

    //! build the transposed of the jacobian inverse and store the volume
    void buildJacobianInverseTransposed () const;

    // template method for map the vertices of EL_INFO to the actual
    // coords with face_,edge_ and vertex_ , needes for operator []
    int mapVertices (int i) const;

    // calculates the volume of the element
    albertCtype elDeterminant () const;

    // temporary need vector
    mutable FieldVector<albertCtype, mydim+1> tmpVec_;

    //! the vertex coordinates
    mutable FieldMatrix<albertCtype,mydim+1,cdim> coord_;

    //! storage for global coords
    mutable FieldVector<albertCtype, cdim> globalCoord_;

    //! storage for local coords
    mutable FieldVector<albertCtype, mydim> localCoord_;

    // make empty EL_INFO
    ALBERTA EL_INFO * makeEmptyElInfo();

    ALBERTA EL_INFO * elInfo_;

    //! Which Face of the Geometry 0...dim+1
    int face_;

    //! Which Edge of the Face of the Geometry 0...dim
    int edge_;

    //! Which Edge of the Face of the Geometry 0...dim-1
    int vertex_;

    enum { matdim = (mydim > 0) ? mydim : 1 };
    mutable FieldMatrix<albertCtype,cdim,matdim> Jinv_; //!< storage for inverse of jacobian
    mutable FieldMatrix<albertCtype,matdim,matdim> Mtmp_;    //!< storage for inverse of jacobian

    mutable FieldMatrix<albertCtype,cdim,mydim> elMat_; //!< storage for mapping matrix
    mutable FieldMatrix<albertCtype,matdim,matdim> elMatT_elMat_; //!< storage for mapping matrix

    //! is true if elMat_ was calced
    mutable bool builtElMat_;
    //! is true if Jinv_ and volume_ is calced
    mutable bool builtinverse_;


    mutable bool calcedDet_; //! true if determinant was calculated
    mutable albertCtype elDet_; //!< storage of element determinant

    // temporary mem for integrationElement with mydim < cdim
    mutable FieldVector<albertCtype,cdim> tmpV_;
    mutable FieldVector<albertCtype,cdim> tmpU_;
    mutable FieldVector<albertCtype,cdim> tmpZ_;

    mutable FieldVector<albertCtype,mydim> AT_x_;
    const GeometryType myGeomType_;
    const albertCtype detFactor_;
  };

}

#endif
