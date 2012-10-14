// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GEOMETRY_HH
#define DUNE_ONE_D_GEOMETRY_HH

#include <dune/geometry/axisalignedcubegeometry.hh>

#include <dune/grid/common/geometry.hh>

#include <dune/grid/onedgrid/onedgridentity.hh>

/** \file
 * \brief The OneDGridElement class and its specializations
 */

namespace Dune {

  // forward declaration
  template <int codim, int dim, class GridImp>
  class OneDGridEntity;

  /** \brief Unspecialized class.  Not used for anything */
  template<int mydim, int coorddim, class GridImp>
  class OneDGridGeometry;

  template<class GridImp>
  class OneDGridGeometry <0, 1, GridImp> :
    public GeometryDefaultImplementation <0, 1, GridImp,OneDGridGeometry>
  {

    template <int codim_, int dim_, class GridImp_>
    friend class OneDGridEntity;
    template <int mydim_, int coorddim_, class GridImp_>
    friend class OneDGridGeometry;

  public:

    OneDGridGeometry() : storeCoordsLocally_(false) {}

    //! return the element type identifier (vertex)
    GeometryType type () const {return GeometryType(0);}

    //! here we have always an affine geometry
    bool affine() const { return true; }

    //! return the number of corners of this element (==1)
    int corners () const {return 1;}

    /** \brief access to coordinates of a corner */
    const FieldVector< typename GridImp::ctype, 1 > corner ( const int i ) const
    {
      return (storeCoordsLocally_ ? pos_ : target_->pos_);
    }

    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    FieldVector<typename GridImp::ctype, 1> global (const FieldVector<typename GridImp::ctype, 0>& local) const {
      return (storeCoordsLocally_) ? pos_ : target_->pos_;
    }

    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    FieldVector<typename GridImp::ctype, 0> local (const FieldVector<typename GridImp::ctype, 1>& global) const {
      FieldVector<typename GridImp::ctype, 0> l;
      return l;
    }

    /** \brief !!!

       This method really doesn't make much sense for a zero-dimensional
       object.  It always returns '1'.
     */
    typename GridImp::ctype integrationElement (const FieldVector<typename GridImp::ctype, 0>& local) const {
      return 1;
    }

    //! The Jacobian matrix of the mapping from the reference element to this element
    const FieldMatrix< typename GridImp::ctype, 0, 1 > &
    jacobianTransposed ( const FieldVector< typename GridImp::ctype, 0 > &local ) const
    {
      return jacTransposed_;
    }

    //! The Jacobian matrix of the mapping from the reference element to this element
    const FieldMatrix<typename GridImp::ctype,1,0>& jacobianInverseTransposed (const FieldVector<typename GridImp::ctype, 0>& local) const {
      return jacInverse_;
    }

    void setPosition(const typename GridImp::ctype& p) {
      storeCoordsLocally_ = true;
      pos_[0] = p;
    }

    //private:
    bool storeCoordsLocally_;

    // Stores the element corner positions if it is returned as geometryInFather
    FieldVector<typename GridImp::ctype,1> pos_;

    OneDEntityImp<0>* target_;

    FieldMatrix< typename GridImp::ctype, 0, 1 > jacTransposed_;
    FieldMatrix<typename GridImp::ctype,1,0> jacInverse_;
  };

  //**********************************************************************
  //
  // --OneDGridGeometry
  /** \brief Defines the geometry part of a mesh entity.
   * \ingroup OneDGrid
   */
  template<int mydim, int coorddim, class GridImp>
  class OneDGridGeometry :
    public AxisAlignedCubeGeometry<typename GridImp::ctype, mydim, coorddim>
  {
  public:
    OneDGridGeometry(const FieldVector<double,1>& left, const FieldVector<double,1>& right)
      : AxisAlignedCubeGeometry<typename GridImp::ctype, mydim, coorddim>(left,right)
    {}
  };

  namespace FacadeOptions
  {

    /** \brief Switch on the new implementation for the Geometry interface class
     * \deprecated Eventually the new implementation will be hardwired,
     *             and this switch may disappear without prior warning!
     */
    template< int mydim, int cdim, class GridImp>
    struct StoreGeometryReference<mydim,cdim,GridImp,OneDGridGeometry>
    {
      static const bool v = false;
    };

  }

}  // namespace Dune

#endif
