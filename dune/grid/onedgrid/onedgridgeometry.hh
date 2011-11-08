// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GEOMETRY_HH
#define DUNE_ONE_D_GEOMETRY_HH

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

    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector<typename GridImp::ctype, 1>& operator[] (int i) const {
      return (storeCoordsLocally_) ? pos_ : target_->pos_;
    }

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
    public GeometryDefaultImplementation <mydim, coorddim, GridImp, OneDGridGeometry>
  {
    template <int codim_, int dim_, class GridImp_>
    friend class OneDGridEntity;

    friend class OneDGrid;

    template <int cc_, int dim_, class GridImp_>
    friend class OneDGridSubEntityFactory;

    template <class GridImp_>
    friend class OneDGridLevelIntersectionIterator;
    template <class GridImp_>
    friend class OneDGridLeafIntersectionIterator;

  public:

    OneDGridGeometry() : storeCoordsLocally_(false) {}

    //! here we have always an affine geometry
    bool affine() const { return true; }

    /** \brief Return the element type identifier
     *
     * OneDGrid obviously supports only lines
     */
    GeometryType type () const {return GeometryType(1);}

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const {return 2;}

    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector<typename GridImp::ctype, coorddim>& operator[](int i) const {
      assert(i==0 || i==1);
      return (storeCoordsLocally_) ? pos_[i] : target_->vertex_[i]->pos_;
    }

    /** \brief access to coordinates of a corner */
    const FieldVector< typename GridImp::ctype, coorddim > corner ( const int i ) const
    {
      assert( (i == 0) || (i == 1) );
      return (storeCoordsLocally_ ? pos_[ i ] : target_->vertex_[ i ]->pos_);
    }

    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    FieldVector<typename GridImp::ctype, coorddim> global (const FieldVector<typename GridImp::ctype, mydim>& local) const {
      FieldVector<typename GridImp::ctype, coorddim> g;
      g[0] = (storeCoordsLocally_)
             ? pos_[0][0] * (1-local[0]) + pos_[1][0] * local[0]
             : target_->vertex_[0]->pos_[0] * (1-local[0]) + target_->vertex_[1]->pos_[0] * local[0];
      return g;
    }

    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    FieldVector<typename GridImp::ctype, mydim> local (const FieldVector<typename GridImp::ctype, coorddim>& global) const {
      FieldVector<typename GridImp::ctype, mydim> l;
      if (storeCoordsLocally_) {
        l[0] = (global[0] - pos_[0][0]) / (pos_[1][0] - pos_[0][0]);
      } else {
        const typename GridImp::ctype& v0 = target_->vertex_[0]->pos_[0];
        const typename GridImp::ctype& v1 = target_->vertex_[1]->pos_[0];
        l[0] = (global[0] - v0) / (v1 - v0);
      }
      return l;
    }

    /** ???
     */
    typename GridImp::ctype integrationElement (const FieldVector<typename GridImp::ctype, mydim>& local) const {
      return (storeCoordsLocally_)
             ? pos_[1][0] - pos_[0][0]
             : target_->vertex_[1]->pos_[0] - target_->vertex_[0]->pos_[0];
    }

    //! The Jacobian matrix of the mapping from the reference element to this element
    const FieldMatrix< typename GridImp::ctype, mydim, mydim > &
    jacobianTransposed ( const FieldVector< typename GridImp::ctype, mydim > &local ) const
    {
      if( storeCoordsLocally_ )
        jacTransposed_[ 0 ][ 0 ] = (pos_[ 1 ][ 0 ] - pos_[ 0 ][ 0 ] );
      else
        jacTransposed_[ 0 ][ 0 ] = target_->vertex_[ 1 ]->pos_[ 0 ] - target_->vertex_[ 0 ]->pos_[ 0 ];

      return jacTransposed_;
    }

    //! The Jacobian matrix of the mapping from the reference element to this element
    const FieldMatrix<typename GridImp::ctype,mydim,mydim>& jacobianInverseTransposed (const FieldVector<typename GridImp::ctype, mydim>& local) const {
      if (storeCoordsLocally_)
        jacInverse_[0][0] = 1 / (pos_[1][0] - pos_[0][0]);
      else
        jacInverse_[0][0] = 1 / (target_->vertex_[1]->pos_[0] - target_->vertex_[0]->pos_[0]);

      return jacInverse_;
    }

    void setPositions(const typename GridImp::ctype& p1, const typename GridImp::ctype& p2) {
      storeCoordsLocally_ = true;
      pos_[0][0] = p1;
      pos_[1][0] = p2;
    }

    //private:
    OneDEntityImp<1>* target_;

    bool storeCoordsLocally_;

    // Stores the element corner positions if it is returned as geometryInFather
    FieldVector<typename GridImp::ctype,coorddim> pos_[2];

    //! jacobian transposed
    mutable FieldMatrix< typename GridImp::ctype, coorddim, coorddim > jacTransposed_;
    //! jacobian inverse
    mutable FieldMatrix<typename GridImp::ctype,coorddim,coorddim> jacInverse_;

  };

}  // namespace Dune

#endif