// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BOUNDARYPROJECTION_HH
#define DUNE_BOUNDARYPROJECTION_HH

//- system includes
#include <cmath>
#include <memory>

//- Dune includes
#include <dune/common/fvector.hh>

#include <dune/geometry/multilineargeometry.hh>

#include <dune/grid/common/boundarysegment.hh>

namespace Dune
{

  /** \brief Interface class for vertex projection at the boundary.
   */
  template <int dimworld>
  struct DuneBoundaryProjection
  {
    //! \brief type of coordinate vector
    typedef FieldVector< double, dimworld> CoordinateType;
    //! \brief destructor
    virtual ~DuneBoundaryProjection() {}

    //! \brief projection operator projection a global coordinate
    virtual CoordinateType operator() (const CoordinateType& global) const = 0;
  };

  template < int dimworld >
  class BoundaryProjectionWrapper
    : public DuneBoundaryProjection< dimworld >
  {
  protected:
    typedef DuneBoundaryProjection< dimworld > BaseType;
    const BaseType& proj_;
  public:
    //! \brief type of coordinate vector
    typedef typename BaseType :: CoordinateType CoordinateType;

    // constructor taking other projection
    BoundaryProjectionWrapper( const BaseType& proje )
      : proj_( proje )
    {}

    //! destructor
    ~BoundaryProjectionWrapper () {}

    //! \brief projection operator projection a global coordinate
    CoordinateType operator() (const CoordinateType& global) const
    {
      return proj_( global );
    }
  };

  // BoundarySegmentWrapper
  // ----------------------

  /** \tparam dim Dimension of the grid */
  template< int dim, int dimworld >
  class BoundarySegmentWrapper
    : public DuneBoundaryProjection< dimworld >
  {
    typedef DuneBoundaryProjection< dimworld > Base;

    typedef MultiLinearGeometry<typename Base::CoordinateType::value_type,dim-1,dimworld> FaceMapping;

  public:
    typedef typename Base::CoordinateType CoordinateType;
    typedef Dune::BoundarySegment< dim, dimworld > BoundarySegment;

    /** constructor
     *
     *  \param[in]  type             geometry type of the boundary face
     *  \param[in]  vertices         vertices of the boundary face
     *  \param[in]  boundarySegment  geometric realization of the shaped boundary
     *
     *  \note The BoundarySegmentWrapper takes control of the boundary segment.
     */
    BoundarySegmentWrapper ( const GeometryType &type,
                             const std::vector< CoordinateType > &vertices,
                             const std::shared_ptr< BoundarySegment > &boundarySegment )
      : faceMapping_( FaceMapping( type, vertices ) ),
        boundarySegment_( boundarySegment )
    {}

    CoordinateType operator() ( const CoordinateType &global ) const
    {
      return boundarySegment() ( faceMapping_.local( global ) );
    }

    const BoundarySegment &boundarySegment () const
    {
      return *boundarySegment_;
    }

  private:
    FaceMapping faceMapping_;
    const std::shared_ptr< BoundarySegment > boundarySegment_;
  };



  //////////////////////////////////////////////////////////////////////
  //
  // Example of boundary projection projection to a circle
  //
  //////////////////////////////////////////////////////////////////////
  template <int dimworld>
  struct CircleBoundaryProjection : public DuneBoundaryProjection< dimworld >
  {
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

  protected:
    //! radius of circ
    const double radius_;
  };

} // end namespace

#endif // #ifndef DUNE_BOUNDARYPROJECTION_HH
