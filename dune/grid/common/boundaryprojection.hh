// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BOUNDARYPROJECTION_HH
#define DUNE_BOUNDARYPROJECTION_HH

//- system includes
#include <cmath>

//- Dune includes
#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/genericgeometry/mappingprovider.hh>
#include <dune/grid/genericgeometry/geometrytraits.hh>

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

  template< int dim, int dimworld >
  class BoundarySegmentWrapper
    : public DuneBoundaryProjection< dimworld >
  {
    typedef BoundarySegmentWrapper< dim, dimworld > This;
    typedef DuneBoundaryProjection< dimworld > Base;

    typedef GenericGeometry::DefaultGeometryTraits< double, dim-1, dimworld > GeometryTraits;
    typedef GenericGeometry::HybridMapping< dim-1, GeometryTraits > FaceMapping;
    typedef GenericGeometry::MappingProvider< FaceMapping, 0 > FaceMappingProvider;

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
                             const shared_ptr< BoundarySegment > boundarySegment )
      : allocator_(),
        faceMapping_( FaceMappingProvider::mapping( GenericGeometry::topologyId( type ), vertices, allocator_ ) ),
        boundarySegment_( boundarySegment )
    {
      faceMapping_->referenceCount = 1;
    }

    BoundarySegmentWrapper ( const This &other )
      : allocator_( other.allocator_ ),
        faceMapping_( other.faceMapping_ ),
        boundarySegment_( other.boundarySegment_ )
    {
      ++(faceMapping_->referenceCount);
    }

    ~BoundarySegmentWrapper ()
    {
      if( --(faceMapping_->referenceCount) == 0 )
        allocator_.destroy( faceMapping_ );
    }

    This &operator= ( const This &other ) const
    {
      ++(other.faceMapping_->referenceCount);
      if( --(faceMapping_->referenceCount == 0) )
        allocator_.destroy( faceMapping_ );
      allocator_ = other.allocator_;
      faceMapping_ = other.faceMapping_;
      boundarySegment_ = other.boundarySegment_;
      return *this;
    }

    CoordinateType operator() ( const CoordinateType &global ) const
    {
      return boundarySegment() ( faceMapping_->local( global ) );
    }

    const BoundarySegment &boundarySegment () const
    {
      return *boundarySegment_;
    }

  private:
    typename GeometryTraits::Allocator allocator_;
    FaceMapping *faceMapping_;
    const shared_ptr< BoundarySegment > boundarySegment_;
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
