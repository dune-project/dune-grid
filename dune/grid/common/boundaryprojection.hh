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
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/io/file/gmshreader.hh>

namespace Dune
{
  /** \brief Interface class for vertex projection at the boundary.
   */
  template <int dimworld>
  struct DuneBoundaryProjection;

  /** \brief Interface class for vertex projection at the boundary.
   */
  template <int dimworld>
  struct DuneBoundaryProjection
    : public BoundarySegmentBackupRestore< DuneBoundaryProjection< dimworld > >
  {
    typedef DuneBoundaryProjection< dimworld > ThisType;
    typedef BoundarySegmentBackupRestore< DuneBoundaryProjection< dimworld > > BaseType;

    using BaseType :: restore;
    using BaseType :: registerFactory;

    //! \brief type of coordinate vector
    typedef FieldVector< double, dimworld> CoordinateType;
    //! \brief destructor
    virtual ~DuneBoundaryProjection() {}

    //! \brief projection operator projection a global coordinate
    virtual CoordinateType operator() (const CoordinateType& global) const = 0;

    /** \brief write DuneBoundaryProjection's data to stream buffer
     *  \param buffer buffer to store data
     */
    virtual void backup( std::stringstream& buffer ) const
    {
      DUNE_THROW(NotImplemented,"DuneBoundaryProjection::backup not overloaded!");
    }

    template <class BufferImp>
    void toBuffer( BufferImp& buffer ) const
    {
      MessageBufferIF< BufferImp > buf( buffer );
      toBuffer( buf );
    }

    template <class BufferImp>
    void toBuffer( MessageBufferIF< BufferImp > & buffer ) const
    {
      std::stringstream str;
      // call virtual interface backup
      backup( str );
      std::string data = str.str();
      const size_t size = data.size();
      buffer.write( size );
      for( size_t i=0; i<size; ++i )
        buffer.write( data[ i ] );
    }

    template <class BufferImp>
    static ThisType* restoreFromBuffer( BufferImp & buffer )
    {
      MessageBufferIF< BufferImp > buf( buffer );
      return restoreFromBuffer( buf );
    }

    template <class BufferImp>
    static ThisType* restoreFromBuffer( MessageBufferIF< BufferImp > & buffer )
    {
      std::string data;
      size_t size = 0;
      buffer.read( size );
      data.resize( size );
      for( size_t i=0; i<size; ++i )
        buffer.read( data[ i ] );

      std::stringstream str;
      str.write( data.c_str(), size );
      return BaseType::restore( str );
    }
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

    typedef GmshReaderQuadraticBoundarySegment< dim, dimworld >  GmshBoundarySegment;
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

    void backup( std::stringstream& buffer ) const
    {
      // write identifier key first
      buffer.write( key(), Base::keyLength );
      // now all data
      GeometryType type = faceMapping_.type();
      buffer.write( (const char *) &type, sizeof(GeometryType) );

      int corners = faceMapping_.corners() ;
      buffer.write( (const char *) &corners, sizeof(int) );

      CoordinateType corner( 0 );
      for( int i=0; i<corners; ++i )
      {
        corner = faceMapping_.corner( i );
        buffer.write( (const char *) &corner[ 0 ], sizeof(double)*CoordinateType::dimension );
      }

      boundarySegment_->backup( buffer );
    }

    static void registerFactory()
    {
      Base::registerFactory( key(), &factory );
    }

  protected:
    static const char* key () { return "bswr"; }

    // create and object of this class from a stream buffer
    static Base* factory( std::stringstream& buffer )
    {
      GeometryType type;
      buffer.read( (char *) &type, sizeof(GeometryType) );
      int corners = 0;
      buffer.read( (char *) &corners, sizeof(int) );
      std::vector< CoordinateType > vertices( corners, CoordinateType(0) );
      for( int i=0; i<corners; ++i )
      {
        buffer.read( (char *) &vertices[ i ][ 0 ], sizeof(double)*CoordinateType::dimension );
      }
      std::shared_ptr< BoundarySegment > ptr( BoundarySegment::restore( buffer ) );
      return new BoundarySegmentWrapper< dim, dimworld >( type, vertices, ptr );
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
