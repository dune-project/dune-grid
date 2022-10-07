// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_GRIDPTR_HH
#define DUNE_DGF_GRIDPTR_HH

#include <cassert>
#include <cctype>

#include <array>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

//- Dune includes
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/dgfparser/entitykey.hh>
#include <dune/grid/io/file/dgfparser/parser.hh>

#include <dune/grid/io/file/gmshreader.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template < class G >
  struct DGFGridFactory;

  template< class GridImp, class IntersectionImp >
  class Intersection;



  // GridPtr
  // -------

  //! \brief Class for constructing grids from DGF files.
  //!
  //! The constructor of the class is given the filename of the DGF file.
  //! From that file a pointer to an instance of type GridType is created by reading
  //! the given file which is translated to the specific format of the given
  //! GridType. The GridPtr class behaves like an auto pointer of GridType.
  //! An auto pointer to a grid of type GridType is constructed
  //! as follows:
  //! @code
  //! GridPtr<GridType> gridptr(filename, MPI_COMM_WORLD );
  //! GridType & grid = *gridptr;
  //! @endcode
  template< class GridType >
  struct GridPtr
  {
    class mygrid_ptr : public std::shared_ptr< GridType >
    {
      typedef std::shared_ptr< GridType > base_t ;
      // empty deleter to avoid deletion on release
      typedef null_deleter< GridType > emptydeleter_t ;

      void removeObj()
      {
        // if use count is only 1 delete object
        if( use_count() == 1 )
        {
          // delete point here, since we use the empty deleter
          GridType* grd = release();
          if( grd ) delete grd ;
        }
      }

      void assignObj( const mygrid_ptr& other )
      {
        removeObj();
        base_t :: operator = ( other );
      }
    public:
      using base_t :: get ;
      using base_t :: swap ;
      using base_t :: use_count  ;

      // default constructor
      mygrid_ptr() : base_t( ( GridType * ) 0, emptydeleter_t() ) {}
      // copy constructor
      mygrid_ptr( const mygrid_ptr& other ) : base_t(nullptr) { assignObj( other ); }
      // constructor taking pointer
      explicit mygrid_ptr( GridType* grd ) : base_t( grd, emptydeleter_t() ) {}

      // destructor
      ~mygrid_ptr() { removeObj(); }

      // assigment operator
      mygrid_ptr& operator = ( const mygrid_ptr& other )
      {
        assignObj( other );
        return *this;
      }

      // release pointer
      GridType* release()
      {
        GridType* grd = this->get();
        base_t ptr(( GridType * ) 0, emptydeleter_t() );
        this->swap( ptr );
        return grd ;
      }
    };

  protected:
    std::string getFileExtension( const std::string& filename ) const
    {
      // extract file extension
      auto extpos = filename.find_last_of(".");
      std::string ext;
      if( extpos != std::string::npos)
        ext = filename.substr( extpos + 1 );

      // convert all letters to lower case
      for( auto& item : ext )
        item = std::tolower( item );
      return ext;
    }

    // read gmsh file if dimension world <= 3
    void readGmsh( const std::string& filename, std::integral_constant< bool, true > )
    {
      GridFactory<GridType> gridFactory;
      std::vector<int> boundaryIDs;
      std::vector<int> elementsIDs;
      GmshReader<GridType>::read(gridFactory,filename,boundaryIDs,elementsIDs);
      initialize( gridFactory, boundaryIDs,elementsIDs);
    }

    // if dimension world > 3 throw GridError
    void readGmsh( const std::string& filename, std::integral_constant< bool, false > )
    {
      DUNE_THROW(GridError, "GmshReader requires dimWorld <= 3." );
    }

  public:

    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    static const int dimension = GridType::dimension;

    //! constructor given the name of a DGF file
    explicit GridPtr ( const std::string &filename,
                       MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : gridPtr_(),
        elParam_(),
        vtxParam_(),
        bndParam_(),
        bndId_(),
        emptyParam_(),
        nofElParam_( 0 ),
        nofVtxParam_( 0 ),
        haveBndParam_( false )
    {
      std::string fileExt = getFileExtension( filename );

      if( fileExt == "dgf" )
      {
        DGFGridFactory< GridType > dgfFactory( filename, comm );
        initialize( dgfFactory );
      }
      else if( fileExt == "msh" )
      {
        // Gmsh reader only compiles for dimworld <= 3
        readGmsh( filename, std::integral_constant< bool, GridType::dimensionworld <= 3 > () );
      }
      else if( fileExt == "amc" || fileExt == "2d" || fileExt == "3d" )
      {
        // TODO: AlbertaReader
        DUNE_THROW( NotImplemented, "GridPtr: file format '" << fileExt << "' not supported yet!" );
      }
      else if( fileExt == "vtu" )
      {
        // TODO: vtu/vtk reader
        DUNE_THROW( NotImplemented, "GridPtr: file format '" << fileExt << "' not supported yet!" );
      }
      else
      {
        DUNE_THROW( NotImplemented, "GridPtr: file format '" << fileExt << "' not supported yet!" );
      }
    }

    //! constructor given a std::istream
    explicit GridPtr ( std::istream &input,
                       MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : gridPtr_(),
        elParam_(),
        vtxParam_(),
        bndParam_(),
        bndId_(),
        emptyParam_(),
        nofElParam_( 0 ),
        nofVtxParam_( 0 ),
        haveBndParam_( false )
    {
      // input stream only works for DGF format right now
      DGFGridFactory< GridType > dgfFactory( input, comm );
      initialize( dgfFactory );
    }

    //! Default constructor, creating empty GridPtr
    GridPtr()
      : gridPtr_(),
        elParam_(),
        vtxParam_(),
        bndParam_(),
        bndId_(),
        emptyParam_(),
        nofElParam_(0),
        nofVtxParam_(0),
        haveBndParam_( false )
    {}

    //! Constructor storing given pointer to internal auto pointer
    explicit GridPtr( GridType *grd )
      : gridPtr_(grd),
        elParam_(),
        vtxParam_(),
        bndParam_(),
        bndId_(),
        emptyParam_(),
        nofElParam_(0),
        nofVtxParam_(0),
        haveBndParam_( false )
    {}

    //! Copy constructor, copies internal auto pointer
    GridPtr( const GridPtr &org ) = default;

    //! assignment of grid pointer
    GridPtr& operator= ( const GridPtr &org )
    {
      gridPtr_    = org.gridPtr_;
      elParam_    = org.elParam_;
      vtxParam_   = org.vtxParam_;
      bndParam_   = org.bndParam_;
      bndId_      = org.bndId_;
      emptyParam_ = org.emptyParam_;

      nofElParam_ = org.nofElParam_;
      nofVtxParam_ = org.nofVtxParam_;
      haveBndParam_ = org.haveBndParam_;
      return *this;
    }

    //! assignment of pointer to internal auto pointer
    GridPtr& operator = (GridType * grd)
    {
      gridPtr_ = mygrid_ptr( grd );
      elParam_.resize(0);
      vtxParam_.resize(0);
      bndParam_.resize(0);
      bndId_.resize(0);
      emptyParam_.resize(0);

      nofVtxParam_ = 0;
      nofElParam_ = 0;
      haveBndParam_ = false;
      return *this;
    }

    //! return reference to GridType instance
    GridType& operator*() {
      return *gridPtr_;
    }

    //! return pointer to GridType instance
    GridType* operator->() {
      return gridPtr_.operator -> ();
    }

    //! return const reference to GridType instance
    const GridType& operator*() const {
      return *gridPtr_;
    }

    //! return const pointer to GridType instance
    const GridType* operator->() const {
      return gridPtr_.operator -> ();
    }

    //! release pointer from internal ownership
    GridType* release () { return gridPtr_.release();  }

    //! get number of parameters defined for a given codimension
    int nofParameters(int cdim) const {
      switch (cdim) {
      case 0 : return nofElParam_; break;
      case GridType::dimension : return nofVtxParam_; break;
      }
      return 0;
    }

    //! get parameters defined for given entity
    template <class Entity>
    int nofParameters ( const Entity & ) const
    {
      return nofParameters( (int) Entity::codimension );
    }

    //! get number of parameters defined for a given intersection
    template< class GridImp, class IntersectionImp >
    int nofParameters ( const Intersection< GridImp, IntersectionImp > & intersection ) const
    {
      return parameters( intersection ).size();
    }

    //! get parameters defined for each codim 0 und dim entity on the grid through the grid file
    template <class Entity>
    const std::vector< double > &parameters ( const Entity &entity ) const
    {
      typedef typename GridType::LevelGridView GridView;
      GridView gridView = gridPtr_->levelGridView( 0 );
      switch( (int)Entity::codimension )
      {
      case 0 :
        if( nofElParam_ > 0 )
        {
          assert( (unsigned int)gridView.indexSet().index( entity ) < elParam_.size() );
          return elParam_[ gridView.indexSet().index( entity ) ];
        }
        break;
      case GridType::dimension :
        if( nofVtxParam_ > 0 )
        {
          assert( (unsigned int)gridView.indexSet().index( entity ) < vtxParam_.size() );
          return vtxParam_[ gridView.indexSet().index( entity ) ];
        }
        break;
      }
      return emptyParam_;
    }

    //! get parameters for intersection
    template< class GridImp, class IntersectionImp >
    const DGFBoundaryParameter::type & parameters ( const Intersection< GridImp, IntersectionImp > & intersection ) const
    {
      // if no parameters given return empty vector
      if ( !haveBndParam_ )
        return DGFBoundaryParameter::defaultValue();

      return bndParam_[ intersection.boundarySegmentIndex() ];
    }

    void communicate ()
    {
      if( gridPtr_->comm().size() > 1 )
      {
        DataHandle dh(*this);
        gridPtr_->levelGridView( 0 ).communicate( dh.interface(), InteriorBorder_All_Interface,ForwardCommunication );
      }
    }

    void loadBalance ()
    {
      if( gridPtr_->comm().size() > 1 )
      {
        DataHandle dh(*this);
        gridPtr_->loadBalance( dh.interface() );
        gridPtr_->levelGridView( 0 ).communicate( dh.interface(), InteriorBorder_All_Interface,ForwardCommunication );
      }
    }

  protected:
    template< class Range >
    static bool isEmpty ( Range &&range )
    {
      return range.begin() == range.end();
    }

    void initialize ( DGFGridFactory< GridType > &dgfFactory )
    {
      gridPtr_ = mygrid_ptr( dgfFactory.grid() );

      const auto gridView = gridPtr_->levelGridView( 0 );
      const auto &indexSet = gridView.indexSet();

      nofElParam_ = dgfFactory.template numParameters< 0 >();
      nofVtxParam_ = dgfFactory.template numParameters< dimension >();
      haveBndParam_ = dgfFactory.haveBoundaryParameters();

      std::array< int, 3 > nofParams = {{ nofElParam_, nofVtxParam_, static_cast< int >( haveBndParam_ ) }};
      gridView.comm().max( nofParams.data(), nofParams.size() );

      // empty grids have no parameters associated
      if( isEmpty( elements( gridView, Partitions::interiorBorder ) ) )
      {
        nofElParam_ = nofParams[ 0 ];
        nofVtxParam_ = nofParams[ 1 ];
      }

      // boundary parameters may be empty
      haveBndParam_ = static_cast< bool >( nofParams[ 2 ] );

      if( (nofElParam_ != nofParams[ 0 ]) || (nofVtxParam_ != nofParams[ 1 ]) )
        DUNE_THROW( DGFException, "Number of parameters differs between processes" );

      elParam_.resize( nofElParam_ > 0 ? indexSet.size( 0 ) : 0 );
      vtxParam_.resize( nofVtxParam_ > 0 ? indexSet.size( dimension ) : 0 );

      bndId_.resize( indexSet.size( 1 ) );
      if( haveBndParam_ )
        bndParam_.resize( gridPtr_->numBoundarySegments() );

      for( const auto &element : elements( gridView, Partitions::interiorBorder ) )
      {
        if( nofElParam_ > 0 )
        {
          std::swap( elParam_[ indexSet.index( element ) ], dgfFactory.parameter( element ) );
          assert( elParam_[ indexSet.index( element ) ].size() == static_cast< std::size_t >( nofElParam_ ) );
        }

        if( nofVtxParam_ > 0 )
        {
          for( unsigned int v = 0, n = element.subEntities( dimension ); v < n; ++v )
          {
            const auto index = indexSet.subIndex( element, v, dimension );
            if( vtxParam_[ index ].empty() )
              std::swap( vtxParam_[ index ], dgfFactory.parameter( element.template subEntity< dimension >( v ) ) );
            assert( vtxParam_[ index ].size() == static_cast< std::size_t >( nofVtxParam_ ) );
          }
        }

        if( element.hasBoundaryIntersections() )
        {
          for( const auto &intersection : intersections( gridView, element ) )
          {
            // dirty hack: check for "none" to make corner point grid work
            if( !intersection.boundary() || intersection.type().isNone() )
              continue;

            const auto k = indexSet.subIndex( element, intersection.indexInInside(), 1 );
            bndId_[ k ] = dgfFactory.boundaryId( intersection );
            if( haveBndParam_ )
              bndParam_[ intersection.boundarySegmentIndex() ] = dgfFactory.boundaryParameter( intersection );
          }
        }
      }
    }

    void initialize ( GridFactory< GridType > &factory,
                      std::vector<int>& boundaryIds,
                      std::vector<int>& elementIds )
    {
      gridPtr_ = mygrid_ptr( factory.createGrid().release() );

      const auto& gridView = gridPtr_->leafGridView();
      const auto& indexSet = gridView.indexSet();

      nofElParam_   = elementIds.empty() ? 0 : 1 ;
      nofVtxParam_  = 0;
      haveBndParam_ = boundaryIds.empty() ? 0 : 1 ;

      std::array< int, 3 > nofParams = {{ nofElParam_, nofVtxParam_, static_cast< int >( haveBndParam_ ) }};
      gridView.comm().max( nofParams.data(), nofParams.size() );

      // empty grids have no parameters associated
      if( isEmpty( elements( gridView, Partitions::interiorBorder ) ) )
      {
        nofElParam_ = nofParams[ 0 ];
      }

      // boundary parameters may be empty
      haveBndParam_ = static_cast< bool >( nofParams[ 2 ] );

      // Reorder boundary IDs according to the insertion index
      if(!boundaryIds.empty() || !elementIds.empty() )
      {
        bndParam_.resize( boundaryIds.size() );
        elParam_.resize( elementIds.size(), std::vector<double>(1) );
        for(const auto& entity : elements( gridView ))
        {
          elParam_[ indexSet.index( entity ) ][ 0 ] = elementIds[ factory.insertionIndex( entity ) ];
          if( haveBndParam_ )
          {
            for(const auto& intersection : intersections( gridView,entity) )
            {
              if(intersection.boundary())
              {
                // DGFBoundaryParameter::type is of type string.
                bndParam_[intersection.boundarySegmentIndex()] = std::to_string(boundaryIds[factory.insertionIndex(intersection)]);
              }
            }
          }
        }
      }
    }

    template <class Entity>
    std::vector< double > &params ( const Entity &entity )
    {
      const auto gridView = gridPtr_->levelGridView( 0 );
      switch( (int)Entity::codimension )
      {
      case 0 :
        if( nofElParam_ > 0 ) {
          if ( gridView.indexSet().index( entity ) >= elParam_.size() )
            elParam_.resize( gridView.indexSet().index( entity ) );
          return elParam_[ gridView.indexSet().index( entity ) ];
        }
        break;
      case GridType::dimension :
        if( nofVtxParam_ > 0 ) {
          if ( gridView.indexSet().index( entity ) >= vtxParam_.size() )
            vtxParam_.resize( gridView.indexSet().index( entity ) );
          return vtxParam_[ gridView.indexSet().index( entity ) ];
        }
        break;
      }
      return emptyParam_;
    }

    void setNofParams( int cdim, int nofP )
    {
      switch (cdim) {
      case 0 : nofElParam_ = nofP; break;
      case GridType::dimension : nofVtxParam_ = nofP; break;
      }
    }

    struct DataHandle
      : public CommDataHandleIF< DataHandle, char >
    {
      explicit DataHandle ( GridPtr &gridPtr )
        : gridPtr_( gridPtr ), idSet_( gridPtr->localIdSet() )
      {
        const auto gridView = gridPtr_->levelGridView( 0 );
        const auto &indexSet = gridView.indexSet();

        for( const auto &element : elements( gridView, Partitions::interiorBorder ) )
        {
          if( gridPtr_.nofElParam_ > 0 )
            std::swap( gridPtr_.elParam_[ indexSet.index( element ) ], elData_[ idSet_.id( element ) ] );

          if( gridPtr_.nofVtxParam_ > 0 )
          {
            for( unsigned int v = 0, n = element.subEntities( dimension ); v < n; ++v )
            {
              const auto index = indexSet.subIndex( element, v, dimension );
              if ( !gridPtr_.vtxParam_[ index ].empty() )
                std::swap( gridPtr_.vtxParam_[ index ], vtxData_[ idSet_.subId( element, v, dimension ) ] );
            }
          }

          if( element.hasBoundaryIntersections() )
          {
            for( const auto &intersection : intersections( gridView, element ) )
            {
              // dirty hack: check for "none" to make corner point grid work
              if( !intersection.boundary() || intersection.type().isNone() )
                continue;

              const int i = intersection.indexInInside();
              auto &bndData = bndData_[ idSet_.subId( element, i, 1 ) ];
              bndData.first = gridPtr_.bndId_[ indexSet.subIndex( element, i, 1 ) ];
              if( gridPtr_.haveBndParam_ )
                std::swap( bndData.second, gridPtr_.bndParam_[ intersection.boundarySegmentIndex() ] );
            }
          }
        }
      }

      DataHandle ( const DataHandle & ) = delete;
      DataHandle ( DataHandle && ) = delete;

      ~DataHandle ()
      {
        const auto gridView = gridPtr_->levelGridView( 0 );
        const auto &indexSet = gridView.indexSet();

        if( gridPtr_.nofElParam_ > 0 )
          gridPtr_.elParam_.resize( indexSet.size( 0 ) );
        if( gridPtr_.nofVtxParam_ > 0 )
          gridPtr_.vtxParam_.resize( indexSet.size( dimension ) );
        gridPtr_.bndId_.resize( indexSet.size( 1 ) );
        if( gridPtr_.haveBndParam_ )
            gridPtr_.bndParam_.resize( gridPtr_->numBoundarySegments() );

        for( const auto &element : elements( gridView, Partitions::all ) )
        {
          if( gridPtr_.nofElParam_ > 0 )
          {
            std::swap( gridPtr_.elParam_[ indexSet.index( element ) ], elData_[ idSet_.id( element ) ] );
            assert( gridPtr_.elParam_[ indexSet.index( element ) ].size() == static_cast< std::size_t >( gridPtr_.nofElParam_ ) );
          }

          if( gridPtr_.nofVtxParam_ > 0 )
          {
            for( unsigned int v = 0; v < element.subEntities( dimension ); ++v )
            {
              const auto index = indexSet.subIndex( element, v, dimension );
              if( gridPtr_.vtxParam_[ index ].empty() )
                std::swap( gridPtr_.vtxParam_[ index ], vtxData_[ idSet_.subId( element, v, dimension ) ] );
              assert( gridPtr_.vtxParam_[ index ].size() == static_cast< std::size_t >( gridPtr_.nofVtxParam_ ) );
            }
          }

          if( element.hasBoundaryIntersections() )
          {
            for( const auto &intersection : intersections( gridView, element ) )
            {
              // dirty hack: check for "none" to make corner point grid work
              if( !intersection.boundary() || intersection.type().isNone() )
                continue;

              const int i = intersection.indexInInside();
              auto &bndData = bndData_[ idSet_.subId( element, i, 1 ) ];
              gridPtr_.bndId_[ indexSet.subIndex( element, i, 1 ) ] = bndData.first;
              if( gridPtr_.haveBndParam_ )
                std::swap( bndData.second, gridPtr_.bndParam_[ intersection.boundarySegmentIndex() ] );
            }
          }
        }
      }

      CommDataHandleIF< DataHandle, char > &interface () { return *this; }

      bool contains ( int dim, int codim ) const
      {
        assert( dim == dimension );
        // do not use a switch statement, because dimension == 1 is possible
        return (codim == 1) || ((codim == dimension) && (gridPtr_.nofVtxParam_ > 0)) || ((codim == 0) && (gridPtr_.nofElParam_ > 0));
      }

      bool fixedSize (int /* dim */, int /* codim */) const { return false; }

      template< class Entity >
      std::size_t size ( const Entity &entity ) const
      {
        std::size_t totalSize = 0;

        // do not use a switch statement, because dimension == 1 is possible
        if( (Entity::codimension == 0) && (gridPtr_.nofElParam_ > 0) )
        {
          assert( elData_[ idSet_.id( entity ) ].size() == static_cast< std::size_t >( gridPtr_.nofElParam_ ) );
          for( double &v : elData_[ idSet_.id( entity ) ] )
            totalSize += dataSize( v );
        }

        if( (Entity::codimension == dimension) && (gridPtr_.nofVtxParam_ > 0) )
        {
          assert( vtxData_[ idSet_.id( entity ) ].size() == static_cast< std::size_t >( gridPtr_.nofVtxParam_ ) );
          for( double &v : vtxData_[ idSet_.id( entity ) ] )
            totalSize += dataSize( v );
        }

        if( Entity::codimension == 1 )
        {
          const auto bndData = bndData_.find( idSet_.id( entity ) );
          if( bndData != bndData_.end() )
            totalSize += dataSize( bndData->second.first ) + dataSize( bndData->second.second );
        }

        return totalSize;
      }

      template< class Buffer, class Entity >
      void gather ( Buffer &buffer, const Entity &entity ) const
      {
        // do not use a switch statement, because dimension == 1 is possible
        if( (Entity::codimension == 0) && (gridPtr_.nofElParam_ > 0) )
        {
          assert( elData_[ idSet_.id( entity ) ].size() == static_cast< std::size_t >( gridPtr_.nofElParam_ ) );
          for( double &v : elData_[ idSet_.id( entity ) ] )
            write( buffer, v );
        }

        if( (Entity::codimension == dimension) && (gridPtr_.nofVtxParam_ > 0) )
        {
          assert( vtxData_[ idSet_.id( entity ) ].size() == static_cast< std::size_t >( gridPtr_.nofVtxParam_ ) );
          for( double &v : vtxData_[ idSet_.id( entity ) ] )
            write( buffer, v );
        }

        if( Entity::codimension == 1 )
        {
          const auto bndData = bndData_.find( idSet_.id( entity ) );
          if( bndData != bndData_.end() )
          {
            write( buffer, bndData->second.first );
            write( buffer, bndData->second.second );
          }
        }
      }

      template< class Buffer, class Entity >
      void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
      {
        // do not use a switch statement, because dimension == 1 is possible
        if( (Entity::codimension == 0) && (gridPtr_.nofElParam_ > 0) )
        {
          auto &p = elData_[ idSet_.id( entity ) ];
          p.resize( gridPtr_.nofElParam_ );
          for( double &v : p )
            read( buffer, v, n );
        }

        if( (Entity::codimension == dimension) && (gridPtr_.nofVtxParam_ > 0) )
        {
          auto &p = vtxData_[ idSet_.id( entity ) ];
          p.resize( gridPtr_.nofVtxParam_ );
          for( double &v : p )
            read( buffer, v, n );
        }

        if( (Entity::codimension == 1) && (n > 0) )
        {
          auto &bndData = bndData_[ idSet_.id( entity ) ];
          read( buffer, bndData.first, n );
          read( buffer, bndData.second, n );
        }

        assert( n == 0 );
      }

    private:
      template< class T >
      static std::enable_if_t< std::is_trivially_copyable< T >::value, std::size_t > dataSize ( const T & /* value */ )
      {
        return sizeof( T );
      }

      static std::size_t dataSize ( const std::string &s )
      {
        return dataSize( s.size() ) + s.size();
      }

      template< class Buffer, class T >
      static std::enable_if_t< std::is_trivially_copyable< T >::value > write ( Buffer &buffer, const T &value )
      {
        std::array< char, sizeof( T ) > bytes;
        std::memcpy( bytes.data(), &value, sizeof( T ) );
        for( char &b : bytes )
          buffer.write( b );
      }

      template< class Buffer >
      static void write ( Buffer &buffer, const std::string &s )
      {
        write( buffer, s.size() );
        for( const char &c : s )
          buffer.write( c );
      }

      template< class Buffer, class T >
      static std::enable_if_t< std::is_trivially_copyable< T >::value > read ( Buffer &buffer, T &value, std::size_t &n )
      {
        assert( n >= sizeof( T ) );
        n -= sizeof( T );

        std::array< char, sizeof( T ) > bytes;
        for( char &b : bytes )
          buffer.read( b );
        std::memcpy( &value, bytes.data(), sizeof( T ) );
      }

      template< class Buffer >
      static void read ( Buffer &buffer, std::string &s, std::size_t &n )
      {
        std::size_t size;
        read( buffer, size, n );
        s.resize( size );

        assert( n >= size );
        n -= size;

        for( char &c : s )
          buffer.read( c );
      }

      GridPtr &gridPtr_;
      const typename GridType::LocalIdSet &idSet_;
      mutable std::map< typename GridType::LocalIdSet::IdType, std::vector< double > > elData_, vtxData_;
      mutable std::map< typename GridType::LocalIdSet::IdType, std::pair< int, DGFBoundaryParameter::type > > bndData_;
    };

    // grid auto pointer
    mutable mygrid_ptr gridPtr_;
    // element and vertex parameters
    std::vector< std::vector< double > > elParam_;
    std::vector< std::vector< double > > vtxParam_;
    std::vector< DGFBoundaryParameter::type > bndParam_;
    std::vector< int > bndId_;
    std::vector< double > emptyParam_;

    int nofElParam_;
    int nofVtxParam_;
    bool haveBndParam_;
  }; // end of class GridPtr

} // end namespace Dune

#endif
