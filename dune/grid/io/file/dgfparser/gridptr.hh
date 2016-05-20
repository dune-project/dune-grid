// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_GRIDPTR_HH
#define DUNE_DGF_GRIDPTR_HH

#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>

//- Dune includes
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/intersection.hh>

#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/dgfparser/entitykey.hh>
#include <dune/grid/io/file/dgfparser/parser.hh>

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
      mygrid_ptr( const mygrid_ptr& other ) { assignObj( other ); }
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
      DGFGridFactory< GridType > dgfFactory( filename, comm );
      initialize( dgfFactory );
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
    GridPtr( const GridPtr &org )
      : gridPtr_(org.gridPtr_),
        elParam_(org.elParam_),
        vtxParam_(org.vtxParam_),
        bndParam_(org.bndParam_),
        bndId_(org.bndId_),
        emptyParam_( org.emptyParam_ ),
        nofElParam_(org.nofElParam_),
        nofVtxParam_(org.nofVtxParam_),
        haveBndParam_(org.haveBndParam_)
    {}

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

    void loadBalance()
    {
      if ( gridPtr_->comm().size() == 1 )
        return;
      int params = nofElParam_ + nofVtxParam_;
      if ( haveBndParam_ )
        params += 1;
      if ( gridPtr_->comm().max( params ) > 0 )
      {
        DataHandle dh(*this);
        gridPtr_->loadBalance( dh.interface() );
        gridPtr_->communicate( dh.interface(), InteriorBorder_All_Interface,ForwardCommunication);
      } else
      {
        gridPtr_->loadBalance();
      }
    }

  protected:
    void initialize ( DGFGridFactory< GridType > &dgfFactory )
    {
      gridPtr_ = mygrid_ptr( dgfFactory.grid() );

      typedef typename GridType::LevelGridView GridView;
      GridView gridView = gridPtr_->levelGridView( 0 );
      const typename GridView::IndexSet &indexSet = gridView.indexSet();

      nofElParam_ = dgfFactory.template numParameters< 0 >();
      nofVtxParam_ = dgfFactory.template numParameters< dimension >();
      haveBndParam_ = dgfFactory.haveBoundaryParameters();

      if ( nofElParam_ > 0 )
        elParam_.resize( indexSet.size(0) );
      if ( nofVtxParam_ > 0 )
        vtxParam_.resize( indexSet.size(dimension) );
      bndId_.resize( indexSet.size(1) );
      if ( haveBndParam_ )
        bndParam_.resize( gridPtr_->numBoundarySegments() );

      const PartitionIteratorType partType = Interior_Partition;
      typedef typename GridView::template Codim< 0 >::template Partition< partType >::Iterator Iterator;
      const Iterator enditer = gridView.template end< 0, partType >();
      for( Iterator iter = gridView.template begin< 0, partType >(); iter != enditer; ++iter )
      {
        const typename Iterator::Entity &el = *iter;
        if ( nofElParam_ > 0 ) {
          std::swap( elParam_[ indexSet.index(el) ], dgfFactory.parameter(el) );
          assert( elParam_[ indexSet.index(el) ].size()  == (size_t)nofElParam_ );
        }
        if ( nofVtxParam_ > 0 )
        {
          const unsigned int subEntities = el.subEntities(dimension);
          for ( unsigned int v = 0; v < subEntities; ++v)
          {
            typename GridView::IndexSet::IndexType index = indexSet.subIndex(el,v,dimension);
            if ( vtxParam_[ index ].empty() )
              std::swap( vtxParam_[ index ], dgfFactory.parameter(el.template subEntity<dimension>(v) ) );
            assert( vtxParam_[ index ].size()  == (size_t)nofVtxParam_ );
          }
        }
        if ( el.hasBoundaryIntersections() )
        {
          typedef typename GridView::IntersectionIterator IntersectionIterator;
          const IntersectionIterator iend = gridView.iend(el);
          for( IntersectionIterator iiter = gridView.ibegin(el); iiter != iend; ++iiter )
          {
            const typename IntersectionIterator::Intersection &inter = *iiter;
            // dirty hack: check for "none" to make corner point grid work
            if ( inter.boundary() && !inter.type().isNone() )
            {
              const int k = indexSet.subIndex(el,inter.indexInInside(),1);
              bndId_[ k ] = dgfFactory.boundaryId( inter );
              if ( haveBndParam_ )
                bndParam_[ inter.boundarySegmentIndex() ] = dgfFactory.boundaryParameter( inter );
            }
          }
        }
      }
    }

    template <class Entity>
    std::vector< double > &params ( const Entity &entity )
    {
      typedef typename GridType::LevelGridView GridView;
      GridView gridView = gridPtr_->levelGridView( 0 );
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

    struct DataHandle : public CommDataHandleIF<DataHandle,double>
    {
      DataHandle( GridPtr& gridPtr) :
        gridPtr_(gridPtr),
        idSet_(gridPtr->localIdSet())
      {
        typedef typename GridType::LevelGridView GridView;
        GridView gridView = gridPtr_->levelGridView( 0 );
        const typename GridView::IndexSet &indexSet = gridView.indexSet();

        const PartitionIteratorType partType = Interior_Partition;
        typedef typename GridView::template Codim< 0 >::template Partition< partType >::Iterator Iterator;
        const Iterator enditer = gridView.template end< 0, partType >();
        for( Iterator iter = gridView.template begin< 0, partType >(); iter != enditer; ++iter )
        {
          const typename Iterator::Entity &el = *iter;
          if ( gridPtr_.nofElParam_ > 0 )
            std::swap( gridPtr_.elParam_[ indexSet.index(el) ], elData_[ idSet_.id(el) ] );
          if ( gridPtr_.nofVtxParam_ > 0 )
          {
            for ( unsigned int v = 0; v < el.subEntities(dimension); ++v)
            {
              typename GridView::IndexSet::IndexType index = indexSet.subIndex(el,v,dimension);
              if ( ! gridPtr_.vtxParam_[ index ].empty() )
                std::swap( gridPtr_.vtxParam_[ index ], vtxData_[ idSet_.subId(el,v,dimension) ] );
            }
          }
        }
      }

      ~DataHandle()
      {
        typedef typename GridType::LevelGridView GridView;
        GridView gridView = gridPtr_->levelGridView( 0 );
        const typename GridView::IndexSet &indexSet = gridView.indexSet();

        if ( gridPtr_.nofElParam_ > 0 )
          gridPtr_.elParam_.resize( indexSet.size(0) );
        if ( gridPtr_.nofVtxParam_ > 0 )
          gridPtr_.vtxParam_.resize( indexSet.size(dimension) );

        const PartitionIteratorType partType = All_Partition;
        typedef typename GridView::template Codim< 0 >::template Partition< partType >::Iterator Iterator;
        const Iterator enditer = gridView.template end< 0, partType >();
        for( Iterator iter = gridView.template begin< 0, partType >(); iter != enditer; ++iter )
        {
          const typename Iterator::Entity &el = *iter;
          if ( gridPtr_.nofElParam_ > 0 )
          {
            std::swap( gridPtr_.elParam_[ indexSet.index(el) ], elData_[ idSet_.id(el) ] );
            assert( gridPtr_.elParam_[ indexSet.index(el) ].size() == (unsigned int)gridPtr_.nofElParam_ );
          }
          if ( gridPtr_.nofVtxParam_ > 0 )
          {
            for ( unsigned int v = 0; v < el.subEntities(dimension); ++v)
            {
              typename GridView::IndexSet::IndexType index = indexSet.subIndex(el,v,dimension);
              if ( gridPtr_.vtxParam_[ index ].empty() )
                std::swap( gridPtr_.vtxParam_[ index ], vtxData_[ idSet_.subId(el,v,dimension) ] );
              assert( gridPtr_.vtxParam_[ index ].size() == (unsigned int)gridPtr_.nofVtxParam_ );
            }
          }
        }
      }

      CommDataHandleIF<DataHandle,double> &interface()
      {
        return *this;
      }

      bool contains (int dim, int codim) const
      {
        return (codim==dim || codim==0);
      }

      bool fixedsize (int dim, int codim) const
      {
        return false;
      }

      template<class EntityType>
      size_t size (const EntityType& e) const
      {
        return gridPtr_.nofParameters( (int) e.codimension);
      }

      template<class MessageBufferImp, class EntityType>
      void gather (MessageBufferImp& buff, const EntityType& e) const
      {
        const std::vector<double> &v = (e.codimension==0) ? elData_[idSet_.id(e)] : vtxData_[idSet_.id(e)];
        const size_t s = v.size();
        for (size_t i=0; i<s; ++i)
          buff.write( v[i] );
        assert( s == (size_t)gridPtr_.nofParameters(e.codimension) );
      }

      template<class MessageBufferImp, class EntityType>
      void scatter (MessageBufferImp& buff, const EntityType& e, size_t n)
      {
        std::vector<double> &v = (e.codimension==0) ? elData_[idSet_.id(e)] : vtxData_[idSet_.id(e)];
        v.resize( n );
        gridPtr_.setNofParams( e.codimension, n );
        for (size_t i=0; i<n; ++i)
          buff.read( v[i] );
      }

    private:
      typedef typename GridType::LocalIdSet IdSet;
      GridPtr &gridPtr_;
      const IdSet &idSet_;
      mutable std::map< typename IdSet::IdType, std::vector<double> > elData_, vtxData_;
    };

    // grid auto pointer
    mutable mygrid_ptr gridPtr_;
    // element and vertex parameters
    std::vector< std::vector< double > > elParam_;
    std::vector< std::vector< double > > vtxParam_;
    std::vector< DGFBoundaryParameter::type > bndParam_;
    std::vector< int > bndId_;
    std::vector < double > emptyParam_;

    int nofElParam_;
    int nofVtxParam_;
    bool haveBndParam_;
  }; // end of class GridPtr

} // end namespace Dune

#endif
