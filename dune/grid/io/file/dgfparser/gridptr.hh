// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_GRIDPTR_HH
#define DUNE_DGF_GRIDPTR_HH

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <assert.h>

//- Dune includes
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/dgfparser/entitykey.hh>
#include <dune/grid/io/file/dgfparser/parser.hh>

#include <dune/grid/common/intersection.hh>

namespace Dune
{
  // forward declarations
  // --------------------
  template < class G >
  struct DGFGridFactory;

  template< class GridImp, template < class > class IntersectionImp >
  class Intersection;

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
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    static const int dimension = GridType::dimension;

    //! constructor given the name of a DGF file
    explicit GridPtr ( const std::string &filename,
                       MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : gridPtr_( 0 ),
        elParam_( 0 ),
        vtxParam_( 0 ),
        bndParam_( 0 ),
        bndId_( 0 ),
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
      : gridPtr_( 0 ),
        elParam_( 0 ),
        vtxParam_( 0 ),
        bndParam_( 0 ),
        bndId_( 0 ),
        nofElParam_( 0 ),
        nofVtxParam_( 0 ),
        haveBndParam_( false )
    {
      DGFGridFactory< GridType > dgfFactory( input, comm );
      initialize( dgfFactory );
    }

    //! Default constructor, creating empty GridPtr
    GridPtr()
      : gridPtr_(0),
        nofElParam_(0),
        nofVtxParam_(0),
        haveBndParam_( false )
    {}

    //! Constructor storing given pointer to internal auto pointer
    GridPtr( GridType *grd )
      : gridPtr_(grd),
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
        nofElParam_(org.nofElParam_),
        nofVtxParam_(org.nofVtxParam_),
        haveBndParam_(org.haveBndParam_)
    {}

    //! assignment of grid pointer
    GridPtr &operator= ( const GridPtr &org )
    {
      gridPtr_ = org.gridPtr_;
      elParam_ = org.elParam_;
      vtxParam_ = org.vtxParam_;
      bndParam_ = org.bndParam_;
      bndId_ = org.bndId_;
      nofElParam_ = org.nofElParam_;
      nofVtxParam_ = org.nofVtxParam_;
      haveBndParam_ = org.haveBndParam_;
      return *this;
    }

    //! assignment of pointer to internal auto pointer
    GridPtr & operator = (GridType * grd)
    {
      gridPtr_ = std::auto_ptr<GridType>(grd);
      nofVtxParam_ = 0;
      nofElParam_ = 0;
      haveBndParam_ = false;
      elParam_.resize(0);
      vtxParam_.resize(0);
      bndParam_.resize(0);
      bndId_.resize(0);
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
    GridType* release () {
      return gridPtr_.release();
    }

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
      return nofParamters( Entity::codimension );
    }

    //! get number of parameters defined for a given intersection
    template< class GridImp, template< class > class IntersectionImp >
    int nofParameters ( const Intersection< GridImp, IntersectionImp > & intersection ) const
    {
      return parameters( intersection ).size();
    }

    //! get parameters defined for each codim 0 und dim entity on the grid through the grid file
    template <class Entity>
    const std::vector< double > &parameters ( const Entity &entity ) const
    {
      switch( (int)Entity::codimension )
      {
      case 0 :
        if( nofElParam_ > 0 )
        {
          assert( (unsigned int)gridPtr_->leafView().indexSet().index( entity ) < elParam_.size() );
          return elParam_[ gridPtr_->leafView().indexSet().index( entity ) ];
        }
        break;
      case GridType::dimension :
        if( nofVtxParam_ > 0 )
        {
          assert( (unsigned int)gridPtr_->leafView().indexSet().index( entity ) < vtxParam_.size() );
          return vtxParam_[ gridPtr_->leafView().indexSet().index( entity ) ];
        }
        break;
      }
      return emptyParam_;
    }

    //! get parameters for intersection
    template< class GridImp, template< class > class IntersectionImp >
    const DGFBoundaryParameter::type & parameters ( const Intersection< GridImp, IntersectionImp > & intersection ) const
    {
      // if no parameters given return empty vecto
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
      gridPtr_ = std::auto_ptr< GridType >( dgfFactory.grid() );

      typedef typename GridType::LeafGridView GridView;
      GridView gridView = gridPtr_->leafView();
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
          for ( int v = 0; v < el.template count<dimension>(); ++v)
          {
            typename GridView::IndexSet::IndexType index = indexSet.subIndex(el,v,dimension);
            if ( vtxParam_[ index ].empty() )
              std::swap( vtxParam_[ index ], dgfFactory.parameter(*el.template subEntity<dimension>(v) ) );
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
      switch( (int)Entity::codimension )
      {
      case 0 :
        if( nofElParam_ > 0 ) {
          if ( gridPtr_->leafView().indexSet().index( entity ) >= elParam_.size() )
            elParam_.resize( gridPtr_->leafView().indexSet().index( entity ) );
          return elParam_[ gridPtr_->leafView().indexSet().index( entity ) ];
        }
        break;
      case GridType::dimension :
        if( nofVtxParam_ > 0 ) {
          if ( gridPtr_->leafView().indexSet().index( entity ) >= vtxParam_.size() )
            vtxParam_.resize( gridPtr_->leafView().indexSet().index( entity ) );
          return vtxParam_[ gridPtr_->leafView().indexSet().index( entity ) ];
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
        typedef typename GridType::LeafGridView GridView;
        GridView gridView = gridPtr_->leafView();
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
            for ( int v = 0; v < el.template count<dimension>(); ++v)
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
        typedef typename GridType::LeafGridView GridView;
        GridView gridView = gridPtr_->leafView();
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
            for ( int v = 0; v < el.template count<dimension>(); ++v)
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
        return gridPtr_.nofParameters(e.codimension);
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
    mutable std::auto_ptr<GridType> gridPtr_;
    // element and vertex parameters
    std::vector< std::vector< double > > elParam_;
    std::vector< std::vector< double > > vtxParam_;
    std::vector< DGFBoundaryParameter::type > bndParam_;
    std::vector< int > bndId_;
    int nofElParam_, nofVtxParam_;
    bool haveBndParam_;
    std::vector < double > emptyParam_;
  }; // end of class GridPtr

} // end namespace Dune

#endif
