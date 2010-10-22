// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MACROGRIDPARSER_HH
#define DUNE_MACROGRIDPARSER_HH

#include <iostream>
#include <fstream>

#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <memory>
#include <map>
#include <assert.h>
#include <cmath>

//- Dune includes
#include <dune/common/mpihelper.hh>
#include <dune/common/stdstreams.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>
//- local includes

#include "dgfexception.hh"
#include "entitykey.hh"
#include "dgfparserblocks.hh"

namespace Dune
{

  class DGFPrintInfo;

  //! \brief The %DuneGridFormatParser class: reads a DGF file and stores
  //! build information in vector structures used by the MacroGrid class.
  class DuneGridFormatParser
  {
    template< class GridType >
    friend class DGFGridFactory;
    template< class GridType >
    friend class DGFBaseFactory;
  public:
    typedef enum {Simplex,Cube,General} element_t;
    typedef enum {counterclockwise=1,clockwise=-1} orientation_t;

    //! constructor
    DuneGridFormatParser ( int rank, int size );

    static bool isDuneGridFormat ( std::istream & );

    //! \brief method which reads the dgf file
    //!
    //! fills the vtx,element, and bound vectors
    //! returns true if reading succeded
    bool readDuneGrid( std::istream &, int dimG, int dimW );

    bool readDuneGrid( std::istream &input, int dimG = -1 ) DUNE_DEPRECATED
    {
      return readDuneGrid( input, dimG, -1 );
    }

    //! method to write in Tetgen/Triangle Poly Format
    void writeTetgenPoly ( const std :: string &, std ::string &, std :: string & );
    void writeTetgenPoly ( std::ostream &out, const bool writeSegments = true );

  protected:
    // dimension of world and problem: set through the readDuneGrid() method
    int dimw, dimgrid;
    // vector of vertex coordinates
    std::vector < std::vector <double> > vtx;
    int nofvtx;
    int vtxoffset;
    double minVertexDistance; // min. L^1 distance of distinct points
    // vector of elements
    std :: vector< std :: vector< unsigned int > > elements;
    int nofelements;
    // vector of boundary segments + identifier
    std::vector < std::vector <int> > bound;
    int nofbound;
    // map to generate and find boundary segments
    typedef std :: map< DGFEntityKey< unsigned int >, int > facemap_t;
    facemap_t facemap;
    // set by generator depending on element type wanted
    element_t element;
    // set by the readDuneGrid method depending
    // on what type the elements were generated
    bool simplexgrid;
    // true if grid is generated using the intervall Block
    bool cube2simplex;
    // parameters on elements
    int nofvtxparams,nofelparams;
    std::vector<std::vector<double> > vtxParams,elParams;
    // write information about generation process
    DGFPrintInfo* info;
    int rank_;
    int size_;

    void generateBoundaries ( std::istream &, bool );

    // call to tetgen/triangle
    void generateSimplexGrid ( std::istream & );
    void readTetgenTriangle ( const std :: string & );

    // helper methods
    void removeCopies ();
    void setOrientation ( int use1, int use2,
                          orientation_t orientation=counterclockwise );
    void setRefinement ( int use1, int use2, int is1=-1, int is2=-1 );
    double testTriang ( int snr );

    std :: vector< double > &getElParam ( int i, std::vector< double > &coord );
    std :: vector< double > &getVtxParam ( int i, std::vector< double > &coord );

    static std::string temporaryFileName ();
  };


  class MacroGrid : protected DuneGridFormatParser
  {
    template< class GridType >
    friend class DGFGridFactory;

  public:
    typedef MPIHelper::MPICommunicator MPICommunicatorType;

  protected:
    //! constructor given the name of a DGF file
    MacroGrid(const char* filename, MPICommunicatorType MPICOMM = MPIHelper::getCommunicator())
      : DuneGridFormatParser( rank(MPICOMM), size(MPICOMM) )
        , filename_(filename)
        , MPICOMM_(MPICOMM) {}

    //! constructor given the name of a DGF file
    MacroGrid(MPICommunicatorType MPICOMM = MPIHelper::getCommunicator())
      : DuneGridFormatParser( rank(MPICOMM), size(MPICOMM) )
        , filename_(0)
        , MPICOMM_(MPICOMM) {}

    //! returns pointer to a new instance of type GridType created from a DGF file
    template <class GridType>
    inline GridType * createGrid ()
    {
      return Impl<GridType>::generate(*this,filename_,MPICOMM_);
    }
  private:
    static int rank( MPICommunicatorType MPICOMM )
    {
      int rank = 0;
#if HAVE_MPI
      MPI_Comm_rank( MPICOMM, &rank );
#endif
      return rank;
    }
    static int size( MPICommunicatorType MPICOMM )
    {
      int size = 1;
#if HAVE_MPI
      MPI_Comm_size( MPICOMM, &size );
#endif
      return size;
    }
    /** \brief container for the actual grid generation method
     *
     *  For each grid implementation to be used with the DGF parser, this class
     *  has to be specialized. It has to contain one static method of the
     *  following prototype:
     *  \code
     *  static GridType *
     *  generate ( MacroGrid &macroGrid, const char *filename,
     *             MPIHelper :: MPICommunicator comm = MPIHelper :: getCommunicator() );
     *  \endcode
     */
    template< class GridType >
    struct Impl
    {
      static GridType* generate(MacroGrid& mg,
                                const char* filename, MPICommunicatorType MPICOMM = MPIHelper::getCommunicator() )
      {
        // make assertion depend on the template argument but always evaluate to false
        dune_static_assert( GridType::dimension<0,"dgf grid factory missing - did you forget to add the corresponding dgf header or dgfgridtype.hh ?");
      }
    };

    const char* filename_;
    MPICommunicatorType MPICOMM_;
  };

  template <class G>
  struct DGFGridFactory
  {
    typedef G Grid;
    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;

  private:
    typedef typename Grid::template Codim< 0 >::Entity Element;
    typedef typename Grid::template Codim< dimension >::Entity Vertex;

  public:
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : macroGrid_( comm )
    {
      DUNE_THROW( DGFException, "DGF factories using old MacroGrid implementation"
                  "don't support creation from std::istream." );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : macroGrid_( filename.c_str(), comm )
    {
      grid_ = macroGrid_.template createGrid< Grid >();

      if( macroGrid_.nofelparams > 0 )
      {
        const size_t nofElements = macroGrid_.elements.size();
        for( size_t i = 0; i < nofElements; ++i )
        {
          std::vector< double > coord;

          DomainType p(0);
          const size_t nofCorners = macroGrid_.elements[i].size();
          for (size_t k=0; k<nofCorners; ++k)
            for (int j=0; j<DomainType::dimension; ++j)
              p[j]+=macroGrid_.vtx[macroGrid_.elements[i][k]][j];
          p/=double(nofCorners);

          elInsertOrder_.insert( std::make_pair( p, i ) );
        }
      }

      if( macroGrid_.nofvtxparams > 0 )
      {
        const size_t nofVertices = macroGrid_.vtx.size();
        for( size_t i = 0; i < nofVertices; ++i )
        {
          std::vector< double > coord;

          DomainType p;
          for( int k = 0; k < DomainType::dimension; ++k )
            p[ k ] = macroGrid_.vtx[i][k];

          vtxInsertOrder_.insert( std::make_pair( p, i ) );
        }
      }
    }
    Grid *grid()
    {
      return grid_;
    }
    template <class Intersection>
    bool wasInserted(const Intersection &intersection) const
    {
      return intersection.boundary();
    }
    template <class Intersection>
    int boundaryId(const Intersection &intersection) const
    {
      return intersection.boundaryId();
    }

    template< int codim >
    int numParameters () const
    {
      if( codim == 0 )
        return macroGrid_.nofelparams;
      else if( codim == dimension )
        return macroGrid_.nofvtxparams;
      else
        return 0;
    }

    std::vector<double>& parameter(const Element &element)
    {
      const typename Element::Geometry &geo = element.geometry();
      DomainType coord( geo.corner( 0 ) );
      for( int i = 1; i < geo.corners(); ++i )
        coord += geo.corner( i );
      coord /= double( geo.corners() );

      InsertOrderIterator it = elInsertOrder_.find( coord );
      if( it != elInsertOrder_.end() )
        return macroGrid_.elParams[ it->second ];
      assert(0);
      return emptyParam;
    }

    std::vector<double>& parameter(const Vertex &vertex)
    {
      const typename Vertex::Geometry &geo = vertex.geometry();
      DomainType coord( geo.corner( 0 ) );

      InsertOrderIterator it = vtxInsertOrder_.find( coord );
      if( it != vtxInsertOrder_.end() )
        return macroGrid_.vtxParams[ it->second ];
      return emptyParam;
    }

  private:
    typedef FieldVector<typename Grid::ctype,Grid::dimensionworld> DomainType;
    struct Compare
    {
      bool operator() ( const DomainType &a, const DomainType &b ) const
      {
        // returns true, if a < b; c[i] < -eps;
        const DomainType c = a - b;
        const double eps = 1e-8;

        for( int i = 0; i < DomainType::dimension; ++i )
        {
          if( c[ i ] <= -eps )
            return true;
          if( c[ i ] >= eps )
            return false;
        }
        return false;
      }
    };
    typedef std::map< DomainType, size_t, Compare > InsertOrderMap;
    typedef typename InsertOrderMap::const_iterator InsertOrderIterator;
    MacroGrid macroGrid_;
    Grid *grid_;
    InsertOrderMap elInsertOrder_;
    InsertOrderMap vtxInsertOrder_;
    std::vector<double> emptyParam;
  };

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
        bndId_( 0 ),
        nofElParam_( 0 ),
        nofVtxParam_( 0 )
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
        bndId_( 0 ),
        nofElParam_( 0 ),
        nofVtxParam_( 0 )
    {
      DGFGridFactory< GridType > dgfFactory( input, comm );
      initialize( dgfFactory );
    }

    //! Default constructor, creating empty GridPtr
    GridPtr()
      : gridPtr_(0),
        nofElParam_(0),
        nofVtxParam_(0)
    {}

    //! Constructor storing given pointer to internal auto pointer
    GridPtr( GridType *grd )
      : gridPtr_(grd),
        nofElParam_(0),
        nofVtxParam_(0)
    {}

    //! Copy constructor, copies internal auto pointer
    GridPtr( const GridPtr &org )
      : gridPtr_(org.gridPtr_),
        elParam_(org.elParam_),
        vtxParam_(org.vtxParam_),
        bndId_(org.bndId_),
        nofElParam_(org.nofElParam_),
        nofVtxParam_(org.nofVtxParam_)
    {}

    //! assignment of grid pointer
    GridPtr &operator= ( const GridPtr &org )
    {
      gridPtr_ = org.gridPtr_;
      elParam_ = org.elParam_;
      vtxParam_ = org.vtxParam_;
      bndId_ = org.bndId_;
      nofVtxParam_ = org.nofVtxParam_;
      nofElParam_ = org.nofElParam_;
      return *this;
    }

    //! assignment of pointer to internal auto pointer
    GridPtr & operator = (GridType * grd)
    {
      gridPtr_ = std::auto_ptr<GridType>(grd);
      nofVtxParam_ = 0;
      nofElParam_ = 0;
      elParam_.resize(0);
      vtxParam_.resize(0);
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

    void loadBalance()
    {
      if ( gridPtr_->comm().size() == 1 )
        return;
      int params = nofElParam_ + nofVtxParam_;
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
      if ( nofElParam_ > 0 )
        elParam_.resize( indexSet.size(0) );
      if ( nofVtxParam_ > 0 )
        vtxParam_.resize( indexSet.size(dimension) );
      bndId_.resize( indexSet.size(1) );

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
              bndId_[ indexSet.subIndex(el,inter.indexInInside(),1) ]
                = dgfFactory.boundaryId( inter );
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
    std::vector< int > bndId_;
    int nofElParam_,nofVtxParam_;
    std::vector < double > emptyParam_;
  }; // end of class GridPtr

}

/**  @addtogroup DuneGridFormatParser

     @brief Classes for reading a macrogrid file in the dune
     macrogrid format (dgf)

     @section General
     <!--=========-->
     The DGF format allows the simple description of macrogrids, which
     can be used to construct Dune grids independent of the underlying
     implementation. Due to the generality of the approach only a subset of
     the description language for each grid implementation is available through
     this interface.

     @section Usage
     <!---------------------------------------------->
     There are two ways of constructing Dune grids using DGF files:

     -# By including the file gridtype.hh from the
        dune/grid/io/file/dgfparser
        directory: \n
        by defining one of the symbols
        \c ALBERTAGRID ,
        \c ALUGRID_CUBE ,
        \c ALUGRID_SIMPLEX ,
        \c ALUGRID_CONFORM ,
        \c SGRID ,
        \c UGGRID , or
        \c YASPGRID
        and the integer
        \c GRIDDIM
        one obtains a definition of the type \c GridType and
        the required header files for the desired grid and the macrogrid parser
        are included.
     -# By directly including one of the files dgf*.hh
        from the dune/grid/io/file/dgfparser directory: \n
        in this case only the required header files for the
        desired grid and the macrogrid parser are included but no typedef is made.

     After a grid type (denoted with \c GridType in the following)
     is selected in some way, the grid can be constructed either by calling

     @code
       Dune::GridPtr<GridType> gridptr(filename, mpiHelper.getCommunicator() );
     @endcode
      or
     @code
       Dune::GridPtr<GridType> gridptr(filename);
     @endcode
     or
     @code
       Dune::GridPtr<GridType> gridptr;
       ...
       gridptr=Dune::GridPtr<GridType>(filename);
     @endcode

     where in the second and third example \c MPIHelper::getCommunicator()
     is selected as default
     value; \c filename is the name of the dgf file. This creates an
     auto pointer like object \c Dune::GridPtr<GridType> holding a pointer to a
     the grid instance described through the dgf file.
     Access to the grid is gained by calling the operator * of \c GridPtr.
     @code
       GridType & grid = *gridptr;
     @endcode
     Like in all auto_ptr like objects the grid instance is destroyed
     if the GridPtr goes out of scope. To be able to destruct the GridPtr
     object without losing the grid call the release method:
     @code
       GridType *grid;
       {
         GridPtr< GridType > gridPtr( filename.c_str(), mpiHelper.getCommunicator() );
         grid = gridPtr.release();
       }
     @endcode

     Remarks:
     -# The last argument in the constructor of the GridPtr should be of the type \c Dune::MPIHelper::MPICommunicator
        which defaults to \c MPI_COMM_WORLD for parallel runs or some default value for serial runs.
     -# If the file given through the first argument is not a dgf file
        a suitable constructor on the \c GridType class is called - if
        one is available.

     @section FORMAT Format Description
     <!--=========-->
     We assume in the following that the type
     \c GridType  is suitably defined, denoting a dune grid type
     in \c dimworld  space dimension.
     A point in dimworld space is called a vector.

     In general dgf files consists of one or more blocks, each starting with a
     keyword and ending with a line starting with a # symbol.
     In a block each full line is parsed from the beginning up to the
     first occupance of a \% symbol, which can be used to include comments.
     Trailing whitespaces are ignored
     during line parsing. Also the parsing is not case sensitive.

     Some example files are given below (\ref EXAMPLES).

     @subsection START First line
     <!---------------------------------------------->
     DGF files must start with the keyword \b DGF. Files passed to the
     grid parser not starting with this keyword are directly passed to a
     suitable constructor in the GridType class.

     @subsection BLOCKS Blocks
     <!---------------------------------------------->
     In the following all blocks are briefly described in alphabetical order.
     The construction of the grid is detailed below. In general lines
     are parse sequentially, with the exception of lines starting
     with a special keyword which are treated separately.

     - \b Boundarydomain  \n
       Each line consists of an integer greater than zero and two
       vectors describing an interval in \c dimworld space. The first entry of
       the line is a boundary identifier.
       - A special keyword is
         \b default
         followed by a positive integer which is to be used
         as a default boundary identifier.
     - \b Boundarysegments \n
       Each line consists of a positive integer denoting a boundary
       id and a set of vertex indices (see \b Vertex block)
       describing one boundary patch of the macrogrid.
       Note that a special ordering of the vertex numbers is not required.
     - \b Cube  \n
       Each line consists of \c dimworld<SUP>2</SUP> vertex indices
       (see \b Vertex block) describing one cube
       of the macrogrid.
       - If the ordering of the local vertices of the cube elements
         does not follow the %Dune
         @link GridReferenceElements reference element@endlink then
         the mapping between the local numbering used and the dune reference
         cube has to be prescribed by a line starting with the keyword \b map
         followed by the permutation of the numbers
         \c 0...2<SUP>dimworld</SUP> describing the mapping from the
         local vertex numbers in the reference elements used in the dgf file
         to the local vertex numbers in dune reference element.
       - By means of the special keyword \b parameters followed by a positive
         integer \c n it is possible to add \c n parameters to each element of the
         grid. These double valued parameters then have to be added to the definition
         lines for the elements behind the vertex numbers.
     - \b Interval \n
       Each interval is described by three lines in this block:
       the first two lines give the lower and upper vector of an interval,
       the third line consists of \c dimworld integers. used to determine the
       initial partition of the interval into equal sized cubes.
       More than one interval can be described in this fashion.
     - \b Simplex \n
       Each line consists of \c dimworld+1 vertex indices (see \b Vertex block)
       describing one simplex
       of the macrogrid.
       - Note that no ordering of local vertices is required.
       - Parameters can be added to each element using the \b parameters keyword
         as described for cube elements.
     - \b Simplexgenerator \n
       Using this block a simplex grid can be automatically generated using
       one of the freely available grid generation tools
       Tetgen (http://tetgen.berlios.de) for \c dimworld=3 or
       Triangle (http://www.cs.cmu.edu/~quake/triangle.html) for \c dimworld=2.
       For more detail see \ref Simplexgeneration.
     - \b Vertex \n
       Each line consists of a vector representing a vertex of the
       macrogrid. The vertices are consecutively numbered.
       leading to vertex indices which can be referenced in the
       \b Simplex, \b Cube, and \b Boundarysegment blocks.
       - By default the numbering starts from zero; to change this behaviour
         the keyword
         \b firstindex followed by a positive integer denoting the number
         of the first vertex can be used.
       - Using the \b parameters keyword it is possible to add a set of parameters
         to each vertex of the grid.
     - \b GridParameter \n
       For each grid implementation there is a set of parameters
       that can be passed via the GridParameter block to the momment of
       grid construction. For the specific parameters see the \ref
       DGFGridParameter section. See also the \b examplegrid5.dgf file
       for examples.
     - \b PeriodicFaceTransformation \n
       Each line describes an affine transformation that shall be used to
       glue grid boundaries together. The transformation is denoted as
       \em matrix + \em shift. The following 2d example describes a shift by
       the first unit vector:
       \code
   1 0, 0 1 + 1 0
       \endcode
     - \b Projection \n
       Each line either declares a function or it assigns one a boundary segment.
       Functions are declared as follows:\n
       <b>function</b> <em>function</em> <b>(</b> <em>variable</em> <b>)</b> <b>=</b> <em>expression</em>\n
       The following expressions are currently implemented:
       - <em>constant</em>
       - <em>variable</em>
       - <b>(</b> <em>expression</em> <b>)</b>
       - <em>function</em> <b>(</b> <em>expression</em> <b>)</b>
       - <em>expression</em> <b>[</b> <em>constant</em> <b>]</b>
       - <b>-</b> <em>expression</em>
       - <b>|</b> <em>expression</em> <b>|</b>
       - <em>expression</em> <b>+</b> <em>expression</em>
       - <em>expression</em> <b>-</b> <em>expression</em>
       - <em>expression</em> <b>*</b> <em>expression</em>
       - <em>expression</em> <b>/</b> <em>expression</em>
       - <em>expression</em> <b>**</b> <em>expression</em>
       - <b>sqrt</b> <em>expression</em>
       - <b>sin</b> <em>expression</em>
       - <b>cos</b> <em>expression</em>
       .
       Functions are assigned to boundary segments as follows:\n
       <b>segment</b> <em>vertex</em> ... <em>vertex</em> <em>function</em>\n
       Here, the given vertices identify the boundary segment.\n
       Another possibility is to set a default, assigned to all non-listed
       boundary segments:\n
       <b>default</b> <em>function</em>\n
       Note: Currently, the attached functions map global coordinates to global
             coordinates. This feature is only available with AlbertaGrid (Version 3.0
                 or above) or with ALUGrid.
     .

     @section CONSTR The Grid Construction Process
     <!---------------------------------------------->
     For simplicity we first describe how the grid is manually constructed,
     i.e., if no \b Simplexgenerator  block is found. Details on how
     Tetgen/Triangle can be used is detailed below (see \ref Simplexgeneration).
     How to access the element and vertex parameter is detailed in
     Section \ref PARAMETERS.
     The details of the construction are logged in the file
     \b dgfparser.log.

     The construction of the grid depends on the type of elements the
     given by \c GridType can handle.

     @subsection CONSTRCART Cartesian grids
     (Dune::SGrid , or
      Dune::YaspGrid )

     The grid is constructed using only the information from the
     first three lines of the \b Interval  block.

     @subsection CONSTRSIMPL Simplex grids
     (Dune::AlbertaGrid, or
      Dune::ALUSimplexGrid<3,3>, and
      Dune::ALUSimplexGrid<2,2>)

     The vertices and elements of the grid are constructed in the
     following three steps:
     -# The file is parsed for
        an \b Interval  block; if present a Cartesian grid is build for each
        interval defined in this block and
        each element is partitioned either into two triangles or into six
        tetrahedron. The user has to make sure that this process leads to
        a conforming grid.
     -# If no \b Interval
        block is found, the grid is generated using the information
        from the \b Vertex and the \b Simplex or \b Cube  blocks.
        If a non-empty \b Simplex block is found the element information is
        only taken from there; in the case where no \b Simplex block is found
        or in the case where it is empty, the \b Cube block is read and each
        cube is partitioned into simplex elements - in 3d this process will
        not always lead to a correct macro triangulation!
        Note that no specific ordering of the local numbering of each simplex
        is required but the cubes must either conform with the Dune reference
        element or the mapping must be prescribed using the map keyword.
     .
     If the simplex generation process was successful, boundary ids are
     assigned to all boundary faces of the macrogrid.
     -# If the macrogrid was constructed in the third step of the generation
        process detailed above (i.e. using the vertex block),
        then the file is first parsed for
        \b Boundarysegment  block. Here boundary ids can be
        individually assigned to each boundary segment of the macrogrid
        by specifying the vertex ids. Boundary segments can either be described
        as simplex elements or as cube elements.
     -# All Boundary segments which have not
        yet been assigned an identifier and lie inside the
        interval defined by the first line of the \b Boundarydomain
        block are assigned the corresponding id. This process is then
        repeated with the remaining boundary segments using the following
        intervals defined in the \b Boundarydomain  block.
        If after this process boundary segments without id remain,
        then the default id is used if one was specified in the
        \b Boundarydomain  block - if no default was given, the
        behavior of the parser is not well defined.
     .
     \b Remark:
     -# \c Bisection: \n
        the refinement edge is always chosen to be the longest edge.

     @subsection CONSTRCUBE Cube grids
     (Dune::ALUCubeGrid<3,3>)

     The grid is constructed using the information from the
     \b Interval  block, if present; otherwise the \b Vertex and \b Cube
     blocks are used.

     The boundary ids are assigned in the same manner as for simplex grids
     described above.

     @subsection CONSTRMIXED Mixed grids
     (\c Dune::UGGrid )

     Note that in this version only grids consisting of one element type are
     constructed even if the implemented grid allows for mixed elements.

     The vertices and elements of the grid are constructed in one of the
     following stages:
     -# The file is first parsed for
        an \b Interval  block; if present a Cartesian grid is build;
        a simplex grid is generated if a \c Simplex block is present
        otherwise a cube grid is constructed.
     -# If no \b Interval
        block is found, the grid is generated using the information
        from the \b Vertex block; cube elements are constructed if the
        \b Cube block is present otherwise the \b Simplex blocks is read.
        If both a \b Cube and a \b Simplex block is found, then only
        the element information from the \b Cube block is used and each
        cube is split into simplex elements so that
        a simplex grid is constructed.
     .
     Boundary ids are assigned in the same manner as described for
     simplex grids.

     @section PARAMETERS Accessing parameters
     In addition to the element and vertex information it is possible to add
     a set of parameters through the DGF file. Using the construction
     mechanism via the Dune::GridPtr<GridType> class these parameters can be
     accessed using the
     Dune::GridPtr<GridType>::parameters(const Entity& en)
     method. Depending on the codimentsion of \c en the method returns either
     the element or vertex parameters of that entity in the DGF file.
     The number of parameters for a given
     codimension can be retrived using the method
     Dune::GridPtr<GridType>::nofParameters.
     In the case of parallel runs only the process with rank zero reads the
     dgf file and thus the parameters are only available on that processor.
     Do not call loadBalance (or for that matter globalRefine, adapt) on the
     grid instance before retrieving the parameters. To simplify the handling
     of parameters in a parallel run the Dune::GridPtr<GridType>::loadBalance
     method can be used to distribute the macro grid and the parameters.
     Afterwards the parameters can be extracted on each processor using the
     method GridPtr<GridType>::parameters.

    <!---------------------------------------------->
   \section Simplexgeneration Using Tetgen/Triangle
       The freely available simplex grid generators are direcltly
       called via system
       call through the dgfparser.
       Therefore one should either add the path containing the executables of
       Triangle and/or Tetgen to the environment variable PATH or use the path
       option described below.
       One can either use a file in Tetgen/Triangle format to directly generate
       a macro grid or one can prescribe vertices and boundary faces
       which will be used to generate the grid directly in the DGF file:
       -# For the first approach use the token \b file to give a filename
          and the type of the file (e.g. node, mesh, ply...). Example:
          \code
            file name mesh
          \endcode
          If the filetype is given, it will be appended to \b name and this will be
          passed to Tetgen/Triangle. Additional parameters for the grid generators
          can be given by the the \b parameter token.
          If no file type is given it is assumed that in a previous run of
          Tetgen/Triangle
          files name.node and name.ele were generated and these will be used
          to described the vertices and elements of the Dune grid.
       -# In the second approach the \b vertex and the \b interval blocks
          (if present)
          are used to generate the vertices of the grid; the \b Cube and \b
          Simplex blocks are used for element information in the interior and
          the \b boundarysegment block is used for boundary information:
          - If only a \b vertex and/or \b interval block is found
            the resulting .node file is passed directly to Tetgen/Triangle and a
            tessellation of
            the convex hull of the points is generated.
          - A more detailed
            description of the domain is possible by using the \b
            Boundarysegment block together with the \b vertex block.
            Planar polyhedral boundary faces of the domain can be
            prescribed in this way.
            Note that in this case the whole
            domain boundary has to be defined (see the description of
            the .poly files in the Tetgen documentation).
            In 3d each polyhedral face (p0,..,pn) is automatically closed by adding
            the segment between p0 and pn. To define the whole boundary
            of a 2d domain using only one polyhedron the closing
            segment has to be added, i.e., (p0,..,pn,p0).
          - If a \b cube or \b simplex block is found the element information
            is also passed to tetgen/triangle together with the parameters - if
            given. Note that triangle can only handle one region atribute in
            its .poly files so that only the first parameter is the \b simplex
            or \b cube block can be retrived.
       .

       Some identifiers can be
       used to influence the quality of the generated mesh:
       -  \b max-area
          followed by a positive real number used as an upper bound for the
          area of all simplicies of the mesh.
       -  \b min-angle
          followed by a positive number. In 2d this limits the angles in
          resulting mesh from below; in 3d this bounds the radius-edge ratio
          from above.
       .
       In this case the grid is constructed using two calls to Tetgen/Triangle.           The details of the call are logged in the \b dgfparser.log file.

      \e Note: vertex parameters are interpolated in triangle but tetgen
                assigns a zero to all newly inserted vertices; therfore quality
                enhancement and vertex parameters should not be combined in 3d.
                On the other hand element parameters and boundary ids can be
                used together with quality enhancement.


       The remaining identifiers are
       -  The identifier \b path
          (followed by a path name) can be used
          to give search path for Triangle/Tetgen.
       -  The identifier \b display
          followed by 1 can be used to get a first impression
          of the resulting mesh using the visualization tools distributed with
          Triangle/Tetgen.
       .
       \e Note that parameters can be attached to the vertices and elements of
       the grid as described in the Triangle/Tetgen manual. Then can be
       retrieved as described in Section \ref PARAMETERS.

       Download
       - Tetgen http://tetgen.berlios.de
       - Triangle http://www.cs.cmu.edu/~quake/triangle.html



     <!---------------------------------------------->

   @section OPEN Work in progress
   -# There should be a mechanism to fix the desired refinement edge
      for simplex grids. An automatic selection is performed using the
      longest edge, but
      it should be possible to turn off this automatic selection.
   -# A callback to the user
      between parsing the DGF file and the construction of the dune grid;
      here e.g. a post-processing of the vertex coordinates could be included.

   @section EXAMPLES Examples
   <!--=========-->
   In two space dimensions:
   - \ref dgfexample1
   - \ref dgfexample2
   - \ref dgfexampleParam
   - \ref dgfexample3
   - \ref dgfexample6

   In three space dimensions:
   - \ref dgfexample4
   - \ref dgfexample5

   \section dgfexample1 Manual Grid Construction
   A tessellation of the unit square into six simplex entities.
   Some boundary segments on the lower and right
   boundary are given their own id the remaining are
   given a default id.

   @include examplegrid1s.dgf

   \image html  examplegrid1s.png "The resulting grid"

   A tessellation into cubes using the same vertices as before,

   @include examplegrid1c.dgf

   \image html  examplegrid1c.png "The resulting grid"

   Using the last input file with a simplex grid or by adding an empty
   \c Simplex block
   @code
   Simplex
 #
   @endcode
   leads to the following macro triangulation

   \image html  examplegrid1cs.png "The resulting grid"

   \section dgfexample2 Automated Grid Construction
   Automatic tessellation using Triangle,
   with vertices defined as in the example \ref dgfexample1 :

   @include examplegrid1gen.dgf

   \image html  examplegrid1gen.png "The resulting grid"

   The quality of the grid can be enhanced by adding the line
   @code
   min-angle 30
   @endcode
   in the \c Simplexgenerator block

   \image html  examplegrid1genangle.png "The resulting grid"

   Automatic tessellation using Triangle,
   with vertices are defined on a Cartesian grid with two additional
   vertices in the top right corner and one vertex outside the unit square.

   All boundary are given a default id.

   @include examplegrid2a.dgf

   \image html  examplegrid2a.png "The resulting grid"

   Adding some quality enhancement.
   The boundaries are numbered counterclockwise starting at the left
   boundary from one to four.

   @include examplegrid2b.dgf

   \image html  examplegrid2b.png "The resulting grid"

   Using both quality enhancement and a maximum area restriction.
   The bottom boundary is given the id 1 all other boundaries have
   id 2; here we do not use a default value.

   @include examplegrid2c.dgf

   \image html  examplegrid2c.png "The resulting grid"

   A similar grid is generated by prescribing the boundary of the domain:

   @include examplegrid2d.dgf

   \image html  examplegrid2d.png "The resulting grid"

   \section dgfexampleParam Using Parameters

   We use the same domain as in the previous example but include vertex and
   element parameters:

   @include examplegrid2e.dgf

   The results in piecewise constant data on the elements and piecewise linear
   date using the vertex parameters:

   \image html  examplegrid2e.png "The resulting grid with element and vertex parameters"

   \section dgfexample3 Interval Domain
   A automatic tessellation of the unit square using
   a Cartesian Grid.
   All boundaries have id 1.

   @include examplegrid5.dgf

   \image html  examplegrid5c.png "The resulting grid using SGrid<2,2>"
   \image html  examplegrid5s.png "The resulting grid using AlbertaGrid<2,2>"

   If UGGrid<2,2> is used the result would be the same as for SGrid<2,2>.
   If an empty \c Simplex Block
   @code
   Simplex
 #
   @endcode
   is added than the same simplex grid as for AlbertaGrid<2,2> would be
   constructed.

     <!---------------------------------------------->
   \section dgfexample6 Boundary Projections
   The following example shows a DGF file that projects 3 sides of a
   quadrilateral onto the surrounding circle:
   @include example-projection.dgf

   \image html  example-projection.png "The resulting grid using AlbertaGrid<2>"

     <!---------------------------------------------->
   \section dgfexample4 Grid Generation in 3d
   An automatic tessellation of the unit square using
   a Cartesian Grid is shown.
   All boundaries have id 1 except boundary segment on the lower boundary
   which have value 2.

   First we use only the interval block:

   @include examplegrid6.dgf

   \image html  examplegrid6c.png "The resulting grid using ALUCubeGrid<3,3>"
   \image html  examplegrid6s.png "The resulting grid using ALUSimplexGrid<3,3>"

   Now the vertices are still defined through the interval block; the simplicies
   are constructed using Tetgen (note the comment symbol \% in the
   \b Simplexgenerator block):

   @include examplegrid7.dgf

   \image html  examplegrid7.png "The resulting grid using ALUSimplexGrid<3,3>"

   Note that in the grid would be the same as in the above example if
   ALUCubeGrid<3,3> where used.

   Now we can also include some quality enhancement:

   First: \b min-angle = 1.2
          (remove the first \% in the \b Simplexgenerator block)
   \image html  examplegrid7angle.png "The resulting grid using ALUSimplexGrid<3,3>"

   Second: \b min-angle = 1.2 and \b max-area = 0.1
           (remove both \% in the \b Simplexgenerator block)
   \image html  examplegrid7area.png "The resulting grid using ALUSimplexGrid<3,3>"
     <!---------------------------------------------->

   This examples show different ways to define a grid for the following
   domain and boundary ids:
   \image html  examplegrid10.png "An example domain"

   \e Manual grid construction: Note that the reference element used
   to define the cube is not the one used in dune so that the
   mapping has to be given; furthermore the vertex index are numbered
   starting with 1.
   If a simplex grid is to be constructed
   some care must be taken in the numbering of the nodes so that a
   conforming grid is constructed.

   @include examplegrid10.dgf

   \image html  examplegrid10c.png "The resulting cube grid"
   \image html  examplegrid10s.png "The resulting simplex grid"

   \e Using multiple intervals: Here the boundary ids are not set
   correctly; this could be done using different boundarydomains.

   @include examplegrid11.dgf

   \image html  examplegrid11a.png "The resulting cube grid"
   \image html  examplegrid11b.png "The resulting simplex grid"

   \e By only defining the boundary faces and using Tetgen:

   @include examplegrid12.dgf

   Without the \b max-area and \b min-angle keywords the following
   grid is constructed:

   \image html  examplegrid12_1.png "The resulting grid"

   With the quality enhancement active:

   \image html  examplegrid12_2.png "The resulting grid"

   Finally we demonstrate the use of parameters. With a simple modification
   of the above example

   @include examplegrid10a.dgf

   we can define parameters on each element:

   \image html  examplegrid10a.png "The resulting grid with element parameters"

   The following example

   @include examplegrid10b.dgf

   defines parameters for each vertex of the grid, leading to a piecewise linear
   function

   \image html  examplegrid10b.png "The resulting grid with vertex parameters"

   <!---------------------------------------------->
   \section dgfexample5 Importing Grids written in a Tetgen/Triangle format

   Here a .mesh file used to generate a 3d simplex grid through Tetgen. The
   mesh file is taken from
   http://www-c.inria.fr/Eric.Saltel/download/download.php

   @include examplegrid9.dgf

   \image html  BBEETH1M.d_cut.png "The resulting grid using ALUSimplexGrid<3,3>"

   \image html  Orb_cut.png "The resulting grid using ALUSimplexGrid<3,3>"

   \image html  bunny.p65.param_skin.png "The resulting grid using ALUSimplexGrid<3,3>"

   \image html  pmdc.png "The resulting grid using ALUSimplexGrid<3,3>"

 **/
/*
      Dune::Alberta with \c dimworld=3: \n
        if Tetgen is used to construct a
        tetrahedral grid for Dune::Alberta then the bisection routine does
        not necessarily terminate. This problem does not occur
        if the grid is constructed using the \b Interval block.
 */

namespace Dune {

  /*! @brief Some simple static information for a given GridType
   */
  template <class GridType>
  struct DGFGridInfo
  {
    //! number of globalRefine steps needed to refuce h by 0.5
    static int refineStepsForHalf();
    //! relation between volume of children to volume of father.
    //! If this is not a constant the return value is -1
    static double refineWeight();
  };

} // end namespace Dune
#endif
