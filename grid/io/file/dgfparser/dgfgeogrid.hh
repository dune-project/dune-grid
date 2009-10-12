// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFGEOGRID_HH
#define DUNE_DGFGEOGRID_HH

#include <dune/common/typetraits.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfprojectionblock.hh>

namespace Dune
{

  /************************************************************************
  * Warning:
  * Reading DGF files directly into a GeometryGrid is a dirty hack for
  * two reasons:
  *   1) The host grid and coordinate function are never deleted (dangling
  *      pointers).
  *   2) The coordinate function has to provide a default constructor
  ************************************************************************/


  // DGFCoordFunction
  // ----------------

  template< int dimD, int dimR >
  class DGFCoordFunction
    : public AnalyticalCoordFunction< double, dimD, dimR, DGFCoordFunction< dimD, dimR > >
  {
    typedef DGFCoordFunction< dimD, dimR > This;
    typedef AnalyticalCoordFunction< double, dimD, dimR, This > Base;

  public:
    typedef typename Base::DomainVector DomainVector;
    typedef typename Base::RangeVector RangeVector;

    typedef dgf::ProjectionBlock::Expression Expression;

    DGFCoordFunction ( const Expression *expression )
      : expression_( expression )
    {}

    void evaluate ( const DomainVector &x, RangeVector &y ) const
    {
      std::vector< double > vx( dimD );
      std::vector< double > vy;
      for( int i = 0; i < dimD; ++i )
        vx[ i ] = x[ i ];
      expression_->evaluate( vx, vy );
      assert( vy.size() == size_t( dimR ) );
      for( int i = 0; i < dimR; ++i )
        y[ i ] = vy[ i ];
    }

  private:
    const Expression *expression_;
  };



  // MacroGrid::Impl for GeomegryGrid
  // --------------------------------

  /** \cond */
  template< class HostGrid, class CoordFunction >
  struct MacroGrid::Impl< GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef MPIHelper::MPICommunicator MPICommunicator;

    static const bool isDiscreteCoordFunction
      = GeoGrid::isDiscreteCoordFunctionInterface< typename CoordFunction::Interface >::value;

  private:
    template< bool >
    struct DiscreteFactory;

    template< bool >
    struct AnalyticalFactory;

    typedef typename SelectType< isDiscreteCoordFunction, DiscreteFactory< true >, AnalyticalFactory< false > >::Type
    Factory;

  public:
    static GeometryGrid< HostGrid, CoordFunction > *
    generate ( MacroGrid &macroGrid, const char *filename,
               MPICommunicator communicator = MPIHelper::getCommunicator() );
  };


  template< class HostGrid, class CoordFunction >
  template< bool >
  struct MacroGrid::Impl< GeometryGrid< HostGrid, CoordFunction > >::DiscreteFactory
  {
    static CoordFunction *create ( const HostGrid &hostGrid )
    {
      return new CoordFunction( hostGrid );
    }
  };


  template< class HostGrid, class CoordFunction >
  template< bool >
  struct MacroGrid::Impl< GeometryGrid< HostGrid, CoordFunction > >::AnalyticalFactory
  {
    static CoordFunction *create ( const HostGrid &hostGrid )
    {
      return new CoordFunction;
    }
  };


  template< class HostGrid, class CoordFunction >
  inline GeometryGrid< HostGrid, CoordFunction > *
  MacroGrid::Impl< GeometryGrid< HostGrid, CoordFunction > >
  ::generate ( MacroGrid &macroGrid, const char *filename, MPICommunicator communicator )
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;
    typedef MacroGrid::Impl< HostGrid > HostImpl;

    HostGrid *hostGrid = HostImpl::generate( macroGrid, filename, communicator );
    assert( hostGrid != 0 );
    CoordFunction *coordFunction = Factory::create( *hostGrid );
    return new Grid( *hostGrid, *coordFunction );
  }



  // MacroGrid::Impl for GeomegryGrid with DGFCoordFunction
  // ------------------------------------------------------

  template< class HostGrid, int dimD, int dimR >
  struct MacroGrid::Impl< GeometryGrid< HostGrid, DGFCoordFunction< dimD, dimR > > >
  {
    typedef MPIHelper::MPICommunicator MPICommunicator;

    static GeometryGrid< HostGrid, DGFCoordFunction< dimD, dimR > > *
    generate ( MacroGrid &macroGrid, const char *filename,
               MPICommunicator communicator = MPIHelper::getCommunicator() );
  };


  template< class HostGrid, int dimD, int dimR >
  inline GeometryGrid< HostGrid, DGFCoordFunction< dimD, dimR > > *
  MacroGrid::Impl< GeometryGrid< HostGrid, DGFCoordFunction< dimD, dimR > > >
  ::generate ( MacroGrid &macroGrid, const char *filename, MPICommunicator communicator )
  {
    typedef DGFCoordFunction< dimD, dimR > CoordFunction;
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;
    typedef MacroGrid::Impl< HostGrid > HostImpl;

    HostGrid *hostGrid = HostImpl::generate( macroGrid, filename, communicator );
    assert( hostGrid != 0 );

    std::ifstream file( filename );
    dgf::ProjectionBlock projectionBlock( file, dimR );
    const typename CoordFunction::Expression *expression = projectionBlock.function( "coordfunction" );
    if( expression == 0 )
      DUNE_THROW( DGFException, "no coordfunction specified in DGF file." );
    CoordFunction *coordFunction = new CoordFunction( expression );
    return new Grid( *hostGrid, *coordFunction );
  }



  // DGFGridInfo for GeometryGrid
  // ----------------------------

  template< class HostGrid, class CoordFunction >
  struct DGFGridInfo< GeometryGrid< HostGrid, CoordFunction > >
  {
    static int refineStepsForHalf ()
    {
      return DGFGridInfo< HostGrid >::refineStepsForHalf();
    }

    static double refineWeight ()
    {
      return -1.0;
    }
  };
  /** \endcond */

}

#endif // #ifndef DUNE_DGFGEOGRID_HH
