// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFGEOGRID_HH
#define DUNE_DGFGEOGRID_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/geometrygrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/blocks/projection.hh>
#include <dune/grid/utility/hostgridaccess.hh>
#include <dune/grid/common/intersection.hh>

#include <dune/grid/io/file/dgfparser/parser.hh>


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

  // External Forward Declarations
  // -----------------------------

  template< class GridImp, class IntersectionImp >
  class Intersection;



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



  // DGFCoordFunctionFactory
  // -----------------------

  template< class HostGrid, class CoordFunction,
      bool discrete = GeoGrid::isDiscreteCoordFunctionInterface< typename CoordFunction::Interface >::value >
  struct DGFCoordFunctionFactory;


  template< class HostGrid, class CoordFunction >
  struct DGFCoordFunctionFactory< HostGrid, CoordFunction, false >
  {
    static CoordFunction *create ( std::istream &, const HostGrid & )
    {
      return new CoordFunction;
    }
  };


  template< class HostGrid, class CoordFunction >
  struct DGFCoordFunctionFactory< HostGrid, CoordFunction, true >
  {
    static CoordFunction *create ( std::istream &, const HostGrid &hostGrid )
    {
      return new CoordFunction( hostGrid );
    }
  };


  template< class HostGrid, int dimD, int dimR >
  struct DGFCoordFunctionFactory< HostGrid, DGFCoordFunction< dimD, dimR >, false >
  {
    typedef DGFCoordFunction< dimD, dimR > CoordFunction;

    static CoordFunction *create ( std::istream &input, const HostGrid &hostGrid )
    {
      dgf::ProjectionBlock projectionBlock( input, dimR );
      const typename CoordFunction::Expression *expression = projectionBlock.function( "coordfunction" );
      if( expression == 0 )
        DUNE_THROW( DGFException, "no coordfunction specified in DGF file." );
      return new CoordFunction( expression );
    }
  };



  // DGFGridFactory for GeometryGrid
  // -------------------------------

  template< class HostGrid, class CoordFunction, class Allocator >
  struct DGFGridFactory< GeometryGrid< HostGrid, CoordFunction, Allocator > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction, Allocator > Grid;

    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicator;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;

    typedef DGFCoordFunctionFactory< HostGrid, CoordFunction > CoordFunctionFactory;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicator comm = MPIHelper::getCommunicator() )
      : dgfHostFactory_( input, comm ),
        grid_( 0 )
    {
      auto hostGrid = std::shared_ptr<HostGrid>(dgfHostFactory_.grid());
      assert( hostGrid != 0 );
      auto coordFunction = std::shared_ptr<CoordFunction>(CoordFunctionFactory::create( input, *hostGrid ));
      grid_ = new Grid( hostGrid, coordFunction );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicator comm = MPIHelper::getCommunicator() )
      : dgfHostFactory_( filename, comm ),
        grid_( 0 )
    {
      auto hostGrid = std::shared_ptr<HostGrid>(dgfHostFactory_.grid());
      assert( hostGrid != 0 );
      std::ifstream input( filename.c_str() );
      auto coordFunction = std::shared_ptr<CoordFunction>(CoordFunctionFactory::create( input, *hostGrid ));
      grid_ = new Grid( hostGrid, coordFunction );
    }

    Grid *grid () const
    {
      return grid_;
    }

    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return dgfHostFactory_.wasInserted( HostGridAccess< Grid >::hostIntersection( intersection ) );
    }

    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      return dgfHostFactory_.boundaryId( HostGridAccess< Grid >::hostIntersection( intersection ) );
    }

    template< int codim >
    int numParameters () const
    {
      return dgfHostFactory_.template numParameters< codim >();
    }

    // return true if boundary parameters found
    bool haveBoundaryParameters () const
    {
      return dgfHostFactory_.haveBoundaryParameters();
    }

    template< class GG, class II >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Dune::Intersection< GG, II > & intersection ) const
    {
      return dgfHostFactory_.boundaryParameter( HostGridAccess< Grid >::hostIntersection( intersection ) );
    }

    template< class Entity >
    std::vector< double > &parameter ( const Entity &entity )
    {
      return dgfHostFactory_.parameter( HostGridAccess< Grid >::hostEntity( entity ) );
    }

  private:
    DGFGridFactory< HostGrid > dgfHostFactory_;
    Grid *grid_;
  };



  // DGFGridInfo for GeometryGrid
  // ----------------------------

  template< class HostGrid, class CoordFunction, class Allocator >
  struct DGFGridInfo< GeometryGrid< HostGrid, CoordFunction, Allocator > >
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

}

#endif // #ifndef DUNE_DGFGEOGRID_HH
