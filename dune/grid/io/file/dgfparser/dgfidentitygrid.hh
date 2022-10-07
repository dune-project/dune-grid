// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSER_DGFIDENTITYGRID_HH
#define DUNE_DGFPARSER_DGFIDENTITYGRID_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/identitygrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/utility/hostgridaccess.hh>

namespace Dune
{

  // DGFGridFactory for IdentityGrid
  // -------------------------------

  template< class HostGrid >
  struct DGFGridFactory< IdentityGrid< HostGrid > >
  {
    typedef IdentityGrid< HostGrid > Grid;

    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicator;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicator comm = MPIHelper::getCommunicator() )
      : dgfHostFactory_( input, comm ),
        grid_( 0 )
    {
      HostGrid *hostGrid = dgfHostFactory_.grid();
      assert( hostGrid != 0 );
      grid_ = new Grid( *hostGrid );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicator comm = MPIHelper::getCommunicator() )
      : dgfHostFactory_( filename, comm ),
        grid_( 0 )
    {
      HostGrid *hostGrid = dgfHostFactory_.grid();
      assert( hostGrid != 0 );
      std::ifstream input( filename.c_str() );
      grid_ = new Grid( *hostGrid );
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

    template< class Entity >
    std::vector< double > &parameter ( const Entity &entity )
    {
      return dgfHostFactory_.parameter( HostGridAccess< Grid >::hostEntity( entity ) );
    }

  private:
    DGFGridFactory< HostGrid > dgfHostFactory_;
    Grid *grid_;
  };



  // DGFGridInfo for IdGrid
  // ----------------------

  template< class HostGrid >
  struct DGFGridInfo< IdentityGrid< HostGrid > >
  {
    static int refineStepsForHalf ()
    {
      return DGFGridInfo< HostGrid >::refineStepsForHalf();
    }

    static double refineWeight ()
    {
      return DGFGridInfo< HostGrid >::refineWeight();
    }
  };

}

#endif // #ifndef DUNE_DGFPARSER_DGFIDENTITYGRID_HH
