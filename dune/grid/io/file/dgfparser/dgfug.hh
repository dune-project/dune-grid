// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_IO_FILE_DGFPARSER_DGFUG_HH
#define DUNE_GRID_IO_FILE_DGFPARSER_DGFUG_HH

//- C++ includes
#include <fstream>
#include <istream>
#include <string>
#include <vector>

//- dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

//- dune-grid includes
#include <dune/grid/common/intersection.hh>
#include <dune/grid/uggrid.hh>

//- local includes
#include "dgfparser.hh"
#include "blocks/gridparameter.hh"


namespace Dune
{

  namespace dgf
  {

    // UGGridParameterBlock
    // --------------------

    struct UGGridParameterBlock
      : public GridParameterBlock
    {
      /** \brief constructor taking istream */
      explicit UGGridParameterBlock ( std::istream &input );

      /** \brief returns true if no closure should be used for UGGrid */
      bool noClosure () const { return noClosure_; }
      /** \brief returns true if no copies are made for UGGrid elements */
      bool noCopy () const { return noCopy_; }
      /** \brief returns heap size used on construction of the grid */
      size_t heapSize () const { return heapSize_; }

    protected:
      bool noClosure_;  // no closure for UGGrid
      bool noCopy_;     // no copies  for UGGrid
      size_t heapSize_; // heap size  for UGGrid
    };

  } // namespace dgf



#if HAVE_DUNE_UGGRID
  template< int dim >
  struct DGFGridInfo< UGGrid< dim > >
  {
    static int refineStepsForHalf ()
    {
      return 1;
    }

    static double refineWeight ()
    {
      return -1.;
    }
  };



  // DGFGridFactory< UGGrid< dim > >
  // -------------------------------

  template< int dim >
  struct DGFGridFactory< UGGrid< dim > >
  {
    /** \brief grid type */
    typedef UGGrid< dim > Grid;
    /** \brief grid dimension */
    static const int dimension = dim;
    /** \brief MPI communicator type */
    typedef MPIHelper::MPICommunicator MPICommunicatorType;

    /** \brief constructor taking istream */
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : grid_( 0 ),
        factory_(),
        dgf_( rank( comm ), size( comm ) )
    {
      generate( input );
    }

    /** \brief constructor taking filename */
    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : grid_( 0 ),
        factory_(),
        dgf_( rank( comm ), size( comm ) )
    {
      std::ifstream input( filename.c_str() );
      if ( !input )
        DUNE_THROW( DGFException, "Error: Macrofile " << filename << " not found" );
      generate( input );
    }

    /** \brief return grid */
    Grid *grid ()
    {
      return grid_;
    }

    /** \brief please doc me */
    template< class GG, class II >
    bool wasInserted ( const Dune::Intersection< GG, II > &intersection ) const
    {
      return factory_.wasInserted( intersection );
    }

    /** \brief will return boundary segment index */
    template< class GG, class II >
    int boundaryId ( const Dune::Intersection< GG, II > &intersection ) const
    {
      return intersection.boundarySegmentIndex();
    }

    /** \brief return number of parameters */
    template< int codim >
    int numParameters () const
    {
      if( codim == 0 )
        return dgf_.nofelparams;
      else if( codim == dimension )
        return dgf_.nofvtxparams;
      else
        return 0;
    }

    /** \brief return number of parameters */
    template< class Entity >
    int numParameters ( const Entity & ) const
    {
      return numParameters< Entity::codimension >();
    }

    /** \brief return parameter for codim 0 entity */
    std::vector< double > &parameter ( const typename Grid::template Codim< 0 >::Entity &element )
    {
      if( numParameters< 0 >() <= 0 )
      {
        DUNE_THROW( InvalidStateException,
                    "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
      }
      return dgf_.elParams[ factory_.insertionIndex( element ) ];
    }

    /** \brief return parameter for vertex */
    std::vector< double > &parameter ( const typename Grid::template Codim< dimension >::Entity &vertex )
    {
      if( numParameters< dimension >() <= 0 )
      {
        DUNE_THROW( InvalidStateException,
                    "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
      }
      return dgf_.vtxParams[ factory_.insertionIndex( vertex ) ];
    }

    /** \brief UGGrid does not support boundary parameters */
    bool haveBoundaryParameters () const
    {
      return dgf_.haveBndParameters;
    }

    /** \brief return invalid value */
    template< class GG, class II >
    const DGFBoundaryParameter::type &boundaryParameter ( const Dune::Intersection< GG, II > &intersection ) const
    {
      typedef Dune::Intersection< GG, II > Intersection;
      typename Intersection::Entity entity = intersection.inside();
      const int face = intersection.indexInInside();

      auto refElem = referenceElement< double, dimension >( entity.type() );
      int corners = refElem.size( face, 1, dimension );
      std::vector< unsigned int > bound( corners );
      for( int i = 0; i < corners; ++i )
      {
        const int k = refElem.subEntity( face, 1, i, dimension );
        bound[ i ] = factory_.insertionIndex( entity.template subEntity< dimension >( k ) );
      }

      DuneGridFormatParser::facemap_t::key_type key( bound, false );
      const DuneGridFormatParser::facemap_t::const_iterator pos = dgf_.facemap.find( key );
      if( pos != dgf_.facemap.end() )
        return dgf_.facemap.find( key )->second.second;
      else
        return DGFBoundaryParameter::defaultValue();
    }

  private:
    // create grid
    void generate ( std::istream &input );

    // return rank
    static int rank( MPICommunicatorType MPICOMM )
    {
      int rank = 0;
#if HAVE_MPI
      MPI_Comm_rank( MPICOMM, &rank );
#endif
      return rank;
    }

    // return size
    static int size( MPICommunicatorType MPICOMM )
    {
      int size = 1;
#if HAVE_MPI
      MPI_Comm_size( MPICOMM, &size );
#endif
      return size;
    }

    Grid *grid_;
    GridFactory< UGGrid< dim > > factory_;
    DuneGridFormatParser dgf_;
  };
#endif // #if HAVE_DUNE_UGGRID

} // namespace Dune

#endif // #ifndef DUNE_GRID_IO_FILE_DGFPARSER_DGFUG_HH
