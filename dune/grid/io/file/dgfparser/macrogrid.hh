// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_MACROGRID_HH
#define DUNE_DGF_MACROGRID_HH


#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/parser.hh>


namespace Dune
{
  // forward declarations
  // --------------------
  class DuneGridFormatParser;

  class MacroGrid
    : protected DuneGridFormatParser
  {
    template< class GridType >
    friend struct DGFGridFactory;

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
    static int rank( [[maybe_unused]] MPICommunicatorType MPICOMM )
    {
      int rank = 0;
#if HAVE_MPI
      MPI_Comm_rank( MPICOMM, &rank );
#endif
      return rank;
    }
    static int size( [[maybe_unused]] MPICommunicatorType MPICOMM )
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
        \code
        static GridType *
        generate ( MacroGrid &macroGrid, const char *filename,
                   MPIHelper :: MPICommunicator comm = MPIHelper :: getCommunicator() );
        \endcode
     */
    template< class GridType >
    struct Impl
    {
      static GridType* generate(MacroGrid& mg,
                                const char* filename, MPICommunicatorType MPICOMM = MPIHelper::getCommunicator() )
      {
        // make assertion depend on the template argument but always evaluate to false
        static_assert( GridType::dimension<0,"DGF grid factory missing. Did you forget to add the corresponding dgf header or config.h?");
      }
    };

    const char* filename_;
    MPICommunicatorType MPICOMM_;
  };

} // end namespace Dune

#endif
