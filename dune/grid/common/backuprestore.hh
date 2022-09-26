// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_BACKUPRESTORE_HH
#define DUNE_GRID_COMMON_BACKUPRESTORE_HH

#include <dune/common/exceptions.hh>

namespace Dune
{

  /** \class BackupRestoreFacility
   *  \brief facility for writing and reading grids
   *
   *  The BackupRestoreFacility allows writing hierarchic grids to disk and
   *  reading them back into another program.
   *
   *  It is guaranteed that all index sets and id sets are preserved by the
   *  backup / restore process.
   *  The result of restore is undefined, if the number of processes in a
   *  parallel program differs from the number of processes used on backup.
   *
   *  There are two pairs of backup / restore methods:
   *  - methods writing into one or more dedicated files,
   *  - methods operating on a std::stream.
   *  .
   *
   *  These techniques may not be mixed, i.e., you cannot write the grid into
   *  files and read it back from a stream or vice versa.
   *  While operating on a std::stream might be convenient, a grid written in
   *  another language than C++ might need to emulate this method by writing
   *  through a temporary file.
   *
   *  \note The backup and restore methods might not be implemented for each grid.
   *        In this case one can catch the Dune::NotImplemented exception and do
   *        something else.
   *
   *  \tparam  Grid  type of grid
   */
  template< class Grid >
  struct BackupRestoreFacility
  {
    /** \brief write a hierarchic grid to disk
     *
     *  \param[in]  grid        grid to write
     *  \param[in]  filename    filename of the file to create
     *
     *  \note This method might create multiple files based on the filename plus internal extension.
     */
    static void backup ( const Grid &grid, const std::string &filename )
    {
      DUNE_THROW( NotImplemented, "backup / restore not implemented." );
    }

    /** \brief write a hierarchic grid into a stream
     *
     *  \param[in]  grid        grid to write
     *  \param[in]  stream      std::stream to write the grid to
     *
     *  \note While operating on a std::ostream might be convenient, a grid
     *        written in another language than C++ might need to emulate this
     *        method by writing through a temporary file.
     */
    static void backup ( const Grid &grid, std::ostream &stream )
    {
      DUNE_THROW( NotImplemented, "backup / restore not implemented." );
    }

    /** \brief read a hierarchic grid from disk
     *
     *  \param[in]  filename    filename of the file to read
     *
     *  \returns a pointer to the grid (allocated by new)
     *
     *  \note This method might require multiple files based on the filename plus some extension.
     */
    static Grid *restore ( const std::string &filename )
    {
      DUNE_THROW( NotImplemented, "backup / restore not implemented." );
    }

    /** \brief read a hierarchic grid from a stream
     *
     *  \param[in]  stream      std::stream to read the grid from
     *
     *  \note While operating on a std::istream might be convenient, a grid
     *        written in another language than C++ might need to emulate this
     *        method by writing through a temporary file.
     */
    static Grid *restore ( std::istream &stream )
    {
      DUNE_THROW( NotImplemented, "backup / restore not implemented." );
    }
  };

  /** \brief BackupRestoreFacility taking const Grid as type and deriving from the
             version with the const.
   */
  template< class Grid >
  struct BackupRestoreFacility< const Grid >
    : public BackupRestoreFacility< Grid >
  {};

} // namespace Dune

#endif // #ifndef DUNE_GRID_COMMON_BACKUPRESTORE_HH
