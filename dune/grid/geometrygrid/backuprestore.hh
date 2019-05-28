// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_BACKUPRESTORE_HH
#define DUNE_GEOGRID_BACKUPRESTORE_HH

#include <type_traits>

#include <dune/common/exception.hh>
#include <dune/grid/common/backuprestore.hh>

#include <dune/grid/geometrygrid/declaration.hh>
#include <dune/grid/geometrygrid/capabilities.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // BackupRestoreFacilities
    // -----------------------

    template< class Grid, bool hasBackupRestoreFacilities = Capabilities::hasBackupRestoreFacilities< Grid > ::v >
    class BackupRestoreFacilities
    {};

    template< class Grid >
    class BackupRestoreFacilities< Grid, true >
    {
      typedef BackupRestoreFacilities< Grid, true > This;

    protected:
      BackupRestoreFacilities ()
      {}

    private:
      BackupRestoreFacilities ( const This & );
      This &operator= ( const This & );

    protected:
      const Grid &asImp () const
      {
        return static_cast< const Grid & >( *this );
      }

      Grid &asImp ()
      {
        return static_cast< Grid & >( *this );
      }
    };

  } // namespace GeoGrid



  // BackupRestoreFacility for GeometryGrid
  // --------------------------------------

  template< class HostGrid, class CoordFunction, class Allocator >
  struct BackupRestoreFacility< GeometryGrid< HostGrid, CoordFunction, Allocator > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction, Allocator > Grid;
    typedef BackupRestoreFacility< HostGrid > HostBackupRestoreFacility;

    template <class Output>
    static void backup ( const Grid &grid, const Output &filename_or_stream )
    {
      // notice: We should also backup the coordinate function
      HostBackupRestoreFacility::backup( grid.hostGrid(), filename_or_stream );
    }

    /// \brief Return the grid from file or stream
    /**
     * NOTE: assumes that the CoordFunction can be default-constructed.
     **/
    template <class Input>
    static Grid *restore ( const Input &filename_or_stream )
    {
      return restore_impl(filename_or_stream, std::is_default_constructible<CoordFunction>{});
    }

    /// \brief Return the grid from file or stream by creating the CoordFunction from a creator
    /**
     * This interface extension of BackupRestoreFacility adds the possibility to create the
     * CoordFunction from the restored grid. The creator should return an object that can be passed
     * to the constructor of GeometryGrid and preserves livetime.
     *
     * Example:
     *
       \code
       using HostGrid = ...;
       using Grid = GeometryGrid<HostGrid, CoordFunction>;
       std::unique_ptr<Grid> grid(BackupRestoreFacility<Grid>::restore("grid.ext",
         [](auto const& grid) { return new CoordFunction(grid.leafGridView()); });
       \endcode
     *
     **/
    template <class Input, class CoordFunctionCreator>
    static Grid *restore ( const Input &filename_or_stream, CoordFunctionCreator creator )
    {
      HostGrid *hostGrid = HostBackupRestoreFacility::restore( filename_or_stream );
      return new Grid( hostGrid, creator(*hostGrid) );
    }


  private:
    template <class Input>
    static Grid *restore_impl ( const Input &filename_or_stream, std::true_type )
    {
      // notice: We should also restore the coordinate function
      HostGrid *hostGrid = HostBackupRestoreFacility::restore( filename_or_stream );
      CoordFunction *coordFunction = new CoordFunction();
      return new Grid( hostGrid, coordFunction );
    }

    template <class Input>
    static Grid *restore_impl ( const Input &filename_stream, std::false_type )
    {
      DUNE_THROW(NotImplemented,
        "Restoring a GeometryGrid with a CoordFunction that is not default-constructible is not implemented. Use restore(input, creator) instead.");
      return nullptr;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_BACKUPRESTORE_HH
