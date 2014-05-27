// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ENTITYCOMMHELPER_HH
#define DUNE_ENTITYCOMMHELPER_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/common/unused.hh>

namespace Dune
{

  template< InterfaceType iftype >
  struct EntityCommHelper;


  template<>
  struct EntityCommHelper< InteriorBorder_InteriorBorder_Interface >
  {
    static bool send ( const PartitionType p )
    {
      //return (p == InteriorEntity) || (p == BorderEntity);
      return (p == BorderEntity);
    }

    static bool receive ( const PartitionType p )
    {
      //return (p == InteriorEntity) || (p == BorderEntity);
      return (p == BorderEntity);
    }
  };


  template<>
  struct EntityCommHelper< InteriorBorder_All_Interface >
  {
    static bool send ( const PartitionType p )
    {
      return (p == InteriorEntity) || (p == BorderEntity);
    }

    static bool receive ( const PartitionType p )
    {
      //return true;
      return (p != InteriorEntity);
    }
  };


  template<>
  struct EntityCommHelper< Overlap_OverlapFront_Interface >
  {
    static bool send ( const PartitionType p )
    {
      //return (p == InteriorEntity) || (p == BorderEntity) || (p == OverlapEntity);
      return (p != FrontEntity) && (p != GhostEntity);
    }

    static bool receive ( const PartitionType p )
    {
      //return (p == InteriorEntity) || (p == BorderEntity) || (p == OverlapEntity) || (p == FrontEntity);
      return (p != GhostEntity);
    }
  };


  template<>
  struct EntityCommHelper< Overlap_All_Interface >
  {
    static bool send ( const PartitionType p )
    {
      //return (p == InteriorEntity) || (p == BorderEntity) || (p == OverlapEntity);
      return (p != FrontEntity) && (p != GhostEntity);
    }

    static bool receive ( const PartitionType p )
    {
      DUNE_UNUSED_PARAMETER(p);
      return true;
    }
  };


  template<>
  struct EntityCommHelper< All_All_Interface >
  {
    static bool send ( const PartitionType p )
    {
      DUNE_UNUSED_PARAMETER(p);
      return true;
    }

    static bool receive ( const PartitionType p )
    {
      DUNE_UNUSED_PARAMETER(p);
      return true;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_ENTITYCOMMHELPER_HH
