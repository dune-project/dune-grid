// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_GRIDENUMS_HH
#define DUNE_GRID_COMMON_GRIDENUMS_HH

#include <iostream>

#include <dune/common/exceptions.hh>

namespace Dune {

  /**
   * \addtogroup gridpartitions Parallel Grid Partitions
   * \{
   */

  /** \brief Attributes used in the generic overlap model
   *
     \code
     #include <dune/grid/common/gridenums.hh>
     \endcode
   *
   * The values are ordered intentionally in order to be able to
   * define ranges of partition types.
   *
   * @ingroup GIRelatedTypes
   */
  enum PartitionType {
    InteriorEntity=0,     //!< all interior entities
    BorderEntity=1  ,     //!< on boundary between interior and overlap
    OverlapEntity=2 ,     //!< all entities lying in the overlap zone
    FrontEntity=3  ,      //!< on boundary between overlap and ghost
    GhostEntity=4         //!< ghost entities
  };

  /** \brief Provide names for the partition types
   *
     \code
     #include <dune/grid/common/gridenums.hh>
     \endcode
   *
   * @ingroup GIRelatedTypes
   */
  inline std::string PartitionName(PartitionType type)
  {
    switch(type) {
    case InteriorEntity :
      return "interior";
    case BorderEntity :
      return "border";
    case OverlapEntity :
      return "overlap";
    case FrontEntity :
      return "front";
    case GhostEntity :
      return "ghost";
    default :
      DUNE_THROW(NotImplemented, "name of unknown partition type requested");
    }
  }

  //! write a PartitionType to a stream
  /**
     \code
     #include <dune/grid/common/gridenums.hh>
     \endcode
   *
   * @ingroup GIRelatedTypes
   */
  inline std::ostream &operator<< ( std::ostream &out, const PartitionType &type )
  {
    return out << PartitionName( type );
  }


  /** \brief Parameter to be used for the communication functions
   *
     \code
     #include <dune/grid/common/gridenums.hh>
     \endcode
   *
   * @ingroup GIRelatedTypes
   */
  enum InterfaceType {
    InteriorBorder_InteriorBorder_Interface=0,     //!< send/receive interior and border entities
    InteriorBorder_All_Interface=1,                //!< send interior and border, receive all entities
    Overlap_OverlapFront_Interface=2,              //!< send overlap, receive overlap and front entities
    Overlap_All_Interface=3,                       //!< send overlap, receive all entities
    All_All_Interface=4                            //!< send all and receive all entities
  };


  //! write an InterfaceType to a stream
  /**
     \code
     #include <dune/grid/common/gridenums.hh>
     \endcode
   *
   * @ingroup GIRelatedTypes
   */
  inline std::ostream &operator<< ( std::ostream &out, const InterfaceType &type )
  {
    switch( type )
    {
    case InteriorBorder_InteriorBorder_Interface :
      return out << "interior-border / interior-border interface";

    case InteriorBorder_All_Interface :
      return out << "interior-border / all interface";

    case Overlap_OverlapFront_Interface :
      return out << "overlap / overlap-front interface";

    case Overlap_All_Interface :
      return out << "overlap / all interface";

    case All_All_Interface :
      return out << "all / all interface";

    default :
      return out << "unknown interface";
    }
  }


  /** \brief Parameter to be used for the parallel level- and leaf iterators
   *
     \code
     #include <dune/grid/common/gridenums.hh>
     \endcode
   *
   * @ingroup GIRelatedTypes
   */
  enum PartitionIteratorType {
    Interior_Partition=0,           //!< only interior entities
    InteriorBorder_Partition=1,     //!< interior and border entities
    Overlap_Partition=2,            //!< interior, border, and overlap entities
    OverlapFront_Partition=3,       //!< interior, border, overlap and front entities
    All_Partition=4,                //!< all entities
    Ghost_Partition=5               //!< only ghost entities
  };


  //! write a PartitionIteratorType to a stream
  /**
     \code
     #include <dune/grid/common/gridenums.hh>
     \endcode
   *
   * @ingroup GIRelatedTypes
   */
  inline std::ostream &operator<< ( std::ostream &out, const PartitionIteratorType &type )
  {
    static std::string name[ 6 ] = { "interior partition", "interior-border partition", "overlap partition",
                                     "overlap-front partition", "all partition", "ghost partition" };
    return out << name[ type ];
  }


  /** \brief Define a type for communication direction parameter
   *
     \code
     #include <dune/grid/common/gridenums.hh>
     \endcode
   *
   * @ingroup GIRelatedTypes
   */
  enum CommunicationDirection {
    ForwardCommunication,     //!< communicate as given in InterfaceType
    BackwardCommunication     //!< reverse communication direction
  };

  /**
   * \}
   */


}
#endif
