// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDENUMS_HH
#define DUNE_GRIDENUMS_HH

#include <iostream>

#include <dune/common/exceptions.hh>

namespace Dune {


  /** \brief Attributes used in the generic overlap model

     The values are ordered intentionally in order to be able to
     define ranges of partition types.

     @ingroup GIRelatedTypes
   */
  enum PartitionType
  {
    InteriorEntity = (1 << 0), //!< all interior entities
    BorderEntity = (1 << 1),   //!< on boundary between interior and overlap
    OverlapEntity = (1 << 2),  //!< all entities lying in the overlap zone
    FrontEntity = (1 << 3),    //!< on boundary between overlap and ghost
    GhostEntity = (1 << 4)     //!< ghost entities
  };

  /** \brief Provide names for the partition types
     @ingroup GIRelatedTypes
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


  inline std::ostream &operator<< ( std::ostream &out, const PartitionType &type )
  {
    return out << PartitionName( type );
  }


  /** \brief Parameter to be used for the communication functions
     @ingroup GIRelatedTypes
   */
  enum InterfaceType {
    InteriorBorder_InteriorBorder_Interface=0,     //!< send/receive interior and border entities
    InteriorBorder_All_Interface=1,                //!< send interior and border, receive all entities
    Overlap_OverlapFront_Interface=2,              //!< send overlap, receive overlap and front entities
    Overlap_All_Interface=3,                       //!< send overlap, receive all entities
    All_All_Interface=4                            //!< send all and receive all entities
  };

  /** \brief Parameter to be used for the parallel level- and leaf iterators
     @ingroup GIRelatedTypes
   */
  enum PartitionIteratorType
  {
    //! only interior entities
    Interior_Partition = InteriorEntity,
    //! interior and border entities
    InteriorBorder_Partition = InteriorEntity | BorderEntity,
    //! interior, border, and overlap entities
    Overlap_Partition = InteriorBorder_Partition | OverlapEntity,
    //! interior, border, overlap, and front entities
    OverlapFront_Partition = Overlap_Partition | FrontEntity,
    //! all entities
    All_Partition = OverlapFront_Partition | GhostEntity,
    //! only ghost entities
    Ghost_Partition = GhostEntity
  };


  inline std::ostream &operator<< ( std::ostream &out, const PartitionIteratorType &type )
  {
    static std::string name[ 6 ] = { "interior partition", "interior-border partition", "overlap partition",
                                     "overlap-front partition", "all partition", "ghost partition" };
    return out << name[ type ];
  }


  /** \brief Define a type for communication direction parameter
     @ingroup GIRelatedTypes
   */
  enum CommunicationDirection {
    ForwardCommunication,         //!< communicate as given in InterfaceType
    BackwardCommunication         //!< reverse communication direction
  };

}
#endif
