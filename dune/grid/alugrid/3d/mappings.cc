// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"
#include "mappings.hh"

namespace Dune
{

  const double TrilinearMapping :: _epsilon = 1.0e-8 ;
  const double BilinearSurfaceMapping :: _epsilon = 1.0e-8 ;
  const double SurfaceNormalCalculator :: _epsilon = 1.0e-8 ;



  // Implementation of NonConformingFaceMapping
  // ------------------------------------------

  template< class Comm >
  void NonConformingFaceMapping< tetra, Comm >
  ::child2parent ( const CoordinateType &childCoordinates,
                   CoordinateType& parentCoordinates ) const
  {
    if (rule_ == RefinementRuleType::nosplit) {
      child2parentNosplit(childCoordinates, parentCoordinates);
    }
    else if (rule_ == RefinementRuleType::e01) {
      child2parentE01(childCoordinates, parentCoordinates);
    }
    else if (rule_ == RefinementRuleType::e12) {
      child2parentE12(childCoordinates, parentCoordinates);
    }
    else if (rule_ == RefinementRuleType::e20) {
      child2parentE20(childCoordinates, parentCoordinates);
    }
    else if (rule_ == RefinementRuleType::iso4) {
      child2parentIso4(childCoordinates, parentCoordinates);
    }
    else {
      // check with cases in Hface3Rule (gitter_sti.h)
      assert(false);
      DUNE_THROW(InvalidStateException,"Invalid refinement rule");
    } // end if
  }


  template< class Comm >
  typename NonConformingFaceMapping< tetra, Comm >::CoordinateType
  NonConformingFaceMapping< tetra, Comm >
  ::child2parent ( const FieldVector< alu3d_ctype, 2 > &childCoordinates ) const
  {
    CoordinateType childBary, parentBary;
    childBary[ 0 ] = 1.0 - childCoordinates[ 0 ] - childCoordinates[ 1 ];
    childBary[ 1 ] = childCoordinates[ 0 ];
    childBary[ 2 ] = childCoordinates[ 1 ];
    child2parent( childBary, parentBary );
    return parentBary;
  }


  template< class Comm >
  void NonConformingFaceMapping< tetra, Comm >
  ::child2parentNosplit ( const CoordinateType& childCoordinates,
                          CoordinateType& parentCoordinates ) const
  {
    parentCoordinates = childCoordinates;
  }


  template< class Comm >
  void NonConformingFaceMapping< tetra, Comm >
  ::child2parentE01 ( const CoordinateType &childCoordinates,
                      CoordinateType &parentCoordinates ) const
  {
    DUNE_THROW(NotImplemented,
               "This refinement rule is currently not supported");
  }


  template< class Comm >
  void NonConformingFaceMapping< tetra, Comm >
  ::child2parentE12 ( const CoordinateType &childCoordinates,
                      CoordinateType &parentCoordinates ) const
  {
    DUNE_THROW(NotImplemented,
               "This refinement rule is currently not supported");
  }


  template< class Comm >
  void NonConformingFaceMapping< tetra, Comm >
  ::child2parentE20 ( const CoordinateType &childCoordinates,
                      CoordinateType &parentCoordinates ) const
  {
    DUNE_THROW(NotImplemented,
               "This refinement rule is currently not supported");
  }


  template< class Comm >
  void NonConformingFaceMapping< tetra, Comm >
  ::child2parentIso4 ( const CoordinateType &childCoordinates,
                       CoordinateType &parentCoordinates ) const
  {
    /*
       // The ordering of the coordinates are according to a Dune reference triangle
       //
       //  NOTE: all coordinates are barycentric (with respect to (P_0, P_1, P_2)
       //
       //                  P_2 = (0,0,1)
       //                   |\
       //                   | \
       //                   |  \      each sub triangle is numbered as used below
       //                   |   \     local numbering is count clockwise
       //                   |    \    starting with the lower left vertex
       //                   |     \   (i.e. child 0 consits of  { P_0, (P_0+P_1)/2 , (P_0+P_2)/2 }  )
       //                   |      \
       //                   |   1   \
       //                   |        \
       //                   |         \
       //      (0.5,0,0.5)  |----------\  (0,0.5,0.5) = (P_1 + P_2)/2
       //   = (P_0 + P_2)/2 |\         |\
       //                   | \        | \
       //                   |  \   3   |  \
       //                   |   \      |   \
       //                   |    \     |    \
       //                   |     \    |     \
       //                   |      \   |      \
       //                   |  0    \  |   2   \
       //                   |        \ |        \
       //                   |         \|         \
       //                   -----------------------
       //         (1,0,0) = P_0   (0.5,0.5,0)    P_1 = (0,1,0)
       //                         = (P_0 + P_1)/
       //
       //  NOTE: the strange numbering of the childs is due to the swap
       //  of the vertex number form ALUGrid to Dune reference triangles faces
       //  This means that in ALUGrid child 1 and 2 are swaped compared to
       //  this example here.
     */

    // this mapping map from the points (P_0,P_1,P_2) to the
    // 4 sub triangles from the picture above
    //
    // TODO: this mapping is static, so store it in an map
    switch(nChild_) {
    case 0 :
      // (1,0,0) --> (1,0,0)
      // (0,1,0) --> (0.5,0,5,0)
      // (0,0,1) --> (0.5,0,0.5)
      parentCoordinates[0] =
        1.0 - 0.5*childCoordinates[1] - 0.5*childCoordinates[2];

      // this rocks , best bug ever
      //parentCoordinates[1] = 0,5*childCoordinates[1];

      parentCoordinates[1] = 0.5*childCoordinates[1];
      parentCoordinates[2] = 0.5*childCoordinates[2];
      break;
    case 1 : // swaped case 1 and case 2
      // (1,0,0) --> (0.5,0,0.5)
      // (0,1,0) --> (0,0,5,0)
      // (0,0,1) --> (0,0,1)
      parentCoordinates[0] = 0.5*childCoordinates[0];
      parentCoordinates[1] = 0.5*childCoordinates[1];
      parentCoordinates[2] =
        1.0 - 0.5*childCoordinates[0] - 0.5*childCoordinates[1];
      break;
    case 2 :
      // (1,0,0) --> (0.5,0,5,0)
      // (0,1,0) --> (0,1,0)
      // (0,0,1) --> (0.5,0.5,0)
      parentCoordinates[0] = 0.5*childCoordinates[0];
      parentCoordinates[1] =
        1.0 - 0.5*childCoordinates[0] - 0.5*childCoordinates[2];
      parentCoordinates[2] = 0.5*childCoordinates[2];
      break;
    case 3 :
      // (1,0,0) --> (0.5,0,0.5)
      // (0,1,0) --> (0.5,0.5,0)
      // (0,0,1) --> (0,0.5,0.5)
      // here swaped all to the next position
      parentCoordinates[1] = 0.5 - 0.5*childCoordinates[0];
      parentCoordinates[2] = 0.5 - 0.5*childCoordinates[1];
      parentCoordinates[0] = 0.5 - 0.5*childCoordinates[2];
      break;
    default :
      DUNE_THROW(RangeError, "Only 4 children on a tetrahedron face (val = "
                 << nChild_ << ")");

    } // end switch
  }


  template< class Comm >
  void NonConformingFaceMapping< hexa, Comm >
  ::child2parent ( const CoordinateType &childCoordinates,
                   CoordinateType& parentCoordinates ) const
  {
    if (rule_ == RefinementRuleType::nosplit) {
      child2parentNosplit(childCoordinates, parentCoordinates);
    }
    else if (rule_ == RefinementRuleType::iso4) {
      child2parentIso4(childCoordinates, parentCoordinates);
    }
    else {
      // check with cases in Hface3Rule (gitter_sti.h)
      assert(false);
      DUNE_THROW(InvalidStateException,"Invalid refinement rule");
    }
  }


  template< class Comm >
  typename NonConformingFaceMapping< hexa, Comm >::CoordinateType
  NonConformingFaceMapping< hexa, Comm >
  ::child2parent ( const FieldVector< alu3d_ctype, 2 > &childCoordinates ) const
  {
    CoordinateType parentCoordinates;
    child2parent( childCoordinates, parentCoordinates );
    return parentCoordinates;
  }




  template< class Comm >
  void NonConformingFaceMapping< hexa, Comm >
  ::child2parentNosplit ( const CoordinateType& childCoordinates,
                          CoordinateType& parentCoordinates ) const
  {
    parentCoordinates = childCoordinates;
  }


  template< class Comm >
  void NonConformingFaceMapping< hexa, Comm >
  ::child2parentIso4 ( const CoordinateType &childCoordinates,
                       CoordinateType &parentCoordinates ) const
  {
    // The ordering of the coordinates are according to a Dune reference elemen
    //
    //
    //   (0,1)                   (1,1)
    //    -------------------------
    //    |           |           |     childs within the reference
    //    |           |           |     quadrilateral of Dune
    //    |    1      |     2     |
    //    |           |           |
    //    |           |           |
    //    |-----------|-----------|
    //    |           |           |
    //    |           |           |
    //    |    0      |     3     |
    //    |           |           |
    //    |           |           |
    //    -------------------------
    //  (0,0)                    (1,0)
    //
    //
    switch(nChild_) {
    case 0 :
      parentCoordinates[0] = 0.5*childCoordinates[0];
      parentCoordinates[1] = 0.5*childCoordinates[1];
      break;
    case 1 :
      parentCoordinates[0] = 0.5*childCoordinates[0];
      parentCoordinates[1] = 0.5*childCoordinates[1] + 0.5;
      break;
    case 2 :
      parentCoordinates[0] = 0.5*childCoordinates[0] + 0.5;
      parentCoordinates[1] = 0.5*childCoordinates[1] + 0.5;
      break;
    case 3 :
      parentCoordinates[0] = 0.5*childCoordinates[0] + 0.5;
      parentCoordinates[1] = 0.5*childCoordinates[1];
      break;
    default :
      DUNE_THROW(RangeError, "Only 4 children on a hexahedron face (val = "
                 << nChild_ << ")");
    } // end switch
  }



  // Explicit Template Instatiation
  // ------------------------------

  template class NonConformingFaceMapping< tetra, No_Comm >;
  template class NonConformingFaceMapping< hexa, No_Comm >;

#if ALU3DGRID_PARALLEL
  template class NonConformingFaceMapping< tetra, MPI_Comm >;
  template class NonConformingFaceMapping< hexa, MPI_Comm >;
#endif // #if ALU3DGRID_PARALLEL

} // end namespace Dune
