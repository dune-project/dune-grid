// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_INCLUDE_HH
#define DUNE_ALU2DGRID_INCLUDE_HH

#include <alugrid_2d.h>
#define ALU2DSPACE ALU2dGridSpace ::

namespace Dune {

  struct ALU2dImplTraits {
    template <int cdim>
    struct Codim;
  };

  template<>
  struct ALU2dImplTraits::Codim<0> {
    typedef ALU2DSPACE Hmesh_basic::helement_t InterfaceType;
  };

  template<>
  struct ALU2dImplTraits::Codim<1> {
    typedef ALU2DSPACE Hmesh_basic::helement_t InterfaceType;
  };

  template <>
  struct ALU2dImplTraits::Codim<2> {
    typedef ALU2DSPACE Vertex InterfaceType;
  };

} //end namespace Dune
#endif
