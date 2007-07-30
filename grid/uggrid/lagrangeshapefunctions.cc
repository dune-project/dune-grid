// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"
#include "lagrangeshapefunctions.hh"

namespace Dune {

  namespace UGShapeFunctions {

    template<typename T, int d>
    LagrangeShapeFunctionSetContainer<T,d> LagrangeShapeFunctions<T,d>::general;

    namespace {

      template <class T, int d>
      struct InitLagrangeShapefunctions
      {
        LagrangeShapeFunctionSetContainer<T,d> & f7;
        InitLagrangeShapefunctions() :
          f7(LagrangeShapeFunctions<T,d>::general)
        {
          InitLagrangeShapefunctions<T,d-1> i;
        };
      };

      template <class T>
      struct InitLagrangeShapefunctions<T,0>
      {
        enum { d=0 };
        InitLagrangeShapefunctions()
        {};
      };

      // force creation of symbols and code ...
      void init_lagrangeshapefunctions()
      {
        InitLagrangeShapefunctions<double, 3> i1;
      }
    }

  }
}
