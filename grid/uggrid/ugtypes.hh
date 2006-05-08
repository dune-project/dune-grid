// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGTYPES_HH
#define DUNE_UGTYPES_HH

/** \file
 * \brief Encapsulates a few types from UG
 */

namespace UG {

  namespace D2 {

    struct multigrid;
    struct grid;

    union element;
    struct node;
    struct edge;
    struct vector;

  }

};

namespace UG {

  namespace D3 {

    struct multigrid;
    struct grid;

    union element;
    struct node;
    struct edge;
    struct vector;

  }

};


namespace Dune {

  template <int dim>
  class UGTypes {};

  template <>
  class UGTypes<2>
  {
  public:
    typedef UG::D2::multigrid MultiGridType;

    typedef UG::D2::grid GridType;

    typedef UG::D2::node Node;

    typedef UG::D2::element Element;
  };

  template <>
  class UGTypes<3>
  {
  public:
    typedef UG::D3::multigrid MultiGridType;

    typedef UG::D3::grid GridType;

    typedef UG::D3::node Node;

    typedef UG::D3::element Element;
  };



  /*****************************************************************/
  /*****************************************************************/
  /*****************************************************************/
  /*****************************************************************/

  template <int dim>
  class UGVectorType {};

  template <>
  class UGVectorType<3>
  {
  public:
    typedef UG::D3::vector T;
  };

  template <>
  class UGVectorType<2>
  {
  public:
    typedef UG::D2::vector T;
  };

  template <int codim, int dim>
  class TargetType {};

  template <>
  class TargetType<0,3>
  {
  public:
    typedef UG::D3::element T;
  };

  template <>
  class TargetType<2,3>
  {
  public:
    typedef UG::D3::edge T;
  };

  template <>
  class TargetType<3,3>
  {
  public:
    typedef UG::D3::node T;
  };

  template <>
  class TargetType<0,2>
  {
  public:
    typedef UG::D2::element T;
  };

  template <>
  class TargetType<1,2>
  {
  public:
    typedef UG::D2::edge T;
  };

  template <>
  class TargetType<2,2>
  {
  public:
    typedef UG::D2::node T;
  };

} // end namespace Dune

#endif
