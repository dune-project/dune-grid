// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/grid/alugrid/3d/topology.hh>

namespace Dune
{

  template <int d1,int d2> class YaspGrid;
  template <int d1,int d2> class AlbertaGrid;
  template< int dim, int dimworld > class ALUSimplexGrid;
  template< int dim, int dimworld, ALU3dGridElementType elType > class ALU3dGrid;
  template <int dim> class UGGrid;



  // Topology
  // --------

  template <class Grid>
  struct Topology;

  template< int d >
  struct Topology< YaspGrid< d, d > >
  {
    enum { basicType = GeometryType :: cube };
    enum { dimension = d };
    enum {correctJacobian=false};

    typedef GenericGeometry :: Convert< (GeometryType :: BasicType) basicType, dimension > Convert;
    typedef typename Convert :: type Type;
  };

  template <int d>
  struct Topology<UGGrid<d> >
  {
    enum { basicType = GeometryType :: cube };
    enum { dimension = d };
    enum {correctJacobian=false};

    typedef GenericGeometry :: Convert< (GeometryType :: BasicType) basicType, dimension > Convert;
    typedef typename Convert :: type Type;
  };

  template <int d>
  struct Topology<AlbertaGrid<d,d> >
  {
    enum { basicType = GeometryType :: simplex };
    enum { dimension = d };
    enum {correctJacobian=false};

    typedef GenericGeometry :: Convert< (GeometryType :: BasicType) basicType, dimension > Convert;
    typedef typename Convert :: type Type;
  };

  template<>
  struct Topology< ALU3dGrid< 3, 3, tetra > >
  {
    enum { basicType = GeometryType :: simplex };
    enum { dimension = 3 };
    enum {correctJacobian=false};

    typedef GenericGeometry :: Convert< (GeometryType :: BasicType) basicType, dimension > Convert;
    typedef Convert :: type Type;
  };

  template<>
  struct Topology< ALUSimplexGrid< 3, 3 > >
  {
    enum { basicType = GeometryType :: simplex };
    enum { dimension = 3 };
    enum {correctJacobian=true};

    typedef GenericGeometry :: Convert< (GeometryType :: BasicType) basicType, dimension > Convert;
    typedef Convert :: type Type;
  };

  template<>
  struct Topology< ALUSimplexGrid< 2, 2 > >
  {
    enum { basicType = GeometryType :: simplex };
    enum { dimension = 2 };
    enum {correctJacobian=true};

    typedef GenericGeometry :: Convert< (GeometryType :: BasicType) basicType, dimension > Convert;
    typedef Convert :: type Type;
  };



  // GeometryTraits
  // --------------

  namespace GenericGeometry
  {

    template< class Grid >
    struct GeometryTraits
      : public DefaultGeometryTraits
        < typename Grid :: ctype, Grid :: dimension, Grid :: dimensionworld >
    {
      enum { hybrid = true };
      //enum { duneType = Topology< Grid > :: basicType };
    };

  }

}
