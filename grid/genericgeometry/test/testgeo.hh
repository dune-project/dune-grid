// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  template <int d1,int d2> class YaspGrid;
  template <int d1,int d2> class AlbertaGrid;
  template< int dim, int dimworld, ALU3dGridElementType elType > class ALU3dGrid;
  template <class ct,int d1,int d2,class CCOMM> class ParallelSimplexGrid;
  template <int dim> class UGGrid;
  template <class Grid>
  struct Topology;
  template <int d>
  struct Topology<YaspGrid<d,d> > {
    typedef typename GenericGeometry::Convert< GeometryType :: cube , d >::type Type;
    static int faceNr(int duneFN) {
      return duneFN;
    }
    enum {correctJacobian=false};
  };
  template <int d>
  struct Topology<UGGrid<d> > {
    typedef typename GenericGeometry::Convert< GeometryType :: cube , d >::type Type;
    static int faceNr(int duneFN) {
      return duneFN;
    }
    enum {correctJacobian=false};
  };
  template <int d>
  struct Topology<AlbertaGrid<d,d> > {
    typedef typename GenericGeometry::Convert< GeometryType :: simplex , d >::type Type;
    static int faceNr(int duneFN) {
      return d-duneFN;
    }
    enum {correctJacobian=false};
  };
  template<>
  struct Topology< ALU3dGrid< 3, 3, tetra > >
  {
    typedef GenericGeometry :: Convert< GeometryType :: simplex, 3 > :: type Type;
    static int faceNr( int duneFN )
    {
      return 3 - duneFN;
    }
    enum {correctJacobian=false};
  };
  template <int d1,int d2,class CCOMM>
  struct Topology<ParallelSimplexGrid<double,d1,d2,CCOMM> > {
    typedef typename GenericGeometry::Convert< GeometryType :: simplex , d1 >::type Type;
    static int faceNr(int duneFN) {
      return d1-duneFN;
    }
    enum {correctJacobian=true};
  };
  namespace GenericGeometry {
    template <>
    struct GeometryTraits<GridType> :
      public DefaultGeometryTraits<
          GridType::ctype,
          GridType::dimension,
          GridType::dimension>
    {};
  }
}
