// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include "../mappings.hh"
#include "../conversion.hh"
#include <dune/common/fmatrix.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

using namespace Dune;
using namespace GenericGeometry;

template <class DuneGeometry>
struct DuneCoordTraits {
  enum {dimW = DuneGeometry::coorddimension};  // world dimension
  enum {dimG = DuneGeometry::mydimension};          // grid dimension
  typedef typename DuneGeometry::ctype field_type;
  // general vector and matrix types
  template <int dim>
  struct Vector {
    typedef FieldVector<field_type,dim> Type;
  };
  template <int dimR,int dimC>
  struct Matrix {
    typedef FieldMatrix<field_type,dimR,dimC> Type;
  };
  // Vector of global vectors denoting the edges of the range
  // domain, used to construct a mapping together with an offset.
  // Corners used are
  // p[offset],...,p[offset+Topology::numCorners]
  typedef DuneGeometry coord_vector;
  // mapping is of the form Ax+b (used untested)
  enum {affine = false};
};

namespace Dune {
  template <int dim> class YaspGrid;
  template <int d1,int d2> class AlbertaGrid;
}

template <class Grid>
struct Topology;
template <int d>
struct Topology<YaspGrid<d> > {
  typedef typename Convert< GeometryType :: cube , d >::type Type;
};
template <int d>
struct Topology<AlbertaGrid<d,d> > {
  typedef typename Convert< GeometryType :: simplex , d >::type Type;
};

int phiErr;
int jTinvJTErr;
int volumeErr;

template <class GridViewType>
void test(const GridViewType& view) {
  typedef typename GridViewType::Grid GridType;
  typedef typename Topology<GridType>::Type TopologyType;
  typedef typename GridViewType ::template Codim<0>::Iterator ElementIterator;
  ElementIterator eIt    = view.template begin<0>();
  ElementIterator eEndIt = view.template end<0>();
  typedef typename ElementIterator::Entity::Geometry GeometryType;
  typedef FieldVector<double, GeometryType::mydimension> LocalType;
  typedef FieldVector<double, GeometryType::coorddimension> GlobalType;
  typedef FieldMatrix<double, GeometryType::mydimension, GeometryType::coorddimension>
  JacobianTType;
  typedef FieldMatrix<double, GeometryType::coorddimension, GeometryType::mydimension>
  JacobianType;
  typedef FieldMatrix<double, GeometryType::coorddimension, GeometryType::coorddimension>
  SquareType;

  phiErr = 0;
  jTinvJTErr = 0;
  volumeErr = 0;

  for (; eIt!=eEndIt; ++eIt) {
    const GeometryType& geo = eIt->geometry();
    Mapping<TopologyType,
        DuneCoordTraits<GeometryType> > map(geo);
    LocalType x(0.1);

    // test phi
    GlobalType y=geo.global(x);
    GlobalType yy;
    map.phi(x,yy);
    y-=yy;
    if (y.two_norm2()>1e-10) {
      phiErr++;
      std::cout << "PHI: \n" << y << std::endl;
    }
    // test jacobian
    JacobianType JTInv = geo.jacobianInverseTransposed(x);
    JacobianTType JT;
    map.jacobianT(x,JT);
    SquareType JTInvJT;
    const int m = GeometryType::coorddimension;
    const int n = GeometryType::mydimension;
    for( int i = 0; i < m; ++i ) {
      for( int j = 0; j < m; ++j ) {
        JTInvJT[ i ][ j ] = 0;
        for( int k = 0; k < n; ++k )
          JTInvJT[ i ][ j ] += JTInv[ i ][ k ] * JT[ k ][ j ];
      }
    }
    if (std::abs(JTInvJT[0][0]-1.)+std::abs(JTInvJT[1][1]-1.)+
        std::abs(JTInvJT[1][0])+std::abs(JTInvJT[0][1])>1e-10) {
      jTinvJTErr++;
      std::cout << "JTINVJT: \n" << JTInvJT << std::endl;
      // std::cout << "***\n" << JT << std::endl;
      // std::cout << "***\n" << JTInv << std::endl;
    }
    // test volume
    double vol = geo.volume();
    //  double volvol = map.volume();
  }
}

int main(int argc, char ** argv, char ** envp)
try {
  // this method calls MPI_Init, if MPI is enabled
  // MPIHelper & mpiHelper = MPIHelper::instance(argc,argv);
  // int myrank = mpiHelper.rank();

  if (argc<2) {
    std::cerr << "supply grid file as parameter!" << std::endl;
    return 1;
  }

  // create Grid from DGF parser
  GridPtr<GridType> grid;
  {
    grid = GridPtr<GridType>(argv[1]); // , mpiHelper.getCommunicator() );
  }

  test(grid->leafView());

  if ( phiErr>0) {
    std::cout << phiErr << " errors occured in mapping.phi!" << std::endl;
  }
  else std::cout << "ZERO ERRORS in mapping.phi!!!!!!!" << std::endl;
  if ( jTinvJTErr>0) {
    std::cout << jTinvJTErr << " errors occured in mapping.jacobianT!" << std::endl;
  }
  else std::cout << "ZERO ERRORS in mapping.jacobianT!!!!!!!" << std::endl;

}
catch (Dune::Exception &e) {
  std::cerr << e << std::endl;
  return 1;
}
catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 1;
}
