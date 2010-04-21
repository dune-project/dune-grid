// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/polyallocator.hh>

#include "../cachedmapping.hh"
#include "../conversion.hh"
#include <dune/common/fmatrix.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfoned.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>

#include <dgfgridtype.hh>

#include "testgeo.hh"

using namespace Dune;
using namespace GenericGeometry;

template< class DuneGeometry >
struct MyGeometryTraits
{
  typedef DuneGeometry DuneGeometryType;

  typedef DuneCoordTraits< typename DuneGeometryType :: ctype > CoordTraits;

  static const int dimGrid = DuneGeometryType :: dimension;
  static const int dimWorld = DuneGeometryType :: coorddimension;

  template< class Topology >
  struct Mapping
  {
    typedef CornerMapping< CoordTraits, Topology, dimWorld > type;
  };

  struct Caching
  {
    static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
    static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
    static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
    static const EvaluationType evaluateNormal = ComputeOnDemand;
  };

  typedef Dune::PolyAllocator Allocator;
};


int phiErr;
int jTinvJTErr;
int volumeErr;
int detErr;
int volErr;
int normalErr;

template< class GridViewType >
void test ( const GridViewType &view )
{
  typedef typename GridViewType :: Grid GridType;
  typedef Dune :: Topology< GridType > ConversionType;
  typedef typename ConversionType :: Type TopologyType;
  typedef typename GridViewType ::template Codim<0>::Iterator ElementIterator;
  typedef typename GridViewType :: IntersectionIterator IIteratorType;
  typedef typename GridViewType :: Intersection IntersectionType;
  ElementIterator eIt    = view.template begin<0>();
  ElementIterator eEndIt = view.template end<0>();
  typedef typename ElementIterator::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;
  typedef FieldVector<double, GeometryType::mydimension> LocalType;
  typedef FieldVector<double, GeometryType::mydimension-1> LocalFaceType;
  typedef FieldVector<double, GeometryType::coorddimension> GlobalType;
  typedef FieldMatrix<double, GeometryType::mydimension, GeometryType::coorddimension>
  JacobianTType;
  typedef FieldMatrix<double, GeometryType::coorddimension, GeometryType::mydimension>
  JacobianType;
  typedef FieldMatrix<double, GeometryType::mydimension, GeometryType::mydimension>
  SquareType;

  phiErr = 0;
  jTinvJTErr = 0;
  volumeErr = 0;
  detErr = 0;
  volErr = 0;

  for (; eIt!=eEndIt; ++eIt) {
    const GeometryType& geo = eIt->geometry();
    std::vector< Dune::FieldVector<double, GeometryType::coorddimension> > corners;
    for (int i=0; i<geo.corners(); i++)
      corners.push_back( geo.corner(i) );
    CachedMapping< TopologyType, MyGeometryTraits< GeometryType > > map( corners );
    LocalType x(0.1);
    for (int i=0; i<10000; ++i) {
      // test phi
      GlobalType y = geo.global(x);
      GlobalType yy = map.global( x );
      y-=yy;
      if (y.two_norm2()>1e-10) {
        phiErr++;
        std::cout << "PHI: \n" << y << std::endl;
      }
      // test jacobian
      JacobianType JTInv = geo.jacobianInverseTransposed(x);
      JacobianTType JT = map.jacobianTransposed( x );
      JacobianType JTInvM = map.jacobianInverseTransposed( x );
      double det = geo.integrationElement(x);
      double detM = map.integrationElement(x);
      if (std::abs(det-detM)>1e-10) {
        std::cout << " Error in det: G = " << det << " M = " << detM << std::endl;
        detErr++;
      }
      SquareType JTInvJT;
      const int n = GeometryType::coorddimension;
      const int m = GeometryType::mydimension;
      for( int i = 0; i < m; ++i ) {
        for( int j = 0; j < m; ++j ) {
          JTInvJT[ i ][ j ] = 0;
          for( int k = 0; k < n; ++k )
            JTInvJT[ i ][ j ] += JT[i][k] * JTInvM[k][j];
        }
      }
      if (std::abs(JTInvJT[0][0]-1.)+std::abs(JTInvJT[1][1]-1.)+
          std::abs(JTInvJT[1][0])+std::abs(JTInvJT[0][1])>1e-10) {
        jTinvJTErr++;
        std::cout << "JTINVJT: \n" << JTInvJT << std::endl;
        // std::cout << "***\n" << JT << std::endl;
        // std::cout << "***\n" << JTInv << std::endl;
      }
      JTInv -= JTInvM;
      if (JTInv.frobenius_norm2()>1e-10) {
        jTinvJTErr++;
        std::cout << "JTINVJT: \n" << JTInv << std::endl;
      }
      // test volume
      double vol  = geo.volume();
      double volM = map.volume();
      if (std::abs(vol-volM)>1e-10) {
        volErr++;
        std::cout << "volume : " << vol << " " << volM << std::endl;
      }
      // test normal
      IIteratorType iiter = view.ibegin(*eIt);
      for (; iiter != view.iend(*eIt); ++ iiter) {
        LocalFaceType xf(0.1);
        GlobalType n = iiter->integrationOuterNormal(xf);
        LocalType xx(iiter->geometryInInside().global(xf));
        const int genericFaceNr = iiter->indexInInside();
        GlobalType nM = map.normal( genericFaceNr, xx );
        if ((n-nM).two_norm2()>1e-10) {
          normalErr++;
          std::cout << nM.two_norm() << " " << n.two_norm() << std::endl;
          std::cout << "normal: "
                    << " ( " << nM << " )   ( " << n << " )   "
                    << n-nM
                    << std::endl;
        }
      }
    }
  }
}


int main ( int argc, char **argv )
try
{
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper :: instance( argc, argv );
  // MPIHelper &mpiHelper = MPIHelper :: instance( argc, argv );
  // int myrank = mpiHelper.rank();

  if (argc<2) {
    std::cerr << "supply grid file as parameter!" << std::endl;
    return 1;
  }

  // create Grid from DGF parser
  GridPtr< GridSelector::GridType > grid( argv[ 1 ] );

  test(grid->leafView());

  if ( phiErr>0) {
    std::cout << phiErr << " errors occured in mapping.phi?" << std::endl;
  }
  else std::cout << "ZERO ERRORS in mapping.phi!!!!!!!" << std::endl;
  if ( detErr>0) {
    std::cout << detErr << " errors occured in mapping.det?" << std::endl;
  }
  else std::cout << "ZERO ERRORS in mapping.det!!!!!!!" << std::endl;
  if ( jTinvJTErr>0) {
    std::cout << jTinvJTErr
              << " errors occured in mapping.jacobianT?" << std::endl;
  }
  else std::cout << "ZERO ERRORS in mapping.jacobianT!!!!!!!" << std::endl;
  if ( volErr>0) {
    std::cout << volErr
              << " errors occured in mapping.volume?" << std::endl;
  }
  else std::cout << "ZERO ERRORS in mapping.volume!!!!!!!" << std::endl;
  if ( normalErr>0) {
    std::cout << normalErr
              << " errors occured in mapping.normal?" << std::endl;
  }
  else std::cout << "ZERO ERRORS in mapping.normal!!!!!!!" << std::endl;
}
catch (Dune::Exception &e) {
  std::cerr << e << std::endl;
  return 1;
}
catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 1;
}
