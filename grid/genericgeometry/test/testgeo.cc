// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
// #define DUNE_THROW(E, m) assert(0)

#include "../geometry.hh"
#include "../conversion.hh"
#include <dune/common/fmatrix.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/psg/dgfgridtype.hh>
// #include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

using namespace Dune;

template <class DuneGeometry>
struct DuneCoordTraits {
  enum {dimCoord = DuneGeometry::coorddimension};  // world dimension
  enum {dimGrid  = DuneGeometry::dimension};       // grid dimension
  typedef typename DuneGeometry::ctype FieldType;
  // general vector and matrix types
  template <int dim>
  struct Vector {
    typedef FieldVector<FieldType,dim> Type;
  };
  template <int dimR,int dimC>
  struct Matrix {
    typedef FieldMatrix<FieldType,dimR,dimC> Type;
  };
  // Vector of global vectors denoting the edges of the range
  // domain, used to construct a mapping together with an offset.
  // Corners used are
  // p[offset],...,p[offset+Topology::numCorners]
  typedef DuneGeometry CoordVector;
  // mapping is of the form Ax+b (used untested)
  enum {affine = false};
};

template <class Traits>
struct DuneCache : public GenericGeometry::ComputeAll<Traits> {
  typedef typename Traits::CoordVector GeometryType;
  DuneCache(const GeometryType&) {}
};

template <class CoordTraits>
struct DuneCache< GenericGeometry::MappingTraits<CoordTraits::CoordVector::mydimension,
        CoordTraits> > {
  enum {jTCompute = GenericGeometry::geoCompute,
        jTInvCompute = GenericGeometry::geoCompute,
        intElCompute = GenericGeometry::geoCompute,
        normalCompute = GenericGeometry::geoCompute};
  typedef typename CoordTraits::CoordVector GeometryType;
  typedef GenericGeometry::MappingTraits<CoordTraits::CoordVector::mydimension,
      CoordTraits> Traits;
  const GeometryType& geo_;
  DuneCache(const GeometryType& geo) : geo_(geo) {}
  void jacobianT(typename Traits::JacobianTransposeType& d) const {}
  void integrationElement(typename Traits::FieldType& intEl) const {
    FieldVector<double,GeometryType::mydimension> x(0);
    intEl = geo_.integrationElement(x);
  }
  void jacobianInverseTransposed(typename Traits::JacobianType& dInv) const {}
  void normal(int face, typename Traits::GlobalCoordType& n) const {}
};


namespace Dune {
  template <int d1,int d2> class YaspGrid;
  template <int d1,int d2> class AlbertaGrid;
  template <class ct,int d1,int d2,class CCOMM> class ParallelSimplexGrid;
  template <int dim> class UGGrid;
}

template <class Grid>
struct Topology;
template <int d>
struct Topology<YaspGrid<d,d> > {
  typedef typename GenericGeometry::Convert< GeometryType :: cube , d >::type Type;
  static int faceNr(int duneFN) {
    return duneFN;
  }
};
template <int d>
struct Topology<UGGrid<d> > {
  typedef typename GenericGeometry::Convert< GeometryType :: cube , d >::type Type;
  static int faceNr(int duneFN) {
    return duneFN;
  }
};
template <int d>
struct Topology<AlbertaGrid<d,d> > {
  typedef typename GenericGeometry::Convert< GeometryType :: simplex , d >::type Type;
  static int faceNr(int duneFN) {
    return d-duneFN;
  }
};
template <int d1,int d2,class CCOMM>
struct Topology<ParallelSimplexGrid<double,d1,d2,CCOMM> > {
  typedef typename GenericGeometry::Convert< GeometryType :: simplex , d1 >::type Type;
  static int faceNr(int duneFN) {
    return d1-duneFN;
  }
};

int phiErr;
int jTinvJTErr;
int volumeErr;
int detErr;
int volErr;
int normalErr;

template <class GridViewType>
void test(const GridViewType& view) {
  typedef typename GridViewType::Grid GridType;
  typedef Topology<GridType> ConversionType;
  typedef typename ConversionType::Type TopologyType;
  typedef typename GridViewType ::template Codim<0>::Iterator ElementIterator;
  typedef typename GridViewType :: IntersectionIterator IIteratorType;
  typedef typename GridViewType :: Intersection IntersectionType;
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
  typedef GenericGeometry::Geometry< TopologyType,
      DuneCoordTraits<GeometryType>,
      DuneCache > GenericGeometryType;

  phiErr = 0;
  jTinvJTErr = 0;
  volumeErr = 0;
  detErr = 0;
  volErr = 0;

  for (int j=0; j<1; ++j) {
    ElementIterator eEndIt = view.template end<0>();
    ElementIterator eIt    = view.template begin<0>();
    for (; eIt!=eEndIt; ++eIt) {
      const GeometryType& geoDune = eIt->geometry();
      GenericGeometryType mapDune(geoDune,typename GenericGeometryType::CachingType(geoDune) );
      {
        LocalType x(0.1);
        GlobalType y=geoDune.global(x);
        GlobalType yM=mapDune.global(x);
      }

      for (int i=0; i<1; ++i) {
        // const GenericGeometryType& geo = mapDune;        // 17.4 (24)
        const GenericGeometryType& map = mapDune;
        const GeometryType& geo = geoDune;            // 15.2
        // const GeometryType& map = geoDune;
        LocalType x(0.1);
        // test phi
        GlobalType y=geo.global(x);
        GlobalType yM=map.global(x);
        if ((y-yM).two_norm2()>1e-10) {
          phiErr++;
          std::cout << "Error in PHI: G = " << y << " M = " << yM << " " << std::endl;
        }

        const JacobianType&  JTInv  = geo.jacobianInverseTransposed(x);
        const JacobianType&  JTInvM = map.jacobianInverseTransposed(x);
        JacobianType testJTInv(JTInv);
        testJTInv -= JTInvM;
        if (testJTInv.frobenius_norm2()>1e-10) {
          jTinvJTErr++;
          std::cout << "JTINVJT: \n" << JTInv << std::endl;
        }
        {
          // test jacobian
          const JacobianTType& JTM    = map.jacobianT(x);
          SquareType JTInvJT;
          const int n = GeometryType::coorddimension;
          const int m = GeometryType::mydimension;
          for( int i = 0; i < m; ++i ) {
            for( int j = 0; j < m; ++j ) {
              JTInvJT[ i ][ j ] = 0;
              for( int k = 0; k < n; ++k )
                JTInvJT[ i ][ j ] += JTM[i][k] * JTInvM[k][j];
            }
          }
          if (std::abs(JTInvJT[0][0]-1.)+std::abs(JTInvJT[1][1]-1.)+
              std::abs(JTInvJT[1][0])+std::abs(JTInvJT[0][1])>1e-10) {
            jTinvJTErr++;
            std::cout << "JTINVJT: \n" << JTInvJT << std::endl;
            // std::cout << "***\n" << JT << std::endl;
            // std::cout << "***\n" << JTInv << std::endl;
          }
        }
        // test integration element
        double det = geo.integrationElement(x);
        double detM = map.integrationElement(x);
        if (std::abs(det-detM)>1e-10) {
          std::cout << " Error in det: G = " << det << " M = " << detM << std::endl;
          detErr++;
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
          LocalType xx(iiter->intersectionSelfLocal().global(xf));
          const GlobalType& n  = iiter->integrationOuterNormal(xf);
          const GlobalType& nM = map.normal(ConversionType::faceNr(iiter->numberInSelf()),xx);
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
  GridPtr<GridType> grid( argv[ 1 ] );

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
/*
   catch( PSG :: Exception &e ) {
   std :: cerr << e << std :: endl;
   }
 */
catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 1;
}
