// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
    \brief A unit test for the BasicGeometry class
 */


#include <config.h>

#include <cstddef>
#include <cstdlib>
#include <ostream>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/ios_state.hh>

#include <dune/grid/genericgeometry/geometry.hh>



using namespace Dune;

void fail(int &result) {
  result = 1;
}
void pass(int &result) {
  if(result == 77) result = 0;
}

/** \brief Test the interface of the given BasicGeometry object
 *
 * \param geometry       The Geometry object to test.
 * \param expectedVolume The volume to expect from the geometry.
 * \param expectedAffine Whether the geometry should be affine.
 * \param result         Collect pass/fail results.
 */
template <class TestGeometry>
void testBasicGeometry(const TestGeometry& geometry,
                       typename TestGeometry::ctype expectedVolume,
                       bool expectedAffine, int &result)
{
  typedef typename TestGeometry::ctype ctype;
  const int dim = TestGeometry::mydimension;

  if(dim > 0) {
    geometry.normal(0, FieldVector<double,dim>(0));
    pass(result);
  }

  if(expectedVolume == expectedVolume) {
    ctype volume = geometry.volume();
    if(std::abs(volume - expectedVolume) > 1e-8) {
      std::cerr << "Volume: " << volume << ", but "
                << expectedVolume << " was expected!" << std::endl;
      fail(result);
    }
    else
      pass(result);
  }
  else
    std::cerr << "Warning: volume check skipped." << std::endl;

  bool affine = geometry.affine();
  if(affine != expectedAffine) {
    Dune::ios_base_all_saver saver(std::cerr);
    std::cerr << std::boolalpha;
    std::cerr << "Affine: \"" << affine << "\", but "
              << "\"" << expectedAffine << "\" was expected!" << std::endl;
    fail(result);
  }
  else
    pass(result);
}

int main (int argc , char **argv) try
{
  // 77 means "SKIP"
  int result = 77;

  double volume;
  bool affine;
  GeometryType gt;

  ////////////////////////////////////////////////////////////////////////
  //
  //  coorddim = 0
  //

  {
    static const std::size_t coorddim = 0;
    std::cout << "== coorddimension = " << coorddim << std::endl;

    std::vector<FieldVector<double, coorddim> > corners;

    //
    //  mydim = 0
    //

    {
      static const std::size_t mydim = 0;
      std::cout << "=== mydimension = " << mydim << std::endl;

      typedef GenericGeometry::BasicGeometry<
          mydim, GenericGeometry::DefaultGeometryTraits<double,coorddim,
              coorddim>
          > ElementGeometry;

      {       // Test a vertex
        std::cout << "==== Testing vertices..." << std::endl;

        affine = true;

        corners.resize(1);

        volume = 1;

        gt.makeVertex();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);
      }       // vertex

    }     // mydim = 0

  }   // coorddim = 0

  ////////////////////////////////////////////////////////////////////////
  //
  //  coorddim = 1
  //

  {
    static const std::size_t coorddim = 1;
    std::cout << "== coorddimension = " << coorddim << std::endl;

    std::vector<FieldVector<double, coorddim> > corners;

    //
    //  mydim = 0
    //

    {
      static const std::size_t mydim = 0;
      std::cout << "=== mydimension = " << mydim << std::endl;

      typedef GenericGeometry::BasicGeometry<
          mydim, GenericGeometry::DefaultGeometryTraits<double,coorddim,
              coorddim>
          > ElementGeometry;

      {       // Test a vertex
        std::cout << "==== Testing vertices..." << std::endl;

        affine = true;

        corners.resize(1);
        corners[0][0] = 0.314;

        volume = 1;

        gt.makeVertex();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);
      }       // vertex

    }     // mydim = 0

    //
    //  mydim = 1
    //

    {
      static const std::size_t mydim = 1;
      std::cout << "=== mydimension = " << mydim << std::endl;

      typedef GenericGeometry::BasicGeometry<
          mydim, GenericGeometry::DefaultGeometryTraits<double,coorddim,
              coorddim>
          > ElementGeometry;

      {       // Test a line segment
        std::cout << "==== Testing line segments..."
                  << std::endl;

        affine = true;

        corners.resize(2);
        corners[0][0] = 0.33;
        corners[1][0] = 0.66;

        volume = 0.33;

        gt.makeLine();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);
      }       // line segment

    }     // mydim = 1

  }   // coorddim = 1

  ////////////////////////////////////////////////////////////////////////
  //
  //  coorddim = 2
  //

  {
    static const std::size_t coorddim = 2;
    std::cout << "== coorddimension = " << coorddim << std::endl;

    std::vector<FieldVector<double, coorddim> > corners;

    //
    //  mydim = 0
    //

    {
      static const std::size_t mydim = 0;
      std::cout << "=== mydimension = " << mydim << std::endl;

      typedef GenericGeometry::BasicGeometry<
          mydim, GenericGeometry::DefaultGeometryTraits<double,coorddim,
              coorddim>
          > ElementGeometry;

      {       // Test a vertex
        std::cout << "==== Testing vertices..." << std::endl;

        affine = true;

        corners.resize(1);
        corners[0][0] = 0.314;  corners[0][1] = 0.27;

        volume = 1;

        gt.makeVertex();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);
      }       // vertex

    }     // mydim = 0

    //
    //  mydim = 1
    //

    {
      static const std::size_t mydim = 1;
      std::cout << "=== mydimension = " << mydim << std::endl;

      typedef GenericGeometry::BasicGeometry<
          mydim, GenericGeometry::DefaultGeometryTraits<double,coorddim,
              coorddim>
          > ElementGeometry;

      {       // Test a line segment
        std::cout << "==== Testing line segments..."
                  << std::endl;

        affine = true;

        corners.resize(2);
        corners[0][0] = 0.33;  corners[0][1] = 0.33;
        corners[1][0] = 0.66;  corners[1][1] = 0.66;

        volume = std::sqrt(2.0)*0.33;

        gt.makeLine();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);
      }       // line segment

    }     // mydim = 1

    //
    //  mydim = 2
    //

    {
      static const std::size_t mydim = 2;
      std::cout << "=== mydimension = " << mydim << std::endl;

      typedef GenericGeometry::BasicGeometry<
          mydim, GenericGeometry::DefaultGeometryTraits<double,coorddim,
              coorddim>
          > ElementGeometry;

      {       // Test a triangle
        std::cout << "==== Testing triangles..."
                  << std::endl;

        affine = true;

        corners.resize(3);
        corners[0][0] = 0.5;   corners[0][1] = 0.0;
        corners[1][0] = 1.0;   corners[1][1] = 0.5;
        corners[2][0] = 0.0;   corners[2][1] = 1.0;

        volume = 3.0/8;

        gt.makeTriangle();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // triangle

      {       // Test an affine quadrilateral
        std::cout << "==== Testing quadrilaterals (affine)..."
                  << std::endl;

        affine = true;

        corners.resize(4);
        corners[0][0] = 0.5;   corners[0][1] = 0.0;
        corners[1][0] = 1.0;   corners[1][1] = 0.5;
        corners[2][0] = 0.0;   corners[2][1] = 0.5;
        corners[3][0] = 0.5;   corners[3][1] = 1.0;

        volume = 0.5;

        gt.makeQuadrilateral();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // affine quadrilateral

      {       // Test a non-affine quadrilateral
        std::cout << "==== Testing quadrilaterals (non-affine)..."
                  << std::endl;

        affine = false;

        corners.resize(4);
        corners[0][0] = 0.5;   corners[0][1] = 0.0;
        corners[1][0] = 1.0;   corners[1][1] = 0.0;
        corners[2][0] = 0.5;   corners[2][1] = 0.25;
        corners[3][0] = 0.75;  corners[3][1] = 0.25;

        volume = 1.5/16;

        gt.makeQuadrilateral();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // non-affine quadrilateral

    }     // mydim = 2

  }   // coorddim = 2

  ////////////////////////////////////////////////////////////////////////
  //
  //  coorddim = 3
  //

  {
    static const std::size_t coorddim = 3;
    std::cout << "== coorddimension = " << coorddim << std::endl;

    std::vector<FieldVector<double, coorddim> > corners;

    //
    //  mydim = 0
    //

    {
      static const std::size_t mydim = 0;
      std::cout << "=== mydimension = " << mydim << std::endl;

      typedef GenericGeometry::BasicGeometry<
          mydim, GenericGeometry::DefaultGeometryTraits<double,coorddim,
              coorddim>
          > ElementGeometry;

      {       // Test a vertex
        std::cout << "==== Testing vertices..." << std::endl;

        affine = true;

        corners.resize(1);
        corners[0][0] = 0.314;  corners[0][1] = 0.27;  corners[0][2] = 0.71;

        volume = 1;

        gt.makeVertex();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);
      }       // vertex

    }     // mydim = 0

    //
    //  mydim = 1
    //

    {
      static const std::size_t mydim = 1;
      std::cout << "=== mydimension = " << mydim << std::endl;

      typedef GenericGeometry::BasicGeometry<
          mydim, GenericGeometry::DefaultGeometryTraits<double,coorddim,
              coorddim>
          > ElementGeometry;

      {       // Test a line segment
        std::cout << "==== Testing line segments..."
                  << std::endl;

        affine = true;

        corners.resize(2);
        corners[0][0] = 0.33; corners[0][1] = 0.33; corners[0][2] = 0.33;
        corners[1][0] = 0.66; corners[1][1] = 0.66; corners[1][2] = 0.66;

        volume = std::sqrt(3.0)*0.33;

        gt.makeLine();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);
      }       // line segment

    }     // mydim = 1

    //
    //  mydim = 2
    //

    {
      static const std::size_t mydim = 2;
      std::cout << "=== mydimension = " << mydim << std::endl;

      typedef GenericGeometry::BasicGeometry<
          mydim, GenericGeometry::DefaultGeometryTraits<double,coorddim,
              coorddim>
          > ElementGeometry;

      {       // Test a triangle
        std::cout << "==== Testing triangles..."
                  << std::endl;

        affine = true;

        corners.resize(3);
        corners[0][0] = 0.5; corners[0][1] = 0.0; corners[0][2] = 1.0;
        corners[1][0] = 1.0; corners[1][1] = 0.5; corners[1][2] = 1.0;
        corners[2][0] = 0.0; corners[2][1] = 1.0; corners[2][2] = 0.0;

        volume = std::sqrt(0.125*(2.25-0.125));

        gt.makeTriangle();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // triangle

      {       // Test an affine quadrilateral
        std::cout << "==== Testing quadrilaterals (affine)..."
                  << std::endl;

        affine = true;

        corners.resize(4);
        corners[0][0] = 0.5; corners[0][1] = 0.0; corners[0][2] = 0.0;
        corners[1][0] = 1.0; corners[1][1] = 0.5; corners[1][2] = 1.0;
        corners[2][0] = 0.0; corners[2][1] = 0.5; corners[2][2] = 0.0;
        corners[3][0] = 0.5; corners[3][1] = 1.0; corners[3][2] = 1.0;

        volume = std::sqrt(0.75);

        gt.makeQuadrilateral();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // affine quadrilateral

      {       // Test a non-affine quadrilateral
        std::cout << "==== Testing quadrilaterals (non-affine)..."
                  << std::endl;

        affine = false;

        corners.resize(4);
        corners[0][0] = 0.0; corners[0][1] = 0.0; corners[0][2] = 0.0;
        corners[1][0] = 1.0; corners[1][1] = 0.0; corners[1][2] = 0.0;
        corners[2][0] = 0.0; corners[2][1] = 1.0; corners[2][2] = 0.0;
        corners[3][0] = 1.0; corners[3][1] = 1.0; corners[3][2] = 1.0;

        // Assuming the surface equation is z=x*y (this satisfies the
        // points above and is bilinear in x and y) I have to
        // integrate
        //
        //   \int_0^1 dx \int_0^1 dy \sqrt(x^2 + y^2 + 1)
        //
        // which evaluates to the value below according to maxima.
        // See
        // http://de.wikipedia.org/wiki/Oberfl%C3%A4chenintegral#Beispiel_2:_Explizite_Darstellung_2
        // and
        // http://de.wikibooks.org/wiki/Diffgeo:_Fl%C3%A4chentheorie:_Fl%C3%A4cheninhalt
        volume = 1.280789271462219;

        gt.makeQuadrilateral();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // non-affine quadrilateral

    }     // mydim = 2

    //
    //  mydim = 3
    //

    {
      static const std::size_t mydim = 3;
      std::cout << "=== mydimension = " << mydim << std::endl;

      typedef GenericGeometry::BasicGeometry<
          mydim, GenericGeometry::DefaultGeometryTraits<double,coorddim,
              coorddim>
          > ElementGeometry;

      {       // Test a tetrahedron
        std::cout << "==== Testing tetrahedron..."
                  << std::endl;

        affine = true;

        corners.resize(4);
        corners[0][0] = 0.5; corners[0][1] = 0.0; corners[0][2] = 1.0;
        corners[1][0] = 1.0; corners[1][1] = 0.5; corners[1][2] = 1.0;
        corners[2][0] = 0.0; corners[2][1] = 1.0; corners[2][2] = 0.0;
        corners[3][0] = 1.0; corners[3][1] = 0.0; corners[3][2] = 0.0;

        // to be determined
        volume = std::numeric_limits<double>::quiet_NaN();

        gt.makeTetrahedron();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // tetrahedron

      {       // Test an affine pyramid
        std::cout << "==== Testing pyramid (affine)..."
                  << std::endl;

        affine = true;

        corners.resize(5);
        corners[0][0] = 0.5; corners[0][1] = 0.0; corners[0][2] = 0.0;
        corners[1][0] = 1.0; corners[1][1] = 0.5; corners[1][2] = 0.0;
        corners[2][0] = 0.0; corners[2][1] = 0.5; corners[2][2] = 0.0;
        corners[3][0] = 0.5; corners[3][1] = 1.0; corners[3][2] = 0.0;
        corners[4][0] = 0.5; corners[4][1] = 0.5; corners[4][2] = 0.5*std::sqrt(0.5);

        volume = std::pow(0.5, 1.5)/6;

        gt.makePyramid();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // affine pyramid

      {       // Test a non-affine pyramid
        std::cout << "==== Testing pyramid (non-affine)..."
                  << std::endl;

        affine = false;

        corners.resize(5);
        corners[0][0] = 0.0; corners[0][1] = 0.0; corners[0][2] = 0.0;
        corners[1][0] = 1.0; corners[1][1] = 0.0; corners[1][2] = 0.0;
        corners[2][0] = 0.0; corners[2][1] = 1.0; corners[2][2] = 0.0;
        corners[3][0] = 1.0; corners[3][1] = 1.0; corners[3][2] = 1.0;
        corners[4][0] = 0.0; corners[4][1] = 0.0; corners[4][2] = 0.5;

        // to be determined
        volume = std::numeric_limits<double>::quiet_NaN();

        gt.makePyramid();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // non-affine pyramid

      {       // Test an affine prism
        std::cout << "==== Testing prism (affine)..."
                  << std::endl;

        affine = true;

        corners.resize(6);
        corners[0][0] = 0.5; corners[0][1] = 0.0; corners[0][2] = 1.0;
        corners[1][0] = 1.0; corners[1][1] = 0.5; corners[1][2] = 1.0;
        corners[2][0] = 0.0; corners[2][1] = 1.0; corners[2][2] = 0.0;
        corners[3][0] = 1.5; corners[3][1] = 1.0; corners[3][2] = 2.0;
        corners[4][0] = 2.0; corners[4][1] = 1.5; corners[4][2] = 2.0;
        corners[5][0] = 1.0; corners[5][1] = 2.0; corners[5][2] = 1.0;

        // spatprodukt
        volume = 0.375;

        gt.makePrism();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // affine prism

      {       // Test a non-affine prism
        std::cout << "==== Testing prism (non-affine)..."
                  << std::endl;

        affine = false;

        corners.resize(6);
        corners[0][0] = 0.5; corners[0][1] = 0.0; corners[0][2] = 0.0;
        corners[1][0] = 1.0; corners[1][1] = 0.5; corners[1][2] = 0.0;
        corners[2][0] = 0.0; corners[2][1] = 1.0; corners[2][2] = 0.0;
        corners[3][0] = 1.0; corners[3][1] = 1.0; corners[3][2] = 1.0;
        corners[4][0] = 0.0; corners[4][1] = 0.5; corners[4][2] = 1.0;
        corners[5][0] = 0.5; corners[5][1] = 0.0; corners[5][2] = 1.0;

        // to be determined
        volume = std::numeric_limits<double>::quiet_NaN();

        gt.makePrism();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // non-affine prism

      {       // Test an affine hexahedron
        std::cout << "==== Testing hexahedron (affine)..."
                  << std::endl;

        affine = true;

        corners.resize(8);
        corners[0][0] = 0.0; corners[0][1] = 0.0; corners[0][2] = 0.0;
        corners[1][0] = 1.0; corners[1][1] = 0.1; corners[1][2] = 0.1;
        corners[2][0] = 0.1; corners[2][1] = 1.0; corners[2][2] = 0.1;
        corners[3][0] = 1.1; corners[3][1] = 1.1; corners[3][2] = 0.2;
        corners[4][0] = 0.1; corners[4][1] = 0.1; corners[4][2] = 1.0;
        corners[5][0] = 1.1; corners[5][1] = 0.2; corners[5][2] = 1.1;
        corners[6][0] = 0.2; corners[6][1] = 1.1; corners[6][2] = 1.1;
        corners[7][0] = 1.2; corners[7][1] = 1.2; corners[7][2] = 1.2;

        // spatprodukt der basisvektoren
        // [(1, .1, .1) x (.1, 1, .1)] * (.1, .1, 1)
        volume = 1.0 - 3*0.1*0.1 + 2*0.1*0.1*0.1;

        gt.makeHexahedron();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // affine hexahedron

      {       // Test a non-affine hexahedron
        std::cout << "==== Testing hexahedron (non-affine)..."
                  << std::endl;

        affine = false;

        corners.resize(8);
        corners[0][0] = 0.0; corners[0][1] = 0.0; corners[0][2] = 0.0;
        corners[1][0] = 1.0; corners[1][1] = 0.0; corners[1][2] = 0.0;
        corners[2][0] = 0.0; corners[2][1] = 1.0; corners[2][2] = 0.0;
        corners[3][0] = 1.0; corners[3][1] = 1.0; corners[3][2] = 0.0;
        // turn the top by 45 degree (and make it smaller)
        corners[4][0] = 0.5; corners[4][1] = 0.0; corners[4][2] = 1.0;
        corners[5][0] = 1.0; corners[5][1] = 0.5; corners[5][2] = 1.0;
        corners[6][0] = 0.0; corners[6][1] = 0.5; corners[6][2] = 1.0;
        corners[7][0] = 0.5; corners[7][1] = 1.0; corners[7][2] = 1.0;

        // to be determined
        volume = std::numeric_limits<double>::quiet_NaN();

        gt.makeHexahedron();
        ElementGeometry insideGeometry( gt, corners );

        testBasicGeometry(insideGeometry, volume, affine, result);

      }       // non-affine hexahedron

    }     // mydim = 2

  }   // coorddim = 3

  return result;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  throw;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  throw;
}
