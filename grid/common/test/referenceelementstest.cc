// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id: mcmgmappertest.cc 2793 2006-07-10 12:01:52Z mblatt $

/** \file
    \brief A unit test for the ReferenceElements
 */

#include <config.h>

#include <iostream>

#include <dune/grid/common/referenceelements.hh>

using namespace Dune;

int main () try
{

  // //////////////////////////////////////////////////////////////////////////
  //   We test the different elements separately.  I will not try to be smart
  //   and write routines that check all types at once.  That would just be
  //   illegible.
  // //////////////////////////////////////////////////////////////////////////

  GeometryType type;

  // //////////////////////////////////////////////////////////////////////////
  //   Test segment
  // //////////////////////////////////////////////////////////////////////////

  type.makeLine();

  const ReferenceElement<double,1>& referenceLine = ReferenceElements<double, 1>::general(type);

  // size(int c)
  assert(referenceLine.size(0)==1);
  assert(referenceLine.size(1)==2);

  // size(int i, int c, int cc)
  assert(referenceLine.size(0,0,0)==1);
  assert(referenceLine.size(0,0,1)==2);
  assert(referenceLine.size(0,1,1)==1);
  assert(referenceLine.size(1,1,1)==1);

  // subEntity(int i, int c, int ii, int cc)

  // type(int i, int c)
  assert(referenceLine.type(0,0).isLine());
  assert(referenceLine.type(0,1).isVertex());
  assert(referenceLine.type(1,1).isVertex());


  // //////////////////////////////////////////////////////////////////////////
  //   Test triangle
  // //////////////////////////////////////////////////////////////////////////

  type.makeTriangle();

  const ReferenceElement<double,2>& referenceTriangle = ReferenceElements<double, 2>::general(type);

  // size(int c)
  assert(referenceTriangle.size(0)==1);
  assert(referenceTriangle.size(1)==3);
  assert(referenceTriangle.size(2)==3);

  // size(int i, int c, int cc)
  assert(referenceTriangle.size(0,0,0)==1);
  assert(referenceTriangle.size(0,0,1)==3);
  assert(referenceTriangle.size(0,0,2)==3);

  assert(referenceTriangle.size(0,1,1)==1);
  assert(referenceTriangle.size(0,1,2)==2);
  assert(referenceTriangle.size(1,1,1)==1);
  assert(referenceTriangle.size(1,1,2)==2);
  assert(referenceTriangle.size(2,1,1)==1);
  assert(referenceTriangle.size(2,1,2)==2);

  assert(referenceTriangle.size(0,2,2)==1);
  assert(referenceTriangle.size(1,2,2)==1);
  assert(referenceTriangle.size(2,2,2)==1);

  // subEntity(int i, int c, int ii, int cc)
  assert(referenceTriangle.subEntity(0,0,0,0)==0);
  assert(referenceTriangle.subEntity(0,0,0,1)==0);
  assert(referenceTriangle.subEntity(0,0,1,1)==1);
  assert(referenceTriangle.subEntity(0,0,2,1)==2);
  assert(referenceTriangle.subEntity(0,0,0,2)==0);
  assert(referenceTriangle.subEntity(0,0,1,2)==1);
  assert(referenceTriangle.subEntity(0,0,2,2)==2);

  assert(referenceTriangle.subEntity(0,1,0,1)==0);
  assert(referenceTriangle.subEntity(1,1,0,1)==1);
  assert(referenceTriangle.subEntity(2,1,0,1)==2);

  assert(referenceTriangle.subEntity(0,1,0,2)==1);
  assert(referenceTriangle.subEntity(0,1,1,2)==2);
  assert(referenceTriangle.subEntity(1,1,0,2)==2);
  assert(referenceTriangle.subEntity(1,1,1,2)==0);
  assert(referenceTriangle.subEntity(2,1,0,2)==0);
  assert(referenceTriangle.subEntity(2,1,1,2)==1);

  assert(referenceTriangle.subEntity(0,2,0,2)==0);
  assert(referenceTriangle.subEntity(1,2,0,2)==1);
  assert(referenceTriangle.subEntity(2,2,0,2)==2);

  // type(int i, int c)
  assert(referenceTriangle.type(0,0).isTriangle());

  assert(referenceTriangle.type(0,1).isLine());
  assert(referenceTriangle.type(1,1).isLine());
  assert(referenceTriangle.type(2,1).isLine());

  assert(referenceTriangle.type(0,2).isVertex());
  assert(referenceTriangle.type(1,2).isVertex());
  assert(referenceTriangle.type(2,2).isVertex());


  // //////////////////////////////////////////////////////////////////////////
  //   Test quadrilateral
  // //////////////////////////////////////////////////////////////////////////

  type.makeQuadrilateral();

  const ReferenceElement<double,2>& referenceQuad = ReferenceElements<double, 2>::general(type);

  // size(int c)
  assert(referenceQuad.size(0)==1);
  assert(referenceQuad.size(1)==4);
  assert(referenceQuad.size(2)==4);

  // size(int i, int c, int cc)
  assert(referenceQuad.size(0,0,0)==1);
  assert(referenceQuad.size(0,0,1)==4);
  assert(referenceQuad.size(0,0,2)==4);

  assert(referenceQuad.size(0,1,1)==1);
  assert(referenceQuad.size(0,1,2)==2);
  assert(referenceQuad.size(1,1,1)==1);
  assert(referenceQuad.size(1,1,2)==2);
  assert(referenceQuad.size(2,1,1)==1);
  assert(referenceQuad.size(2,1,2)==2);
  assert(referenceQuad.size(3,1,1)==1);
  assert(referenceQuad.size(3,1,2)==2);

  assert(referenceQuad.size(0,2,2)==1);
  assert(referenceQuad.size(1,2,2)==1);
  assert(referenceQuad.size(2,2,2)==1);
  assert(referenceQuad.size(3,2,2)==1);

  // subEntity(int i, int c, int ii, int cc)
  assert(referenceQuad.subEntity(0,0,0,0)==0);
  assert(referenceQuad.subEntity(0,0,0,1)==0);
  assert(referenceQuad.subEntity(0,0,1,1)==1);
  assert(referenceQuad.subEntity(0,0,2,1)==2);
  assert(referenceQuad.subEntity(0,0,3,1)==3);
  assert(referenceQuad.subEntity(0,0,0,2)==0);
  assert(referenceQuad.subEntity(0,0,1,2)==1);
  assert(referenceQuad.subEntity(0,0,2,2)==2);
  assert(referenceQuad.subEntity(0,0,3,2)==3);

  assert(referenceQuad.subEntity(0,1,0,1)==0);
  assert(referenceQuad.subEntity(1,1,0,1)==1);
  assert(referenceQuad.subEntity(2,1,0,1)==2);
  assert(referenceQuad.subEntity(3,1,0,1)==3);

  assert(referenceQuad.subEntity(0,1,0,2)==2);
  assert(referenceQuad.subEntity(0,1,1,2)==0);
  assert(referenceQuad.subEntity(1,1,0,2)==1);
  assert(referenceQuad.subEntity(1,1,1,2)==3);
  assert(referenceQuad.subEntity(2,1,0,2)==0);
  assert(referenceQuad.subEntity(2,1,1,2)==1);
  assert(referenceQuad.subEntity(3,1,0,2)==3);
  assert(referenceQuad.subEntity(3,1,1,2)==2);

  assert(referenceQuad.subEntity(0,2,0,2)==0);
  assert(referenceQuad.subEntity(1,2,0,2)==1);
  assert(referenceQuad.subEntity(2,2,0,2)==2);
  assert(referenceQuad.subEntity(3,2,0,2)==3);

  // type(int i, int c)
  assert(referenceQuad.type(0,0).isQuadrilateral());

  assert(referenceQuad.type(0,1).isLine());
  assert(referenceQuad.type(1,1).isLine());
  assert(referenceQuad.type(2,1).isLine());
  assert(referenceQuad.type(3,1).isLine());

  assert(referenceQuad.type(0,2).isVertex());
  assert(referenceQuad.type(1,2).isVertex());
  assert(referenceQuad.type(2,2).isVertex());
  assert(referenceQuad.type(3,2).isVertex());


  // //////////////////////////////////////////////////////////////////////////
  //   Test tetrahedron
  // //////////////////////////////////////////////////////////////////////////

  type.makeTetrahedron();

  const ReferenceElement<double,3>& referenceTetra = ReferenceElements<double, 3>::general(type);

  // size(int c)
  assert(referenceTetra.size(0)==1);
  assert(referenceTetra.size(1)==4);
  assert(referenceTetra.size(2)==6);
  assert(referenceTetra.size(3)==4);

  // size(int i, int c, int cc)
  assert(referenceTetra.size(0,0,0)==1);
  assert(referenceTetra.size(0,0,1)==4);
  assert(referenceTetra.size(0,0,2)==6);
  assert(referenceTetra.size(0,0,3)==4);

  for (int i=0; i<referenceTetra.size(1); i++) {
    assert(referenceTetra.size(i,1,1)==1);
    assert(referenceTetra.size(i,1,2)==3);
    assert(referenceTetra.size(i,1,3)==3);
  }

  for (int i=0; i<referenceTetra.size(2); i++) {
    assert(referenceTetra.size(i,2,2)==1);
    assert(referenceTetra.size(i,2,3)==2);
  }

  for (int i=0; i<referenceTetra.size(3); i++)
    assert(referenceTetra.size(i,3,3)==1);

  // subEntity(int i, int c, int ii, int cc)


  // type(int i, int c)
  assert(referenceTetra.type(0,0).isTetrahedron());

  for (int i=0; i<referenceTetra.size(1); i++)
    assert(referenceTetra.type(i,1).isTriangle());

  for (int i=0; i<referenceTetra.size(2); i++)
    assert(referenceTetra.type(i,2).isLine());

  for (int i=0; i<referenceTetra.size(3); i++)
    assert(referenceTetra.type(i,3).isVertex());

  // //////////////////////////////////////////////////////////////////////////
  //   Test pyramid
  // //////////////////////////////////////////////////////////////////////////

  type.makePyramid();

  const ReferenceElement<double,3>& referencePyramid = ReferenceElements<double, 3>::general(type);

  // size(int c)
  assert(referencePyramid.size(0)==1);
  assert(referencePyramid.size(1)==5);
  assert(referencePyramid.size(2)==8);
  assert(referencePyramid.size(3)==5);

  // size(int i, int c, int cc)
  assert(referencePyramid.size(0,0,0)==1);
  assert(referencePyramid.size(0,0,1)==5);
  assert(referencePyramid.size(0,0,2)==8);
  assert(referencePyramid.size(0,0,3)==5);

  assert(referencePyramid.size(0,1,1)==1);
  assert(referencePyramid.size(0,1,2)==4);
  assert(referencePyramid.size(0,1,3)==4);
  assert(referencePyramid.size(1,1,1)==1);
  assert(referencePyramid.size(1,1,2)==3);
  assert(referencePyramid.size(1,1,3)==3);
  assert(referencePyramid.size(2,1,1)==1);
  assert(referencePyramid.size(2,1,2)==3);
  assert(referencePyramid.size(2,1,3)==3);
  assert(referencePyramid.size(3,1,1)==1);
  assert(referencePyramid.size(3,1,2)==3);
  assert(referencePyramid.size(3,1,3)==3);
  assert(referencePyramid.size(4,1,1)==1);
  assert(referencePyramid.size(4,1,2)==3);
  assert(referencePyramid.size(4,1,3)==3);

  for (int i=0; i<referencePyramid.size(2); i++) {
    assert(referencePyramid.size(i,2,2)==1);
    assert(referencePyramid.size(i,2,3)==2);
  }

  for (int i=0; i<referencePyramid.size(3); i++)
    assert(referencePyramid.size(i,3,3)==1);

  // subEntity(int i, int c, int ii, int cc)

  // type(int i, int c)
  assert(referencePyramid.type(0,0).isPyramid());

  assert(referencePyramid.type(0,1).isQuadrilateral());
  assert(referencePyramid.type(1,1).isTriangle());
  assert(referencePyramid.type(2,1).isTriangle());
  assert(referencePyramid.type(3,1).isTriangle());
  assert(referencePyramid.type(4,1).isTriangle());

  for (int i=0; i<referencePyramid.size(2); i++)
    assert(referencePyramid.type(i,2).isLine());

  for (int i=0; i<referencePyramid.size(3); i++)
    assert(referencePyramid.type(i,3).isVertex());


  // //////////////////////////////////////////////////////////////////////////
  //   Test prism
  // //////////////////////////////////////////////////////////////////////////

  type.makePrism();

  const ReferenceElement<double,3>& referencePrism = ReferenceElements<double, 3>::general(type);

  // size(int c)
  assert(referencePrism.size(0)==1);
  assert(referencePrism.size(1)==5);
  assert(referencePrism.size(2)==9);
  assert(referencePrism.size(3)==6);

  // size(int i, int c, int cc)
  assert(referencePrism.size(0,0,0)==1);
  assert(referencePrism.size(0,0,1)==5);
  assert(referencePrism.size(0,0,2)==9);
  assert(referencePrism.size(0,0,3)==6);

  assert(referencePrism.size(0,1,1)==1);
  assert(referencePrism.size(0,1,2)==3);
  assert(referencePrism.size(0,1,3)==3);
  assert(referencePrism.size(1,1,1)==1);
  assert(referencePrism.size(1,1,2)==4);
  assert(referencePrism.size(1,1,3)==4);
  assert(referencePrism.size(2,1,1)==1);
  assert(referencePrism.size(2,1,2)==4);
  assert(referencePrism.size(2,1,3)==4);
  assert(referencePrism.size(3,1,1)==1);
  assert(referencePrism.size(3,1,2)==4);
  assert(referencePrism.size(3,1,3)==4);
  assert(referencePrism.size(4,1,1)==1);
  assert(referencePrism.size(4,1,2)==3);
  assert(referencePrism.size(4,1,3)==3);

  for (int i=0; i<referencePrism.size(2); i++) {
    assert(referencePrism.size(i,2,2)==1);
    assert(referencePrism.size(i,2,3)==2);
  }

  for (int i=0; i<referencePrism.size(3); i++)
    assert(referencePrism.size(i,3,3)==1);

  // subEntity(int i, int c, int ii, int cc)

  // type(int i, int c)
  assert(referencePrism.type(0,0).isPrism());

  assert(referencePrism.type(0,1).isTriangle());
  assert(referencePrism.type(1,1).isQuadrilateral());
  assert(referencePrism.type(2,1).isQuadrilateral());
  assert(referencePrism.type(3,1).isQuadrilateral());
  assert(referencePrism.type(4,1).isTriangle());

  for (int i=0; i<referencePrism.size(2); i++)
    assert(referencePrism.type(i,2).isLine());

  for (int i=0; i<referencePrism.size(3); i++)
    assert(referencePrism.type(i,3).isVertex());

  // //////////////////////////////////////////////////////////////////////////
  //   Test hexahedron
  // //////////////////////////////////////////////////////////////////////////

  type.makeHexahedron();

  const ReferenceElement<double,3>& referenceHexa = ReferenceElements<double, 3>::general(type);

  // size(int c)
  assert(referenceHexa.size(0)==1);
  assert(referenceHexa.size(1)==6);
  assert(referenceHexa.size(2)==12);
  assert(referenceHexa.size(3)==8);

  // size(int i, int c, int cc)
  assert(referenceHexa.size(0,0,0)==1);
  assert(referenceHexa.size(0,0,1)==6);
  assert(referenceHexa.size(0,0,2)==12);
  assert(referenceHexa.size(0,0,3)==8);

  for (int i=0; i<referenceHexa.size(1); i++) {
    assert(referenceHexa.size(i,1,1)==1);
    assert(referenceHexa.size(i,1,2)==4);
    assert(referenceHexa.size(i,1,3)==4);
  }

  for (int i=0; i<referenceHexa.size(2); i++) {
    assert(referenceHexa.size(i,2,2)==1);
    assert(referenceHexa.size(i,2,3)==2);
  }

  for (int i=0; i<referenceHexa.size(3); i++)
    assert(referenceHexa.size(i,3,3)==1);

  // subEntity(int i, int c, int ii, int cc)

  // type(int i, int c)
  assert(referenceHexa.type(0,0).isHexahedron());

  for (int i=0; i<referenceHexa.size(1); i++)
    assert(referenceHexa.type(i,1).isQuadrilateral());

  for (int i=0; i<referenceHexa.size(2); i++)
    assert(referenceHexa.type(i,2).isLine());

  for (int i=0; i<referenceHexa.size(3); i++)
    assert(referenceHexa.type(i,3).isVertex());


  return 0;

}
catch (Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
