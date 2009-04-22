// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <set>

#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/genericreferenceelements.hh>

using namespace Dune;

// Given a gt this returns the generic index of the subsubentity (j,c)
// of subentity (i,codim) of gt relative to the subentity (i,codim)
template<class T, int dim, int CODIM=dim>
struct SubEntityNumbering
{
  static unsigned int dune2generic(Dune::GeometryType gt, int i, int codim, int j, int c)
  {
    if (codim==CODIM)
    {
      const ReferenceElement<T, dim> &ref = ReferenceElements<T,dim>::general(gt);

      typedef Dune::GenericGeometry::MapNumberingProvider<dim-CODIM> Numbering;
      const unsigned int subEntityGtId = Dune::GenericGeometry::topologyId(ref.type(i,codim));

      return Numbering::dune2generic(subEntityGtId, j, (dim-codim)-(dim-c));
    }
    return SubEntityNumbering<T,dim,CODIM-1>::dune2generic(gt, i, codim, j, c);
  }
};

template<class T, int dim>
struct SubEntityNumbering<T,dim,-1>
{
  static unsigned int dune2generic(Dune::GeometryType gt, int i, int codim, int j, int c)
  {
    DUNE_THROW(RangeError, "Wrong codim");
  }
};


template<class T, int dim>
void testWithGeometryType(Dune::GeometryType gt, bool checkSubentityOrientation)
{
  typedef Dune::GenericGeometry::MapNumberingProvider<dim> Numbering;
  const unsigned int gtId = Dune::GenericGeometry::topologyId(gt);


  const ReferenceElement<T, dim> &ref = ReferenceElements<T,dim>::general(gt);
  const GenericReferenceElement<T, dim> &gref = GenericReferenceElements<T,dim>::general(gt);

  std::cout << ">>> Testing for " << gt << std::endl;

  // test volume
  if (ref.volume() != gref.volume())
    std::cout << "ref.volume() != gref.volume()" << std::endl;

  for(int codim=0; codim<=dim; ++codim)
  {
    std::cout << "codim " << codim << std::endl;

    //test size
    if (ref.size(codim) != gref.size(codim))
      std::cout << "  ref.size(" << codim << ") != gref.size(" << codim << ")" << std::endl;

    for(int i=0; i<ref.size(codim); ++i)
    {
      int gi = Numbering::dune2generic(gtId, i, codim);

      std::cout << "  testing methods for subentity(" << i << "," << codim << ") / (" << gi << "," << codim << ")" << std::endl;

      //test position
      if (ref.position(i, codim) != gref.position(gi, codim))
        std::cout << "    ref.position(" << i << "," << codim <<") != gref.position(" << gi << "," << codim << ")" << std::endl;

      //test codim
      if (ref.type(i, codim) != gref.type(gi, codim))
      {
        std::cout << "    ref.type(" << i << "," << codim << ") != gref.type(" << gi << "," << codim << ")" << std::endl;
        std::cout << "    ref.type(" << i << "," << codim << ") = " << ref.type(i, codim) << std::endl;
        std::cout << "    gref.type(" << gi << "," << codim << ") = " << gref.type(gi, codim) << std::endl;
      }

      for(int c=codim+1; c<=dim; ++c)
      {
        std::multiset<unsigned int> subSubEntities;
        std::multiset<unsigned int> genericSubSubEntities;

        //test size
        if (ref.size(i, codim, c) != gref.size(gi, codim, c))
          std::cout << "      ref.size(" <<i << "," << codim<< "," << c << ") != gref.size(" << gi << "," << codim<< "," << c << ")" << std::endl;

        for(int j=0; j<ref.size(i, codim, c); ++j)
        {
          int gj = SubEntityNumbering<T,dim,dim>::dune2generic(gt, i, codim, j, c);

          unsigned int subEntityIndex = ref.subEntity(i, codim, j, c);
          unsigned int genericSubEntityIndex = gref.subEntity(gi, codim, gj, c);

          subSubEntities.insert(ref.subEntity(i, codim, j, c));
          genericSubSubEntities.insert(Numbering::generic2dune(gtId, gref.subEntity(gi, codim, gj, c), c));

          if (checkSubentityOrientation)
          {
            // test subEntity uses the same orientation
            if (ref.subEntity(i, codim, j, c) != Numbering::generic2dune(gtId, gref.subEntity(gi, codim, gj, c), c))
            {
              std::cout << "        ref.subEntity(" <<i << "," << codim<< "," << j << "," << c << ") != Numbering::generic2dune(gtId, gref.subEntity(" << gi << "," << codim << "," << gj << "," << c << ")," << c << ")" << std::endl;
              std::cout << "        ref.subEntity(" << i << "," << codim << "," << j << "," << c << ") = " << ref.subEntity(i, codim, j, c) << std::endl;
              std::cout << "        gref.subEntity(" << gi << "," << codim << "," << gj << "," << c << ") = " << gref.subEntity(gi, codim, gj, c) << std::endl;
              std::cout << "        Numbering::generic2dune(gtId, gref.subEntity(" << gi << "," << codim << "," << gj << "," << c << ")," << c << ") = " << Numbering::generic2dune(gtId, gref.subEntity(gi, codim, gj, c), c) << std::endl;
            }
          }
        }

        //test if subEntity returns the same subsubentities
        if (subSubEntities!=genericSubSubEntities)
          std::cout << "      subsubentities of subentity differ" << std::endl;
      }
    }
  }

  std::cout << std::endl;
}


int main ( int argc, char **argv )
{
#if 0
  const GeometryType type( GeometryType::simplex, 3 );
  const GenericReferenceElement< double, 3 > &refElement
    = GenericReferenceElements< double, 3 >::general( type );
#endif

  typedef Dune::GeometryType GT;

  //    bool checkSubentityOrientation = (argc>1);
  bool checkSubentityOrientation = true;

  testWithGeometryType<double,1>(GT(GT::simplex, 1), checkSubentityOrientation);
  testWithGeometryType<double,2>(GT(GT::simplex, 2), checkSubentityOrientation);
  testWithGeometryType<double,3>(GT(GT::simplex, 3), checkSubentityOrientation);

  testWithGeometryType<double,1>(GT(GT::cube, 1), checkSubentityOrientation);
  testWithGeometryType<double,2>(GT(GT::cube, 2), checkSubentityOrientation);
  testWithGeometryType<double,3>(GT(GT::cube, 3), checkSubentityOrientation);

  testWithGeometryType<double,3>(GT(GT::prism, 3), checkSubentityOrientation);

  testWithGeometryType<double,3>(GT(GT::pyramid, 3), checkSubentityOrientation);

  return 0;
}
