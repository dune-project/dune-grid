// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/common/forloop.hh>

#if 0
namespace Dune
{
  template<>
  ReferenceCubeContainer< double, 4 > ReferenceElements< double, 4 >::cube;
  template<>
  ReferenceSimplexContainer< double, 4 > ReferenceElements< double, 4 >::simplices;
  template<>
  ReferenceElementContainer< double, 4 > ReferenceElements< double, 4 >::general;
}
#endif


template< int dim >
struct GenericReferenceElements
{
  static const unsigned int dimension = dim;

  typedef Dune::GenericReferenceElement< double, dim > ReferenceElement;

  static const ReferenceElement &get ( const Dune::GeometryType &type )
  {
    return Dune::GenericReferenceElements< double, dim >::general( type );
  }
};

template< int dim >
struct ReferenceElements
{
  static const unsigned int dimension = dim;

  typedef Dune::ReferenceElement< double, dim > ReferenceElement;

  static const ReferenceElement &get ( const Dune::GeometryType &type )
  {
    return Dune::ReferenceElements< double, dim >::general( type );
  }
};


template< template< int > class RefElements, int dim >
struct CheckReferenceNumbering
{
  typedef typename RefElements< dim >::ReferenceElement ReferenceElement;

  template< int codim >
  struct Codim
  {
    typedef typename RefElements< dim-codim >::ReferenceElement SubReferenceElement;

    static void apply ( const ReferenceElement &refElement )
    {
      const int count = refElement.size( codim );
      for( int i = 0; i < count; ++i )
      {
        const Dune::GeometryType type = refElement.type( i, codim );
        const SubReferenceElement &subRefElement = RefElements< dim-codim >::get( type );

        for( int cc = codim; cc <= dim; ++cc )
        {
          const int subCount = refElement.size( i, codim, cc );
          if( subCount != subRefElement.size( cc-codim ) )
          {
            std::cerr << "Codim " << codim << ", subentity " << i
                      << ": Size of subcodim " << cc << " is inconsistent."
                      << std::endl;
            continue;
          }

          for( int j = 0; j < subCount; ++j )
          {
            const int k = refElement.subEntity( i, codim, j, cc );

            const int vertexCount = refElement.size( k, cc, dim );
            if( vertexCount != subRefElement.size( j, cc-codim, dim-codim ) )
            {
              std::cerr << "Codim " << codim << ", subentity " << i
                        << ": Inconsistent number of vertices for "
                        << "(" << j << ", " << cc << ")." << std::endl;
              continue;
            }

            for( int v = 0; v < vertexCount; ++v )
            {
              const int v1 = refElement.subEntity( k, cc, v, dim );
              const int sv = subRefElement.subEntity( j, cc-codim, v, dim-codim );
              const int v2 = refElement.subEntity( i, codim, sv, dim );
              if( v1 != v2 )
              {
                std::cerr << "Codim " << codim << ", subentity " << i
                          << ": Vertex " << v << " inconsistent vertices for "
                          << "(" << j << ", " << cc << "): "
                          << v1 << " != " << v2 << std::endl;
              }
            }
          }
        }
      }
    }
  };

  static void apply ( const Dune::GeometryType &type )
  {
    const ReferenceElement &refElement = RefElements< dim >::get( type );
    Dune::ForLoop< Codim, 0, dim >::apply( refElement );
    std::cout << std::endl;
  }
};


int main ( int argv, char **argc )
{
  {
    Dune::GeometryType type( Dune::GeometryType::simplex, 1 );
    std::cout << ">>> Checking Dune Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< ReferenceElements, 1 >::apply( type );
    std::cout << ">>> Checking Generic Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< GenericReferenceElements, 1 >::apply( type );
  }

  {
    Dune::GeometryType type( Dune::GeometryType::cube, 1 );
    std::cout << ">>> Checking Dune Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< ReferenceElements, 1 >::apply( type );
    std::cout << ">>> Checking Generic Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< GenericReferenceElements, 1 >::apply( type );
  }

  {
    Dune::GeometryType type( Dune::GeometryType::simplex, 2 );
    std::cout << ">>> Checking Dune Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< ReferenceElements, 2 >::apply( type );
    std::cout << ">>> Checking Generic Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< GenericReferenceElements, 2 >::apply( type );
  }

  {
    Dune::GeometryType type( Dune::GeometryType::cube, 2 );
    std::cout << ">>> Checking Dune Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< ReferenceElements, 2 >::apply( type );
    std::cout << ">>> Checking Generic Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< GenericReferenceElements, 2 >::apply( type );
  }

  {
    Dune::GeometryType type( Dune::GeometryType::simplex, 3 );
    std::cout << ">>> Checking Dune Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< ReferenceElements, 3 >::apply( type );
    std::cout << ">>> Checking Generic Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< GenericReferenceElements, 3 >::apply( type );
  }

  {
    Dune::GeometryType type( Dune::GeometryType::cube, 3 );
    std::cout << ">>> Checking Dune Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< ReferenceElements, 3 >::apply( type );
    std::cout << ">>> Checking Generic Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< GenericReferenceElements, 3 >::apply( type );
  }

  {
    Dune::GeometryType type( Dune::GeometryType::prism, 3 );
    std::cout << ">>> Checking Dune Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< ReferenceElements, 3 >::apply( type );
    std::cout << ">>> Checking Generic Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< GenericReferenceElements, 3 >::apply( type );
  }

  {
    Dune::GeometryType type( Dune::GeometryType::pyramid, 3 );
    std::cout << ">>> Checking Dune Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< ReferenceElements, 3 >::apply( type );
    std::cout << ">>> Checking Generic Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< GenericReferenceElements, 3 >::apply( type );
  }

  {
    Dune::GeometryType type( Dune::GeometryType::simplex, 4 );
    std::cout << ">>> Checking Generic Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< GenericReferenceElements, 4 >::apply( type );
  }

  {
    Dune::GeometryType type( Dune::GeometryType::cube, 4 );
    //std::cout << ">>> Checking Dune Reference Elements for " << type << std::endl;
    //CheckReferenceNumbering< ReferenceElements, 4 >::apply( type );
    std::cout << ">>> Checking Generic Reference Elements for " << type << std::endl;
    CheckReferenceNumbering< GenericReferenceElements, 4 >::apply( type );
  }
}
