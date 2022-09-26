// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef  BASICUNITCUBE_HH
#define  BASICUNITCUBE_HH

// StructuredGridFactory cannot replace BasicUnitCube when creating a projected unit cube in test-alberta.cc

#include <dune/grid/common/gridfactory.hh>

// declaration of a basic unit cube that uses the GridFactory
template< int dim >
struct BasicUnitCube;

template<>
struct BasicUnitCube< 1 >
{
  template< class Grid >
  static void insertVertices ( Dune::GridFactory< Grid > &factory, const double left = 0.0, const double right = 1.0 )
  {
    Dune::FieldVector<double,1> pos;

    pos[0] = left;
    factory.insertVertex(pos);

    pos[ 0 ] = right;
    factory.insertVertex(pos);
  }

  template< class Grid >
  static void insertSimplices ( Dune::GridFactory< Grid > &factory )
  {
    std::vector< unsigned int > cornerIDs( 2 );

    cornerIDs[0] = 0;  cornerIDs[1] = 1;
    factory.insertElement( Dune::GeometryTypes::simplex( 1 ), cornerIDs );
  }

  template< class Grid >
  static void insertCubes ( Dune::GridFactory< Grid > &factory )
  {
    std::vector< unsigned int > cornerIDs( 2 );

    cornerIDs[0] = 0;  cornerIDs[1] = 1;
    factory.insertElement( Dune::GeometryTypes::cube( 1 ), cornerIDs );
  }
};

// unit cube in two dimensions with 2 variants: triangle and rectangle elements
template<>
struct BasicUnitCube< 2 >
{
  template< class Grid >
  static void insertVertices ( Dune::GridFactory< Grid > &factory, const double left = 0.0, const double right = 1.0 )
  {
    Dune::FieldVector<double,2> pos;

    pos[0] = right;  pos[1] = left;
    factory.insertVertex(pos);                         /*@\label{uc:iv}@*/

    pos[0] = left;  pos[1] = left;
    factory.insertVertex(pos);

    pos[0] = right;  pos[1] = right;
    factory.insertVertex(pos);

    pos[0] = left;  pos[1] = right;
    factory.insertVertex(pos);
  }

  template< class Grid >
  static void insertSimplices ( Dune::GridFactory< Grid > &factory )
  {
    std::vector< unsigned int > cornerIDs( 3 );

    cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;
    factory.insertElement( Dune::GeometryTypes::simplex( 2 ), cornerIDs );          /*@\label{uc:ie}@*/

    cornerIDs[0] = 2;  cornerIDs[1] = 1;  cornerIDs[2] = 3;
    factory.insertElement( Dune::GeometryTypes::simplex( 2 ), cornerIDs );
  }

  template< class Grid >
  static void insertCubes ( Dune::GridFactory< Grid > &factory )
  {
    std::vector< unsigned int > cornerIDs( 4 );
    for( int i = 0; i < 4; ++i )
      cornerIDs[ i ] = i;
    factory.insertElement( Dune::GeometryTypes::cube( 2 ), cornerIDs );
  }
};

// unit cube in 3 dimensions with two variants: tetraheda and hexahedra
template<>
struct BasicUnitCube< 3 >
{
  template< class Grid >
  static void insertVertices ( Dune::GridFactory< Grid > &factory, const double left = 0.0, const double right = 1.0 )
  {
    Dune::FieldVector< double, 3 > pos;

    pos[0] = right; pos[1] = left;  pos[2] = left;    factory.insertVertex(pos);
    pos[0] = left;  pos[1] = left;  pos[2] = left;    factory.insertVertex(pos);
    pos[0] = right; pos[1] = right; pos[2] = left;    factory.insertVertex(pos);
    pos[0] = left;  pos[1] = right; pos[2] = left;    factory.insertVertex(pos);
    pos[0] = right; pos[1] = left;  pos[2] = right;   factory.insertVertex(pos);
    pos[0] = left;  pos[1] = left;  pos[2] = right;   factory.insertVertex(pos);
    pos[0] = right; pos[1] = right; pos[2] = right;   factory.insertVertex(pos);
    pos[0] = left;  pos[1] = right; pos[2] = right;   factory.insertVertex(pos);
  }

  template< class Grid >
  static void insertSimplices ( Dune::GridFactory< Grid > &factory )
  {
    std::vector< unsigned int > cornerIDs( 4 );

    cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;  cornerIDs[3] = 4;
    factory.insertElement( Dune::GeometryTypes::simplex( 3 ), cornerIDs );

    cornerIDs[0] = 1;  cornerIDs[1] = 3;  cornerIDs[2] = 2;  cornerIDs[3] = 7;
    factory.insertElement( Dune::GeometryTypes::simplex( 3 ), cornerIDs );

    cornerIDs[0] = 1;  cornerIDs[1] = 7;  cornerIDs[2] = 2;  cornerIDs[3] = 4;
    factory.insertElement( Dune::GeometryTypes::simplex( 3 ), cornerIDs );

    cornerIDs[0] = 1;  cornerIDs[1] = 7;  cornerIDs[2] = 4;  cornerIDs[3] = 5;
    factory.insertElement( Dune::GeometryTypes::simplex( 3 ), cornerIDs );

    cornerIDs[0] = 4;  cornerIDs[1] = 7;  cornerIDs[2] = 2;  cornerIDs[3] = 6;
    factory.insertElement( Dune::GeometryTypes::simplex( 3 ), cornerIDs );
  }

  template< class Grid >
  static void insertCubes ( Dune::GridFactory< Grid > &factory )
  {
    std::vector< unsigned int > cornerIDs( 8 );
    for( int i = 0; i < 8; ++i )
      cornerIDs[ i ] = i;
    factory.insertElement( Dune::GeometryTypes::cube( 3 ), cornerIDs );
  }
};

#endif  /*BASICUNITCUBE_HH*/
