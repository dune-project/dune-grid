// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CHECK_GEOMETRY_CC
#define DUNE_CHECK_GEOMETRY_CC

#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/quadraturerules.hh>

namespace Dune
{

  template< int mydim, int cdim, class Grid, template< int, int, class > class Imp >
  void checkGeometry ( const Geometry< mydim, cdim, Grid, Imp > &geometry )
  {
    const GenericReferenceElement< double, mydim > &refElement = GenericReferenceElements< double, mydim >::general( geometry.type() );
    if( refElement.size( mydim ) == geometry.corners() )
    {
      for( int i = 0; i < geometry.corners(); ++i )
      {
        if( (geometry.corner( i ) - geometry.global( refElement.position( i, mydim ) )).two_norm() > 1e-8 )
          std::cerr << "Error: Methods corner and global are inconsistent." << std::endl;
      }
    }
    else
      std::cerr << "Error: Incorrect number of corners (" << geometry.corners() << ", should be " << refElement.size( mydim ) << ")." << std::endl;

    const QuadratureRule< double, mydim > &quadrature = QuadratureRules< double, mydim >::rule( geometry.type(), 2 );
    for( size_t i = 0; i < quadrature.size(); ++i )
    {
      const FieldVector< double, mydim > &x = quadrature[ i ].position();

      if( (x - geometry.local( geometry.global( x ) )).two_norm() > 1e-8 )
        std::cerr << "Error: global and local are not inverse to each other." << std::endl;

      const FieldMatrix< double, mydim, cdim > &jt = geometry.jacobianTransposed( x );
      const FieldMatrix< double, cdim, mydim > &jit = geometry.jacobianInverseTransposed( x );

      FieldMatrix< double, mydim, mydim > id;
      FMatrixHelp::multMatrix( jt, jit, id );
      bool isId = true;
      for( int i = 0; i < mydim; ++i )
        for( int j = 0; j < mydim; ++j )
          isId &= (std::abs( id[ i ][ j ] - (i == j ? 1 : 0) ) < 1e-8);
      if( !isId )
        std::cerr << "Error: jacobianTransposed and jacobianInverseTransposed are not inverse to each other." << std::endl;

      if( geometry.integrationElement( x ) < 0 )
        std::cerr << "Error: Negative integrationElement found." << std::endl;

      FieldMatrix< double, mydim, mydim > jtj( 0 );
      for( int i = 0; i < mydim; ++i )
        for( int j = 0; j < mydim; ++j )
          for( int k = 0; k < cdim; ++k )
            jtj[ i ][ j ] += jt[ i ][ k ] * jt[ j ][ k ];
      if( std::abs( std::sqrt( jtj.determinant() ) - geometry.integrationElement( x ) ) > 1e-8 )
        std::cerr << "Error: integrationElement is not consistent with jacobianTransposed." << std::endl;
    }
  }


  template< class VT >
  void checkGeometry ( const GridView< VT > &gridView )
  {
    typedef typename GridView< VT >::template Codim< 0 >::Iterator Iterator;

    const Iterator end = gridView.template end< 0 >();
    for( Iterator it = gridView.template begin< 0 >(); it != end; ++it )
    {
      const typename Iterator::Entity &entity = *it;
      const typename Iterator::Entity::Geometry &geometry = entity.geometry();

      if( entity.type() != geometry.type() )
        std::cerr << "Error: Entity and geometry report different geometry types." << std::endl;

      checkGeometry( geometry );
    }
  }

}

#endif // #ifndef DUNE_CHECK_GEOMETRY_CC
