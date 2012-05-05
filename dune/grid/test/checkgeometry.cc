// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_CHECK_GEOMETRY_CC
#define DUNE_CHECK_GEOMETRY_CC

#include <limits>

#include <dune/common/forloop.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules/gaussquadrature.hh>

#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridview.hh>

namespace Dune
{

  template< int mydim, int cdim, class Grid, template< int, int, class > class Imp >
  void checkGeometry ( const Geometry< mydim, cdim, Grid, Imp > &geometry )
  {
    typedef typename Grid :: ctype ctype ;
    typedef Dune::Geometry< mydim, cdim, Grid, Imp > Geometry;
    const GenericReferenceElement< ctype, mydim > &refElement = GenericReferenceElements< ctype, mydim >::general( geometry.type() );
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

    typedef Dune::GenericGeometry::GaussPoints< ctype > Points;
    typedef Dune::GenericGeometry::GenericQuadratureFactory< mydim, ctype, Points > QuadratureFactory;
    const typename QuadratureFactory::Object &quadrature = *QuadratureFactory::create( geometry.type(), 2 );
    for( size_t i = 0; i < quadrature.size(); ++i )
    {
      const typename Geometry::LocalCoordinate &x = quadrature[ i ].position();

      if( (x - geometry.local( geometry.global( x ) )).two_norm() > 1e-8 )
        std::cerr << "Error: global and local are not inverse to each other." << std::endl;

      const FieldMatrix< ctype, mydim, cdim > &jt = geometry.jacobianTransposed( x );
      const FieldMatrix< ctype, cdim, mydim > &jit = geometry.jacobianInverseTransposed( x );

      FieldMatrix< ctype, mydim, mydim > id;
      FMatrixHelp::multMatrix( jt, jit, id );
      bool isId = true;
      for( int j = 0; j < mydim; ++j )
        for( int k = 0; k < mydim; ++k )
          isId &= (std::abs( id[ j ][ k ] - (j == k ? 1 : 0) ) < 1e-8);
      if( !isId )
      {
        std::cerr << "Error: jacobianTransposed and jacobianInverseTransposed are not inverse to each other." << std::endl;
        std::cout << "       id != [ ";
        for( int j = 0; j < mydim; ++j )
          std::cout << (j > 0 ? " | " : "") << id[ j ];
        std::cout << " ]" << std::endl;
      }

      if( geometry.integrationElement( x ) < 0 )
        std::cerr << "Error: Negative integrationElement found." << std::endl;

      FieldMatrix< ctype, mydim, mydim > jtj( 0 );
      for( int i = 0; i < mydim; ++i )
        for( int j = 0; j < mydim; ++j )
          for( int k = 0; k < cdim; ++k )
            jtj[ i ][ j ] += jt[ i ][ k ] * jt[ j ][ k ];
      if( std::abs( std::sqrt( jtj.determinant() ) - geometry.integrationElement( x ) ) > 1e-8 )
        std::cerr << "Error: integrationElement is not consistent with jacobianTransposed." << std::endl;
      if (geometry.affine())
        if( std::abs( geometry.volume() - refElement.volume()*geometry.integrationElement( x ) ) > 1e-8 )
          std::cerr << "Error: volume is not consistent with jacobianTransposed." << std::endl;
    }
    QuadratureFactory::release( &quadrature );

    {
      // get reference element
      typedef typename Grid::ctype ctype;
      GeometryType type = geometry.type();
      const GenericReferenceElement< ctype , mydim > & refElement =
        GenericReferenceElements< ctype, mydim >::general(type);
      // center is (for now) the centroid of the reference element mapped to
      // this geometry.
      const FieldVector<ctype, cdim> center = geometry.global(refElement.position(0,0));
      if( std::abs( (geometry.center() - center).two_norm() ) > 1e-8 )
        DUNE_THROW(Exception, "center() is not consistent with global(refElem.position(0,0)).");
    }
  }


  /** \param geometry The local geometry to be tested
   * \param type The type of the element that the local geometry is embedded in
   * \param geoName Helper string that will appear in the error message
   */
  template< int mydim, int cdim, class Grid, template< int, int, class > class Imp >
  void checkLocalGeometry ( const Geometry< mydim, cdim, Grid, Imp > &geometry,
                            GeometryType type, const std::string &geoName = "local geometry" )
  {
    checkGeometry( geometry );

    // check that corners are within the reference element of the given type
    assert( type.dim() == cdim );
    const GenericReferenceElement< typename Grid::ctype, cdim > &refElement
      = GenericReferenceElements< typename Grid::ctype, cdim >::general( type );

    const int numCorners = geometry.corners();
    for( int i = 0; i < numCorners; ++i )
    {
      if( !refElement.checkInside( geometry.corner( i ) ) )
      {
        std::cerr << "Corner " << i
                  << " of " << geoName << " not within reference element: "
                  << geometry.corner( i ) << "." << std::endl;
      }
    }
  }


  template <int codim>
  struct CheckSubEntityGeometry
  {
    template <int dim,class GI,template <int,int,class> class EI>
    static void apply(const Entity<0,dim,GI,EI> &entity)
    {
      integral_constant<
          bool, Dune::Capabilities::hasEntity<GI,codim>::v
          > capVar;
      check(capVar,entity);
    }
    template <class Entity>
    static void check(const true_type&, const Entity &entity)
    {
      for (int i=0; i<entity.template count<codim>(); ++i)
      {
        typedef typename Entity::template Codim< codim >::EntityPointer SubEP;
        const SubEP subEP = entity.template subEntity<codim>(i);
        const typename SubEP::Entity &subEn = *subEP;
        const typename SubEP::Entity::Geometry &subGeo = subEn.geometry();

        if( subEn.type() != subGeo.type() )
          std::cerr << "Error: Entity and geometry report different geometry types on codimension " << codim << "." << std::endl;
        checkGeometry(subGeo);
      }
    }
    template <class Entity>
    static void check(const false_type&, const Entity &entity)
    {}
  };

  template<typename GV>
  void checkGeometryLifetime (const GV &gridView)
  {
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename GV::template Codim<0>::Geometry Geometry;
    typedef typename GV::ctype ctype;
    enum { dim  = GV::dimension };
    enum { dimw = GV::dimensionworld };

    const FieldVector<ctype, dim> pos(0.2);
    const FieldVector<ctype, dimw> glob =
      gridView.template begin<0>()->geometry().global(pos);

    Iterator it = gridView.template begin<0>();
    const Geometry geomCopy = it->geometry();

    const Iterator end = gridView.template end<0>();
    checkGeometry ( geomCopy );

    for(it = gridView.template begin<0>();
        it != end; ++it )
    {
      // due to register/memory differences we might have
      // errors < 1e-16
      assert (std::abs((geomCopy.global(pos) - glob).one_norm()) < std::numeric_limits<ctype>::epsilon());
    }
  }

  template< class VT >
  void checkGeometry ( const GridView< VT > &gridView )
  {
    typedef typename GridView< VT >::template Codim<0>::Iterator Iterator;
    typedef typename GridView< VT >::template Codim<0>::Geometry Geometry;

    const Iterator end = gridView.template end<0>();
    Iterator it = gridView.template begin<0>();

    for( ; it != end; ++it )
    {
      ForLoop<CheckSubEntityGeometry,0,GridView<VT>::dimension>::apply(*it);
    }

  }

}

#endif // #ifndef DUNE_CHECK_GEOMETRY_CC
