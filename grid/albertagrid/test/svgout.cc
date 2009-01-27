// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <fstream>
#include <iostream>
#include <sstream>

#include <dune/grid/albertagrid.hh>


class SVGWriter
{
  typedef Dune::FieldVector< double, 2 > Coordinate;

  std::ofstream csvg_;
  double shrink_;
  double scale_;
  Coordinate origin_;

public:
  SVGWriter ( const std::string &filename )
    : csvg_( filename.c_str() ),
      shrink_( 1.0 ),
      scale_( 300 ),
      origin_( 0.0 )
  {
    csvg_ << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
    csvg_ << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << std::endl;
    csvg_ << std::endl;
    csvg_ << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" version=\"1.1\" baseProfile=\"full\" >" << std::endl;
  }

  ~SVGWriter ()
  {
    csvg_ << std::endl;
    csvg_ << "</svg>" << std::endl;
    csvg_.close();
  }

  void shrink( double shrink )
  {
    shrink_ = shrink;
  }

  void origin ( double x, double y )
  {
    origin_[ 0 ] = x;
    origin_[ 1 ] = y;
  }

  void scale ( double scale )
  {
    scale_ = scale;
  }

  template< class Entity >
  void write ( const Entity &entity, bool dash = false )
  {
    const typename Entity::Geometry &geometry = entity.geometry();

    Coordinate bary( 0.0 );
    for( int i = 0; i < geometry.corners(); ++i )
      bary += geometry.corner( i );
    bary *= ((1.0 - shrink_) / geometry.corners());

    //<polygon points="10 10 100 10 70 70" fill="white" fill-opacity="0" stroke="black" stroke-dasharray="5,2.5" />
    csvg_ << "<polygon points=\"";
    for( int i = 0; i < geometry.corners(); ++i )
    {
      Coordinate corner = geometry.corner( i );
      corner *= shrink_;
      corner += bary;
      Coordinate point = map( corner );
      csvg_ << (i > 0 ? "  " : "") << point[ 0 ] << " " << point[ 1 ];
    }
    csvg_ << "\" fill=\"white\" fill-opacity=\"0\" stroke=\"black\"";
    if( dash )
      csvg_ << " stroke-dasharray=\"5,2.5\"";
    csvg_ << "/>" << std::endl;
  }

private:
  Coordinate map( const Coordinate &coord )
  {
    Coordinate y = coord;
    y -= origin_;
    y *= scale_;
    return y;
  }
};


int main ( int argc, char **argv )
try
{
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <gridfile>" << std::endl;
    return 1;
  }

  typedef Dune::AlbertaGrid< 2 > Grid;
  typedef Grid::Codim< 0 >::LeafIterator Iterator;
  typedef Grid::Codim< 0 >::Entity Entity;
  typedef Grid::Traits::LeafIntersectionIterator IntersectionIterator;
  typedef Grid::Traits::LeafIntersection Intersection;

  Grid grid( argv[ 1 ] );
  grid.globalRefine( 1 );
  SVGWriter svg( "grid.svg" );
  svg.origin( -1.1, -1.1 );

  const Iterator end = grid.leafend< 0 >();
  for( Iterator it = grid.leafbegin< 0 >(); it != end; ++it )
  {
    const Entity &entity = *it;
    svg.shrink( 0.95 );
    svg.write( entity, false );

    if( !entity.hasBoundaryIntersections() )
      continue;

    svg.shrink( 0.9 );
    const IntersectionIterator iend = entity.ileafend();
    for( IntersectionIterator iit = entity.ileafbegin(); iit != iend; ++iit )
    {
      const Intersection &intersection = *iit;
      if( intersection.boundary() && intersection.neighbor() )
        svg.write( *(intersection.outside()), true );
    }
  }
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch( ... )
{
  std::cerr << "Unkwnown exception raised." << std::endl;
  return 2;
}
