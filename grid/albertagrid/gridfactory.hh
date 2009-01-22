// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GRIDFACTORY_HH
#define DUNE_ALBERTA_GRIDFACTORY_HH

#include <dune/grid/common/gridfactory.hh>

#include <dune/grid/albertagrid/agrid.hh>


namespace Dune
{

  template< int dim, int dimworld >
  class GridFactory< AlbertaGrid< dim, dimworld > >
    : public GridFactoryInterface< AlbertaGrid< dim, dimworld > >
  {
    typedef GridFactory< AlbertaGrid< dim, dimworld > > This;

  public:
    typedef AlbertaGrid< dim, dimworld > Grid;

    typedef typename Grid::ctype ctype;

    static const int dimension = Grid::dimension;
    static const int dimensionworld = Grid::dimensionworld;

    static const bool supportsBoundaryIds = (DUNE_ALBERTA_VERSION >= 0x200);

    typedef FieldVector< ctype, dimensionworld > Coordinate;

  private:
    static const int numVertices
      = Alberta::NumSubEntities< dimension, dimension >::value;

    typedef Alberta::MacroData< dimension > MacroData;
    typedef Alberta::NumberingMap< dimension > NumberingMap;

    MacroData macroData_;
    NumberingMap numberingMap_;

  public:
    GridFactory ()
    {
      macroData_.create();
    }

    ~GridFactory ()
    {
      macroData_.release();
    }

    virtual void insertVertex ( const Coordinate &coord )
    {
      macroData_.insertVertex( coord );
    }

    virtual void insertElement ( const GeometryType &type,
                                 const std::vector< unsigned int > &vertices )
    {
      if( (int)type.dim() != dimension )
        DUNE_THROW( AlbertaError, "Inserting element of wrong dimension: " << type.dim() );
      if( !type.isSimplex() )
        DUNE_THROW( AlbertaError, "Alberta supports only simplices." );

      if( vertices.size() != (size_t)numVertices )
        DUNE_THROW( AlbertaError, "Wrong number of vertices passed: " << vertices.size() << "." );

      int array[ numVertices ];
      for( int i = 0; i < numVertices; ++i )
        array[ i ] = vertices[ numberingMap_.alberta2dune( dimension, i ) ];
      macroData_.insertElement( array );
    }

    virtual void insertBoundary ( int element, int face, int id )
    {
      if( (id <= 0) || (id > 127) )
        DUNE_THROW( AlbertaError, "Invalid boundary id: " << id << "." );
      macroData_.boundaryId( element, numberingMap_.dune2alberta( 1, face ) ) = id;
    }

    Grid *createGrid ( const std::string &gridName, bool markLongestEdge = false )
    {
      macroData_.finalize();
      if( markLongestEdge )
        macroData_.markLongestEdge();
      return new Grid( macroData_, gridName );
    }

    virtual Grid *createGrid ()
    {
      return createGrid( "AlbertaGrid", false );
    }

    static void destroyGrid ( Grid *grid )
    {
      delete grid;
    }

    virtual bool write ( const std::string &filename, bool binary = false )
    {
      macroData_.finalize();
      return macroData_.write( filename, binary );
    }
  };

}

#endif
