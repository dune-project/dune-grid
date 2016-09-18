// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ALBERTAREADER_HH
#define DUNE_ALBERTA_ALBERTAREADER_HH

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/gridreader.hh>

#include <dune/grid/albertagrid/macrodata.hh>

#if HAVE_ALBERTA

namespace Dune
{

  template< class Grid >
  class AlbertaReader
      : public GridReader<Grid, AlbertaReader<Grid>>
  {
    typedef AlbertaReader< Grid > This;
    typedef GridReader<Grid, AlbertaReader<Grid>> Base;

  public:
    typedef Dune::GridFactory< Grid > GridFactory;

    typedef typename Grid::ctype ctype;

    static const int dimension = Grid::dimension;
    static const int dimensionworld = Grid::dimensionworld;

    static void read(GridFactory &factory, const std::string &fileName)
    {
      AlbertaReader<Grid> reader;
      reader.readGrid(fileName, factory);
    }

    using Base::read;

  private:
    static_assert(dimensionworld == Alberta::dimWorld,
                  "AlbertaReader: world dimension must match ALBERTA's world dimension.");

    typedef Alberta::MacroData< dimension > MacroData;

    MacroData macroData_;

    AlbertaReader ( const This & );
    This &operator= ( const This & );

//   private:
  public:
    AlbertaReader ()
    {}

    void readGrid ( const std::string &fileName, GridFactory &factory )
    {
      // read ALBERTA macro triangulation
      macroData_.read( fileName, false );

      // insert all vertices into the factory
      const int numVertices = macroData_.vertexCount();
      for( int i = 0; i < numVertices; ++i )
      {
        FieldVector< ctype, dimensionworld > v;
        const Alberta::GlobalVector &coords = macroData_.vertex( i );
        for( int j = 0; j < dimensionworld; ++j )
          v[ j ] = coords[ j ];
        factory.insertVertex( v );
      }

      // insert all elements into the factory
      std::vector< unsigned int > vertices( dimension+1 );
      const int numElements = macroData_.elementCount();
      for( int i = 0; i < numElements; ++i )
      {
        const typename MacroData::ElementId &id = macroData_.element( i );
        for( int j = 0; j <= dimension; ++j )
          vertices[ j ] = id[ j ];
        typedef typename GenericGeometry::SimplexTopology< dimension >::type Topology;
        factory.insertElement( GeometryType( Topology() ), vertices );
      }

      // release ALBERTA macro data
      macroData_.release();
    }
  };

}

#endif // #if HAVE_ALBERTA

#endif
