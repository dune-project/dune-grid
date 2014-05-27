// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ALBERTAREADER_HH
#define DUNE_ALBERTA_ALBERTAREADER_HH

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/grid/utility/grapedataioformattypes.hh>

#include <dune/grid/albertagrid/macrodata.hh>

#if HAVE_ALBERTA

namespace Dune
{

  template< class Grid >
  class AlbertaReader
  {
    typedef AlbertaReader< Grid > This;

  public:
    typedef Dune::GridFactory< Grid > GridFactory;

    typedef typename Grid::ctype ctype;

    static const int dimension = Grid::dimension;
    static const int dimensionworld = Grid::dimensionworld;

  private:
    static_assert(dimensionworld == Alberta::dimWorld,
                  "AlbertaReader: world dimension must match ALBERTA's world dimension.");

    typedef Alberta::MacroData< dimension > MacroData;

    MacroData macroData_;

    AlbertaReader ( const This & );
    This &operator= ( const This & );

  public:
    AlbertaReader ()
    {}

    template< GrapeIOFileFormatType type >
    void readGrid ( const std::string &fileName, GridFactory &factory )
    {
      static_assert(type != pgm, "AlbertaReader: reading pgm format is not supported.");

      // read ALBERTA macro triangulation
      macroData_.read( fileName, (type == xdr) );

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

    void readGrid ( const std::string &filename, GridFactory &factory )
    {
      readGrid< ascii >( filename, factory );
    }
  };

}

#endif // #if HAVE_ALBERTA

#endif
