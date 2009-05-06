// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_REFERENCETOPOLOGIES_HH
#define DUNE_GENERICGEOMETRY_REFERENCETOPOLOGIES_HH

#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/subtopologies.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // ReferenceTopology
    // -----------------

    template< unsigned int dim >
    class ReferenceTopology
    {
      typedef ReferenceTopology< dim > This;

      class SubEntityInfo;
      template< class Topology > struct Initialize;

    public:
      static const unsigned int dimension = dim;

      unsigned int size ( unsigned int codim ) const
      {
        assert( codim <= dimension );
        return info_[ codim ].size();
      }

      unsigned int
      size ( unsigned int codim, unsigned int i, unsigned int subcodim ) const
      {
        assert( (codim <= dimension) && (i < info_[ codim ].size()) );
        return info_[ codim ][ i ].size( subcodim );
      }

      unsigned int subEntity ( unsigned int codim, unsigned int i,
                               unsigned int subcodim, unsigned int j ) const
      {
        assert( (codim <= dimension) && (i < info_[ codim ].size()) );
        return info_[ codim ][ i ].number( subcodim, j );
      }

      unsigned int topologyId ( unsigned int codim, unsigned int i ) const
      {
        assert( (codim <= dimension) && (i < info_[ codim ].size()) );
        return info_[ codim ][ i ].topologyId();
      }

      const GeometryType &type ( unsigned int codim, unsigned int i ) const
      {
        assert( (codim <= dimension) && (i < info_[ codim ].size()) );
        return info_[ codim ][ i ].type();
      }

      template< class Topology >
      void initialize ()
      {
        typedef Initialize< Topology > Init;
        ForLoop< Init::template Codim, 0, dimension >::apply( info_ );
      }

    private:
      std::vector< SubEntityInfo > info_[ dimension+1 ];
    };



    // ReferenceTopology::SubEntityInfo
    // --------------------------------

    template< unsigned int dim >
    class ReferenceTopology< dim >::SubEntityInfo
    {
      template< class Topology, unsigned int codim > struct Initialize
      {
        template< int subcodim > struct SubCodim;
      };

    public:
      unsigned int size ( unsigned int subcodim ) const
      {
        return numbering_[ subcodim ].size();
      }

      unsigned int number ( unsigned int subcodim, unsigned int j ) const
      {
        return numbering_[ subcodim ][ j ];
      }

      unsigned int topologyId () const
      {
        return topologyId_;
      }

      const GeometryType &type () const
      {
        return type_;
      }

      template< class Topology, unsigned int codim, unsigned int i >
      void initialize ()
      {
        typedef Initialize< Topology, codim > Init;
        typedef typename GenericGeometry::SubTopology< Topology, codim, i >::type SubTopology;

        codim_ = codim;
        topologyId_ = SubTopology::id;
        type_ = DuneGeometryType< SubTopology, GeometryType::simplex >::type();
        numbering_.resize( SubTopology::dimension+1 );

        const unsigned int iVariable = i;
        ForLoop< Init::template SubCodim, 0, SubTopology::dimension >::apply( iVariable, numbering_ );
      }

    private:
      int codim_;
      unsigned int topologyId_;
      GeometryType type_;
      std::vector< std::vector< unsigned int > > numbering_;
    };


    template< unsigned int dim >
    template< class Topology, unsigned int codim >
    template< int subcodim >
    struct ReferenceTopology< dim >::SubEntityInfo::Initialize< Topology, codim >::SubCodim
    {
      typedef SubTopologySize< Topology, codim, subcodim > Size;
      typedef SubTopologyNumbering< Topology, codim, subcodim > Numbering;

      static void
      apply ( unsigned int i, std::vector< std::vector< unsigned int > > &numbering )
      {
        const unsigned int size = Size::size( i );
        numbering[ subcodim ].resize( size );
        for( unsigned int j = 0; j < size; ++j )
          numbering[ subcodim ][ j ] = Numbering::number( i, j );
      }
    };


    // ReferenceTopology::Initialize
    // -----------------------------

    template< unsigned int dim >
    template< class Topology >
    struct ReferenceTopology< dim >::Initialize
    {
      template< int codim >
      struct Codim
      {
        template< int i >
        struct SubTopology
        {
          static void apply ( std::vector< SubEntityInfo > &info )
          {
            info[ i ].template initialize< Topology, codim, i >();
          }
        };

        static void apply ( std::vector< SubEntityInfo > (&info)[ dim+1 ] )
        {
          const unsigned int size = Size< Topology, codim >::value;
          info[ codim ].resize( size );
          ForLoop< SubTopology, 0, size-1 >::apply( info[ codim ] );
        }
      };
    };


    // ReferenceTopologyContainer
    // --------------------------

    template< unsigned int dim >
    class ReferenceTopologies
    {
      typedef ReferenceTopologies< dim > This;

      template< int topologyId >
      struct Init;

    public:
      static const unsigned int dimension = dim;
      static const unsigned int numTopologies = (1 << dimension);

      typedef GenericGeometry::ReferenceTopology< dimension > ReferenceTopology;

      static const ReferenceTopology &get ( const unsigned int topologyId )
      {
        assert( topologyId < numTopologies );
        return instance().refTopology_[ topologyId ];
      }

    private:
      ReferenceTopologies ()
      {
        ForLoop< Init, 0, numTopologies-1 >::apply( refTopology_ );
      }

      ReferenceTopologies ( const This & );
      This &operator= ( const This & );

      static const This &instance ()
      {
        static This instance;
        return instance;
      }

      ReferenceTopology refTopology_[ numTopologies ];
    };


    template< unsigned int dim >
    template< int topologyId >
    struct ReferenceTopologies< dim >::Init
    {
      static void apply ( ReferenceTopology (&refTopology)[ numTopologies ] )
      {
        typedef typename GenericGeometry::Topology< topologyId, dimension >::type Topology;
        refTopology[ topologyId ].template initialize< Topology >();
      }
    };

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_REFERENCETOPOLOGIES_HH
