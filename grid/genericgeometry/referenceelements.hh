// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_REFERENCEELEMENT_HH
#define DUNE_GENERICGEOMETRY_REFERENCEELEMENT_HH

#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/subtopologies.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Topology, class ctype >
    class ReferenceElement;



    // Corner
    // ------

    template< class Topology >
    class Corner;

    template<>
    class Corner< Point >
    {
      typedef Point Topology;

      template< class > friend class Corner;
      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        assert( i < Topology :: numCorners );
      }

    public:
      enum { numCorners = Topology :: numCorners };
      enum { dimension = Topology :: dimension };

      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        evaluate_( i, x );
      }
    };

    template< class BaseTopology >
    class Corner< Prism< BaseTopology > >
    {
      typedef Prism< BaseTopology > Topology;

      template< class > friend class Corner;
      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &x )
      {
        assert( i < Topology :: numCorners );
        const unsigned int j = i % BaseTopology :: numCorners;
        Corner< BaseTopology > :: evaluate_( j, x );
        if( i >= BaseTopology :: numCorners )
          x[ dimension - 1 ] = ctype( 1 );
      }

    public:
      enum { numCorners = Topology :: numCorners };
      enum { dimension = Topology :: dimension };

      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        evaluate_( i, x );
      }
    };

    template< class BaseTopology >
    class Corner< Pyramid< BaseTopology > >
    {
      typedef Prism< BaseTopology > Topology;

      template< class > friend class Corner;
      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &x )
      {
        assert( i < Topology :: numCorners );
        if( i < BaseTopology :: numCorners )
          Corner< BaseTopology > :: evaluate_( i, x );
        else
          x[ dimension - 1 ] = ctype( 1 );
      }

    public:
      enum { numCorners = Topology :: numCorners };
      enum { dimension = Topology :: dimension };

      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        evaluate_( i, x );
      }
    };



    // IntegrationOuterNormal
    // ----------------------

    template< class Topology >
    class IntegrationOuterNormal;

    template<>
    class IntegrationOuterNormal< Point >
    {
      typedef Point Topology;

      enum { dimension = Topology :: dimension };

    public:
      enum { numNormals = 0 };

      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
      }
    };

    template< class BaseTopology >
    class IntegrationOuterNormal< Prism< BaseTopology > >
    {
      typedef Prism< BaseTopology > Topology;

      enum { dimension = Topology :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        if( i >= Size< BaseTopology, 1 > :: value )
        {
          const unsigned int j = i - Size< BaseTopology, 1 > :: value;
          n[ dimension - 1 ] = (j == 0 ? ctype( -1 ) : ctype( 1 ));
        }
        else
          IntegrationOuterNormal< BaseTopology > :: evaluate_( i, n );
      }

    public:
      enum { numNormals = Size< Topology, 1 > :: value };

      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };

    template<>
    class IntegrationOuterNormal< Prism< Point > >
    {
      typedef Prism< Point > Topology;

      enum { dimension = Topology :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        n[ dimension - 1 ] = (i == 0 ? ctype( -1 ) : ctype( 1 ));
      }

    public:
      enum { numNormals = Size< Topology, 1 > :: value };

      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };

    template< class BaseTopology >
    class IntegrationOuterNormal< Pyramid< BaseTopology > >
    {
      typedef Pyramid< BaseTopology > Topology;

      enum { dimension = Topology :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        typedef SubTopologyNumbering< BaseTopology,1,dimension-2 > Numbering;
        const unsigned int m = Size< BaseTopology, 0 > :: value;
        if( i < m )
        {
          n[ dimension -1 ] = ctype( -1 );
        }
        else
        {
          const unsigned int j = Numbering :: number( i-m, 0 );
          FieldVector< ctype, dim > x( ctype( 0 ) );
          Corner< BaseTopology > :: evaluate_( j, x );

          IntegrationOuterNormal< BaseTopology > :: evaluate_( i-m, n );
          n[ dimension - 1 ] = (x * n);
        }
      }

    public:
      enum { numNormals = Size< Topology, 1 > :: value };

      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };

    template<>
    class IntegrationOuterNormal< Pyramid< Point > >
    {
      typedef Pyramid< Point > Topology;

      enum { dimension = Topology :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        n[ dimension - 1 ] = (i == 0 ? ctype( -1 ) : ctype( 1 ));
      }

    public:
      enum { numNormals = Size< Topology, 1 > :: value };

      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };



    // Volume
    // ------

    template< class Topology >
    struct Volume;

    template<>
    struct Volume< Point >
    {
      template< class ctype >
      static ctype evaluate ()
      {
        return ctype( 1 );
      }
    };

    template< class BaseTopology >
    struct Volume< Prism< BaseTopology > >
    {
      template< class ctype >
      static ctype evaluate ()
      {
        return Volume< BaseTopology > :: template evaluate< ctype >();
      }
    };

    template< class BaseTopology >
    struct Volume< Pyramid< BaseTopology > >
    {
      typedef Pyramid< BaseTopology > Topology;
      template< class ctype >
      static ctype evaluate ()
      {
        const ctype baseVolume
          = Volume< BaseTopology > :: template evaluate< ctype >();
        return baseVolume / ctype( Topology::dimension );
      }
    };



    // ReferenceElement
    // ----------------

    template< class Topology, class ctype >
    struct ReferenceElement
    {
      enum { dimension = Topology :: dimension };

      enum { numCorners = Topology :: numCorners };
      enum { numNormals = IntegrationOuterNormal< Topology > :: numNormals };

      typedef FieldVector< ctype, dimension > CoordinateType;

      template< unsigned int codim >
      struct Codim
      {
        enum { size = Size< Topology, codim > :: value };
      };

      template< unsigned int codim, unsigned int subcodim >
      static unsigned int subNumbering ( unsigned int i, unsigned int j )
      {
        return SubTopologyNumbering< Topology, codim, subcodim > :: number( i, j );
      }

      template< unsigned int codim, unsigned int subcodim >
      static unsigned int size ( unsigned int i )
      {
        return SubTopologySize< Topology, codim, subcodim > :: size( i );
      }

      template< unsigned int codim >
      static const FieldVector< ctype, dimension > &
      baryCenter ( unsigned int i )
      {
        Int2Type< codim > codimVariable;
        return instance().baryCenters_[ codimVariable ][ i ];
      }

      static const CoordinateType &corner ( unsigned int i )
      {
        assert( i < numCorners );
        return instance().corners_[ i ];
      }

      static const CoordinateType &
      integrationOuterNormal ( unsigned int i )
      {
        assert( i < numNormals );
        return instance().normals_[ i ];
      }

      static ctype volume ()
      {
        return Volume< Topology > :: template evaluate< ctype >();
      }

    private:
      template< int codim >
      class BaryCenterArray;

      ReferenceElement ()
      {
        for( unsigned int i = 0; i < numCorners; ++i )
          Corner< Topology > :: evaluate( i, corners_[ i ] );
        for( unsigned int i = 0; i < numNormals; ++i )
          IntegrationOuterNormal< Topology > :: evaluate( i, normals_[ i ] );
      }

      static const ReferenceElement &instance ()
      {
        static ReferenceElement inst;
        return inst;
      }

      CoordinateType corners_[ numCorners ];
      CodimTable< BaryCenterArray, dimension > baryCenters_;
      CoordinateType normals_[ numNormals ];
    };



    template< class Topology, class ctype >
    template< int codim >
    class ReferenceElement< Topology, ctype > :: BaryCenterArray
    {
      enum { Size = GenericGeometry :: Size< Topology, codim > :: value };

      typedef FieldVector< ctype, dimension > CoordinateType;

      template< int i >
      struct Builder;

      CoordinateType baryCenters_[ Size ];

    public:
      BaryCenterArray ()
      {
        ForLoop< Builder, 0, Size-1 > :: apply( baryCenters_ );
      }

      const CoordinateType &operator[] ( unsigned int i ) const
      {
        assert( i < Size );
        return baryCenters_[ i ];
      }

      static unsigned int size ()
      {
        return Size;
      }
    };

    template< class Topology, class ctype >
    template< int codim >
    template< int i >
    struct ReferenceElement< Topology, ctype > :: BaryCenterArray< codim > :: Builder
    {
      static void apply ( CoordinateType (&baryCenters)[ Size ] )
      {
        typedef SubTopologyNumbering< Topology, codim, dimension - codim > Numbering;

        CoordinateType &x = baryCenters[ i ];
        Corner< Topology > :: evaluate( 0, x );

        const unsigned int numCorners = Corner< Topology > :: numCorners;
        for( unsigned int k = 1; k < numCorners; ++k )
        {
          unsigned int j = Numbering :: number( i, k );

          CoordinateType y;
          Corner< Topology > :: evaluate( j, y );
          x += y;
        }
        x *= ctype( 1 ) / ctype( numCorners );
      }
    };

  }



  template< class ctype, int dim >
  class GenericReferenceElement
  {
    class SubEntityInfo;
    template< int codim > struct Initialize;

    std :: vector< SubEntityInfo > info_[ dim+1 ];
    double volume_;

  public:
    int size ( int c ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ].size();
    }

    int size ( int i, int c, int cc )
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ][ i ].size( cc );
    }

    int subEntity ( int i, int c, int ii, int cc ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ][ i ].number( ii, cc );
    }

    const FieldVector< ctype, dim > &position( int i, int c ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ][ i ].position();
    }

    template< int codim >
    FieldVector< ctype, dim >
    global( const FieldVector< ctype, dim-codim > &local, int i, int c ) const
    {
      // todo: implement this function
      return position( i, c );
    }

    GeometryType type ( int i, int c ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ][ i ].type();
    }

    double volume () const
    {
      return volume_;
    }

    template< class Topology >
    void initialize ( const GenericGeometry :: ReferenceElement< Topology, ctype > &refElement )
    {
      GenericGeometry :: ForLoop< Initialize, 0, dim > :: apply( refElement, info_ );
      volume_ = GenericGeometry :: Volume< Topology > :: template evaluate< double >();
    }
  };

  template< class ctype, int dim >
  class GenericReferenceElement< ctype, dim > :: SubEntityInfo
  {
    template< int codim > struct Initialize
    {
      template< int subcodim > struct SubCodim;
    };

    int codim_;
    std :: vector< int > numbering_[ dim+1 ];
    FieldVector< ctype, dim > baryCenter_;
    GeometryType type_;

  public:
    int size ( int cc ) const
    {
      assert( (cc >= codim_) && (cc <= dim) );
      return numbering_[ cc ].size();
    }

    int number ( int ii, int cc ) const
    {
      assert( (cc >= codim_) && (cc <= dim) );
      return numbering_[ cc ][ ii ];
    }

    const FieldVector< ctype, dim > &position () const
    {
      return baryCenter_;
    }

    GeometryType type () const
    {
      return type_;
    }

    template< class Topology, int codim >
    void
    initialize ( const GenericGeometry :: ReferenceElement< Topology, ctype > &refElement,
                 unsigned int i )
    {
      GenericGeometry :: ForLoop< Initialize< codim > :: template SubCodim, 0, dim-codim >
      :: apply( refElement, i, numbering_ );
      baryCenter_ = refElement.template baryCenter< codim >( i );
      type_ = GenericGeometry :: DuneGeometryType< Topology, GeometryType :: simplex > :: type();
    }
  };

  template< class ctype, int dim >
  template< int codim >
  template< int subcodim >
  struct GenericReferenceElement< ctype, dim > :: SubEntityInfo :: Initialize< codim > :: SubCodim
  {
    template< class Topology >
    static void
    apply ( const GenericGeometry :: ReferenceElement< Topology, ctype > &refElement,
            int i,
            std :: vector< int > (&numbering)[ dim+1 ] )
    {
      const unsigned int size = refElement.template size< codim, subcodim >( i );
      numbering[ codim+subcodim ].resize( size );
      for( unsigned int j = 0; j < size; ++j )
      {
        numbering[ codim+subcodim ][ j ]
          = refElement.template subNumbering< codim, subcodim >[ i ][ j ];
      }
    }
  };

  template< class ctype, int dim >
  template< int codim >
  struct GenericReferenceElement< ctype, dim > :: Initialize
  {
    template< class Topology >
    static void
    apply ( const GenericGeometry :: ReferenceElement< Topology, ctype > &refElement,
            std :: vector< SubEntityInfo > (&info)[ dim+1 ] )
    {
      const unsigned int size = GenericGeometry :: Size< Topology, codim > :: value;
      info[ codim ].resize( size );
      for( unsigned int i = 0; i < size; ++i )
        info[ codim ][ i ].template initialize< Topology, codim >( refElement, i );
    }
  };

}

#endif
