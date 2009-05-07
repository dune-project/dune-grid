// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICREFERENCEELEMENTS_HH
#define DUNE_GENERICREFERENCEELEMENTS_HH

// #include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/geometry.hh>

namespace Dune
{

  // GenericReferenceElement
  // -----------------------

  template< class ctype, int dim >
  class GenericReferenceElement
  {
    typedef GenericReferenceElement< ctype, dim > This;

    class SubEntityInfo;
    template< class Topology > class CornerStorage;
    template< class Topology > struct Initialize;

    struct GeometryTraits
      : public GenericGeometry::DefaultGeometryTraits< ctype, dim, dim >
    {
      typedef GenericGeometry::DefaultGeometryTraits< ctype, dim, dim > Base;

      typedef typename Base::CoordTraits CoordTraits;

      template< class Topology >
      struct Mapping
      {
        typedef GenericGeometry::CornerMapping< CoordTraits, Topology, dim, CornerStorage< Topology >, true > type;
      };

      struct Caching
      {
        static const GenericGeometry::EvaluationType evaluateJacobianTransposed = GenericGeometry::PreCompute;
        static const GenericGeometry::EvaluationType evaluateJacobianInverseTransposed = GenericGeometry::PreCompute;
        static const GenericGeometry::EvaluationType evaluateIntegrationElement = GenericGeometry::PreCompute;
        static const GenericGeometry::EvaluationType evaluateNormal = GenericGeometry::PreCompute;
      };
    };

  public:
    template< int codim >
    struct Codim
    {
      typedef GenericGeometry::HybridMapping< dim-codim, GeometryTraits > Mapping;
    };

  private:
    template< int codim >
    class MappingArray
      : public std::vector< typename Codim< codim >::Mapping * >
    {};

    typedef GenericGeometry::CodimTable< MappingArray, dim > MappingsTable;

    std::vector< SubEntityInfo > info_[ dim+1 ];
    double volume_;
    MappingsTable mappings_;

  public:
    int size ( int c ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ].size();
    }

    int size ( int i, int c, int cc ) const
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
      if( c != codim )
        DUNE_THROW( Exception, "Local Coordinate Type does not correspond to codimension c." );
      assert( c == codim );
      return mapping< codim >( i ).global( local );
    }

    template< int codim >
    FieldVector< ctype, dim >
    global( const FieldVector< ctype, dim-codim > &local, int i ) const
    {
      return mapping< codim >( i ).global( local );
    }

    template< int codim >
    typename Codim< codim >::Mapping &mapping( int i ) const
    {
      Int2Type< codim > codimVariable;
      return *(mappings_[ codimVariable ][ i ]);
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

    template< GeometryType::BasicType geoType >
    void initialize ()
    {
      typedef typename GenericGeometry::Convert< geoType, dim >::type Topology;
      typedef Initialize< Topology > Init;
      typedef GenericGeometry::VirtualMapping< Topology, GeometryTraits > VirtualMapping;

      Int2Type< 0 > codim0Variable;
      mappings_[ codim0Variable ].resize( 1 );
      mappings_[ codim0Variable ][ 0 ]  = new VirtualMapping( codim0Variable );

      GenericGeometry::ForLoop< Init::template Codim, 0, dim >::apply( info_, mappings_ );
      volume_ = GenericGeometry::ReferenceDomain< Topology >::template volume< double >();
    }
  };


  template< class ctype, int dim >
  class GenericReferenceElement< ctype, dim >::SubEntityInfo
  {
    template< class Topology, int codim > struct Initialize
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

    template< class Topology, unsigned int codim, unsigned int i >
    void initialize ()
    {
      typedef Initialize< Topology, codim > Init;
      typedef GenericGeometry::ReferenceElement< Topology, ctype > RefElement;

      const unsigned int iVariable = i;
      GenericGeometry::ForLoop< Init::template SubCodim, 0, dim-codim >::apply( iVariable, numbering_ );
      baryCenter_ = RefElement::template baryCenter< codim >( i );

      typedef typename GenericGeometry::SubTopology< Topology, codim, i >::type SubTopology;
      type_ = GenericGeometry::DuneGeometryType< SubTopology, GeometryType::simplex >::type();
    }
  };


  template< class ctype, int dim >
  template< class Topology >
  class GenericReferenceElement< ctype, dim >::CornerStorage
  {
    typedef GenericGeometry::ReferenceElement< Topology, ctype > RefElement;

  public:
    static const unsigned int size = Topology::numCorners;

    template< class SubTopology >
    struct SubStorage
    {
      typedef CornerStorage< SubTopology > type;
    };

    explicit CornerStorage ( const Int2Type< 0 > & )
    {
      for( unsigned int i = 0; i < size; ++i )
        coords_[ i ] = RefElement::corner( i );
    }

    template< class Mapping, unsigned int codim >
    explicit
    CornerStorage ( const GenericGeometry::SubMappingCoords< Mapping, codim > &coords )
    {
      for( unsigned int i = 0; i < size; ++i )
        coords_[ i ] = coords[ i ];
    }

    const FieldVector< ctype, dim > &operator[] ( unsigned int i ) const
    {
      return coords_[ i ];
    }

  private:
    FieldVector< ctype, dim > coords_[ size ];
  };


  template< class ctype, int dim >
  template< class Topology, int codim >
  template< int subcodim >
  struct GenericReferenceElement< ctype, dim >::SubEntityInfo::Initialize< Topology, codim >::SubCodim
  {
    typedef GenericGeometry::ReferenceElement< Topology, ctype > RefElement;

    static void apply ( unsigned int i, std::vector< int > (&numbering)[ dim+1 ] )
    {
      const unsigned int size = RefElement::template size< codim, subcodim >( i );
      numbering[ codim+subcodim ].resize( size );
      for( unsigned int j = 0; j < size; ++j )
        numbering[ codim+subcodim ][ j ] = RefElement::template subNumbering< codim, subcodim >( i, j );
    }
  };


  template< class ctype, int dim >
  template< class Topology >
  struct GenericReferenceElement< ctype, dim >::Initialize
  {
    typedef Dune::GenericReferenceElement< ctype, dim > GenericReferenceElement;

    typedef typename GenericReferenceElement::template Codim< 0 >::Mapping ReferenceMapping;

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

      static void apply ( std::vector< SubEntityInfo > (&info)[ dim+1 ],
                          MappingsTable &mappings )
      {
        const unsigned int size = GenericGeometry::Size< Topology, codim >::value;
        info[ codim ].resize( size );
        GenericGeometry::ForLoop< SubTopology, 0, size-1 >::apply( info[ codim ] );
        /*
           for( unsigned int i = 0; i < size; ++i )
           info[ codim ][ i ].template initialize< Topology, codim >( i );
         */

        if( codim > 0 )
        {
          Int2Type< 0 > codim0Variable;
          const ReferenceMapping &refMapping = *(mappings[ codim0Variable ][ 0 ]);

          Int2Type< codim > codimVariable;
          mappings[ codimVariable ].resize( size );
          for( unsigned int i = 0; i < size; ++i )
            mappings[ codimVariable ][ i ] = refMapping.template trace< codim >( i );
        }
      }
    };
  };



  // GenericReferenceElementContainer
  // --------------------------------

  template< class ctype, int dim >
  struct GenericReferenceElementContainer
  {
    typedef GenericReferenceElement< ctype, dim > value_type;

    const value_type &operator() ( const GeometryType &type ) const
    {
      assert( type.dim() == dim );
      switch( type.basicType() )
      {
      case GeometryType::simplex :
        return simplex_;

      case GeometryType::cube :
        return cube_;

      case GeometryType::pyramid :
        return pyramid_;

      case GeometryType::prism :
        return prism_;

      default :
        DUNE_THROW( RangeError, "Unknown geometry type: " << type );
      }
    }

    const value_type &simplex () const
    {
      return simplex_;
    }

    const value_type &cube () const
    {
      return cube_;
    }

    static const GenericReferenceElementContainer &instance ()
    {
      static GenericReferenceElementContainer inst;
      return inst;
    }

  private:
    GenericReferenceElementContainer ()
    {
      simplex_.template initialize< GeometryType::simplex >();
      cube_.template initialize< GeometryType::cube >();
      pyramid_.template initialize< GeometryType::pyramid >();
      prism_.template initialize< GeometryType::prism >();
    }

    value_type simplex_;
    value_type cube_;
    value_type pyramid_;
    value_type prism_;
  };

  template< class ctype >
  struct GenericReferenceElementContainer< ctype, 2 >
  {
    typedef GenericReferenceElement< ctype, 2 > value_type;

    const value_type &operator() ( const GeometryType &type ) const
    {
      assert( type.dim() == 2 );
      switch( type.basicType() )
      {
      case GeometryType::simplex :
        return simplex_;

      case GeometryType::cube :
        return cube_;

      case GeometryType::pyramid :
      case GeometryType::prism :
        DUNE_THROW( RangeError, "Invalid geometry type: " << type );

      default :
        DUNE_THROW( RangeError, "Unknown geometry type: " << type );
      }
    }

    const value_type &simplex () const
    {
      return simplex_;
    }

    const value_type &cube () const
    {
      return cube_;
    }

    static const GenericReferenceElementContainer &instance ()
    {
      static GenericReferenceElementContainer inst;
      return inst;
    }

  private:
    GenericReferenceElementContainer ()
    {
      simplex_.template initialize< GeometryType::simplex >();
      cube_.template initialize< GeometryType::cube >();
    }

    value_type simplex_;
    value_type cube_;
  };

  template< class ctype >
  struct GenericReferenceElementContainer< ctype, 1 >
  {
    typedef GenericReferenceElement< ctype, 1 > value_type;

    const value_type &operator() ( const GeometryType &type ) const
    {
      assert( type.dim() == 1 );
      return line_;
    }

    const value_type &simplex () const
    {
      return line_;
    }

    const value_type &cube () const
    {
      return line_;
    }

    static const GenericReferenceElementContainer &instance ()
    {
      static GenericReferenceElementContainer inst;
      return inst;
    }

  private:
    GenericReferenceElementContainer ()
    {
      line_.template initialize< GeometryType::simplex >();
    }

    value_type line_;
  };

  template< class ctype >
  struct GenericReferenceElementContainer< ctype, 0 >
  {
    typedef GenericReferenceElement< ctype, 0 > value_type;

    const value_type &operator() ( const GeometryType &type ) const
    {
      assert( type.dim() == 0 );
      return point_;
    }

    const value_type &simplex () const
    {
      return point_;
    }

    const value_type &cube () const
    {
      return point_;
    }

    static
    const GenericReferenceElementContainer & instance()
    {
      static GenericReferenceElementContainer inst;
      return inst;
    }

  private:
    GenericReferenceElementContainer ()
    {
      point_.template initialize< GeometryType::simplex >();
    }

    value_type point_;
  };



  // GenericReferenceElements
  // ------------------------

  template< class ctype, int dim >
  struct GenericReferenceElements
  {
    static const GenericReferenceElement< ctype, dim > &
    general ( const GeometryType &type )
    {
      return GenericReferenceElementContainer< ctype, dim >::instance() ( type );
    }

    static const GenericReferenceElement< ctype, dim > &simplex ()
    {
      return GenericReferenceElementContainer< ctype, dim >::instance().simplex();
    }

    static const GenericReferenceElement< ctype, dim > &cube ()
    {
      return GenericReferenceElementContainer< ctype, dim >::instance().cube();
    }
  };

}

#endif
