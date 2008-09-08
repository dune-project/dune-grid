// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GENERICREFERENCEELEMENT_HH
#define DUNE_GENERICGEOMETRY_GENERICREFERENCEELEMENT_HH

#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>

namespace Dune
{

  // GenericReferenceElement
  // -----------------------

  template< class ctype, int dim >
  class GenericReferenceElement
  {
    class SubEntityInfo;
    template< class Topology > struct Initialize;

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

    template< GeometryType :: BasicType type >
    void initialize ()
    {
      typedef typename GenericGeometry :: Convert< type, dim > :: type Topology;
      typedef Initialize< Topology > Init;
      GenericGeometry :: ForLoop< Init :: template Codim, 0, dim > :: apply( info_ );
      volume_ = GenericGeometry :: Volume< Topology > :: template evaluate< double >();
    }
  };

  template< class ctype, int dim >
  class GenericReferenceElement< ctype, dim > :: SubEntityInfo
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

    template< class Topology, int codim >
    void initialize ( unsigned int i )
    {
      typedef Initialize< Topology, codim > Init;
      typedef GenericGeometry :: ReferenceElement< Topology, ctype > RefElement;

      GenericGeometry :: ForLoop< Init :: template SubCodim, 0, dim-codim > :: apply( i, numbering_ );
      baryCenter_ = RefElement :: template baryCenter< codim >( i );
      type_ = GenericGeometry :: DuneGeometryType< Topology, GeometryType :: simplex > :: type();
    }
  };

  template< class ctype, int dim >
  template< class Topology, int codim >
  template< int subcodim >
  struct GenericReferenceElement< ctype, dim >
  :: SubEntityInfo :: Initialize< Topology, codim > :: SubCodim
  {
    typedef GenericGeometry :: ReferenceElement< Topology, ctype > RefElement;

    static void apply ( unsigned int i, std :: vector< int > (&numbering)[ dim+1 ] )
    {
      const unsigned int size = RefElement :: template size< codim, subcodim >( i );
      numbering[ codim+subcodim ].resize( size );
      for( unsigned int j = 0; j < size; ++j )
      {
        numbering[ codim+subcodim ][ j ]
          = RefElement ::  template subNumbering< codim, subcodim >[ i ][ j ];
      }
    }
  };

  template< class ctype, int dim >
  template< class Topology >
  struct GenericReferenceElement< ctype, dim > :: Initialize
  {
    template< int codim >
    struct Codim
    {
      static void apply ( std :: vector< SubEntityInfo > (&info)[ dim+1 ] )
      {
        const unsigned int size = GenericGeometry :: Size< Topology, codim > :: value;
        info[ codim ].resize( size );
        for( unsigned int i = 0; i < size; ++i )
          info[ codim ][ i ].template initialize< Topology, codim >( i );
      }
    };
  };



  // GenericReferenceElementContainer
  // --------------------------------

  template< class ctype, int dim >
  struct GenericReferenceElementContainer
  {
    typedef GenericReferenceElement< ctype, dim > value_type;

  private:
    value_type simplex;
    value_type cube;
    value_type pyramid;
    value_type prism;

    GenericReferenceElementContainer ()
    {
      simplex.template initialize< GeometryType :: simplex >();
      cube.template initialize< GeometryType :: cube >();
      pyramid.template initialize< GeometryType :: pyramid >();
      prism.template initialize< GeometryType :: prism >();
    }

  public:
    const value_type &operator() ( const GeometryType &type ) const
    {
      assert( type.dim() == dim );
      switch( type.basicType() )
      {
      case GeometryType :: simplex :
        return simplex;

      case GeometryType :: cube :
        return cube;

      case GeometryType :: pyramid :
        return pyramid;

      case GeometryType :: prism :
        return prism;

      default :
        DUNE_THROW( RangeError, "Unknown geometry type: " << type );
      }
    }
  };

  template< class ctype >
  struct GenericReferenceElementContainer< ctype, 2 >
  {
    typedef GenericReferenceElement< ctype, 2 > value_type;

  private:
    value_type simplex;
    value_type cube;

    GenericReferenceElementContainer ()
    {
      simplex.template initialize< GeometryType :: simplex >();
      cube.template initialize< GeometryType :: cube >();
    }

  public:
    const value_type &operator() ( const GeometryType &type ) const
    {
      assert( type.dim() == 2 );
      switch( type.basicType() )
      {
      case GeometryType :: simplex :
        return simplex;

      case GeometryType :: cube :
        return cube;

      case GeometryType :: pyramid :
      case GeometryType :: prism :
        DUNE_THROW( RangeError, "Invalid geometry type: " << type );

      default :
        DUNE_THROW( RangeError, "Unknown geometry type: " << type );
      }
    }
  };

  template< class ctype >
  struct GenericReferenceElementContainer< ctype, 1 >
  {
    typedef GenericReferenceElement< ctype, 1 > value_type;

  private:
    value_type line;

    GenericReferenceElementContainer ()
    {
      line.template initialize< GeometryType :: simplex >();
    }

  public:
    const value_type &operator() ( const GeometryType &type ) const
    {
      assert( type.dim() == 1 );
      return line;
    }
  };

  template< class ctype >
  struct GenericReferenceElementContainer< ctype, 0 >
  {
    typedef GenericReferenceElement< ctype, 0 > value_type;

  private:
    value_type point;

    GenericReferenceElementContainer ()
    {
      point.template initialize< GeometryType :: simplex >();
    }

  public:
    const value_type &operator() ( const GeometryType &type ) const
    {
      assert( type.dim() == 0 );
      return point;
    }
  };

}

#endif
