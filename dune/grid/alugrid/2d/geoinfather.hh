// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_GEOINFATHER_HH
#define DUNE_ALU2DGRID_GEOINFATHER_HH

#include <dune/common/array.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>

#include <dune/grid/alugrid/2d/localgeometry.hh>

namespace Dune
{

  // ALU2DGeometryInFatherStorage
  // ----------------------------

  template< ALU2DSPACE ElementType eltype, class LocalGeometryImpl >
  class ALU2DGeometryInFatherStorage;

  template< class LocalGeometryImpl >
  class ALU2DGeometryInFatherStorage< ALU2DSPACE triangle, LocalGeometryImpl >
  {
    typedef ALU2DGeometryInFatherStorage< ALU2DSPACE triangle, LocalGeometryImpl > This;

    ALU2DGeometryInFatherStorage ( const This & );
    const This &operator= ( const This & );

  public:
    ALU2DGeometryInFatherStorage ();
    ~ALU2DGeometryInFatherStorage ();

    const LocalGeometryImpl &operator() ( int corners, bool conforming, int child ) const
    {
      assert( corners == 3 );
      assert( (child >= 0) && (child < (conforming ? 2 : 4)) );
      return *geo_[ conforming ? child : 2+child ];
    }

  private:
    LocalGeometryImpl *geo_[ 6 ];
  };

  template< class LocalGeometryImpl >
  struct ALU2DGeometryInFatherStorage< ALU2DSPACE quadrilateral, LocalGeometryImpl >
  {
    typedef ALU2DGeometryInFatherStorage< ALU2DSPACE quadrilateral, LocalGeometryImpl > This;

    ALU2DGeometryInFatherStorage ( const This & );
    const This &operator= ( const This & );

  public:
    ALU2DGeometryInFatherStorage ();
    ~ALU2DGeometryInFatherStorage ();

    const LocalGeometryImpl &operator() ( int corners, bool conforming, int child ) const
    {
      assert( corners == 4 );
      assert( !conforming );
      assert( (child >= 0) && (child < 4) );
      return *geo_[ child ];
    }

  private:
    LocalGeometryImpl *geo_[ 4 ];
  };

  template< class LocalGeometryImpl >
  struct ALU2DGeometryInFatherStorage< ALU2DSPACE mixed, LocalGeometryImpl >
  {
    typedef ALU2DGeometryInFatherStorage< ALU2DSPACE mixed, LocalGeometryImpl > This;

    ALU2DGeometryInFatherStorage ( const This & );
    const This &operator= ( const This & );

  public:
    ALU2DGeometryInFatherStorage ();
    ~ALU2DGeometryInFatherStorage ();

    const LocalGeometryImpl &operator() ( int corners, bool conforming, int child ) const
    {
      assert( (corners == 3) || (corners == 4) );
      assert( !conforming );
      assert( (child >= 0) && (child < 4) );
      return *geo_[ corners-3 ][ child ];
    }

  private:
    LocalGeometryImpl *geo_[ 2 ][ 4 ];
  };



  // Library Functions
  // -----------------

  void getALU2DGeometryInFather ( bool conforming, int child, array< FieldVector< alu2d_ctype, 2 >, 3 > &corners );
  void getALU2DGeometryInFather ( int child, array< FieldVector< alu2d_ctype, 2 >, 4 > &corners );



  // Implementation of ALU2DGeometryInFatherStorage
  // ----------------------------------------------

  template< class LocalGeometryImpl >
  inline ALU2DGeometryInFatherStorage< ALU2DSPACE triangle, LocalGeometryImpl >
  ::ALU2DGeometryInFatherStorage ()
  {
    GeometryType triangle( GenericGeometry::SimplexTopology< 2 >::type::id, 2 );
    array< FieldVector< alu2d_ctype, 2 >, 3 > corners;
    for( int child = 0; child < 2; ++child )
    {
      getALU2DGeometryInFather( true, child, corners );
      geo_[ child ] = new LocalGeometryImpl( triangle, corners );
    }
    for( int child = 0; child < 4; ++child )
    {
      getALU2DGeometryInFather( false, child, corners );
      geo_[ 2+child ] = new LocalGeometryImpl( triangle, corners );
    }
  }

  template< class LocalGeometryImpl >
  inline ALU2DGeometryInFatherStorage< ALU2DSPACE triangle, LocalGeometryImpl >
  ::~ALU2DGeometryInFatherStorage ()
  {
    for( int i = 0; i < 6; ++i )
      delete geo_[ i ];
  }


  template< class LocalGeometryImpl >
  inline ALU2DGeometryInFatherStorage< ALU2DSPACE quadrilateral, LocalGeometryImpl >
  ::ALU2DGeometryInFatherStorage ()
  {
    GeometryType quadrilateral( GenericGeometry::CubeTopology< 2 >::type::id, 2 );
    array< FieldVector< alu2d_ctype, 2 >, 4 > corners;
    for( int child = 0; child < 4; ++child )
    {
      getALU2DGeometryInFather( child, corners );
      geo_[ child ] = new LocalGeometryImpl( quadrilateral, corners );
    }
  }

  template< class LocalGeometryImpl >
  inline ALU2DGeometryInFatherStorage< ALU2DSPACE quadrilateral, LocalGeometryImpl >
  ::~ALU2DGeometryInFatherStorage ()
  {
    for( int i = 0; i < 4; ++i )
      delete geo_[ i ];
  }


  template< class LocalGeometryImpl >
  inline ALU2DGeometryInFatherStorage< ALU2DSPACE mixed, LocalGeometryImpl >
  ::ALU2DGeometryInFatherStorage ()
  {
    GeometryType triangle( GenericGeometry::SimplexTopology< 2 >::type::id, 2 );
    GeometryType quadrilateral( GenericGeometry::CubeTopology< 2 >::type::id, 2 );
    array< FieldVector< alu2d_ctype, 2 >, 3 > corners3;
    array< FieldVector< alu2d_ctype, 2 >, 4 > corners4;
    for( int child = 0; child < 4; ++child )
    {
      getALU2DGeometryInFather( false, child, corners3 );
      geo_[ 0 ][ child ] = new LocalGeometryImpl( triangle, corners3 );

      getALU2DGeometryInFather( child, corners4 );
      geo_[ 1 ][ child ] = new LocalGeometryImpl( quadrilateral, corners4 );
    }
  }

  template< class LocalGeometryImpl >
  inline ALU2DGeometryInFatherStorage< ALU2DSPACE mixed, LocalGeometryImpl >
  ::~ALU2DGeometryInFatherStorage ()
  {
    for( int i = 0; i < 4; ++i )
    {
      delete geo_[ 0 ][ i ];
      delete geo_[ 1 ][ i ];
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_ALU2DGRID_GEOINFATHER_HH
