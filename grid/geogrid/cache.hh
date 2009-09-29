// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_CACHE_HH
#define DUNE_GEOGRID_CACHE_HH

#include <dune/common/fvector.hh>

namespace Dune
{

  template< class Object >
  class GeometryGridCache
  {
    unsigned int size_;
    Object *objects_;

  public:
    GeometryGridCache ()
      : size_( 0 ),
        objects_( 0 )
    {}

    GeometryGridCache ( const GeometryGridCache &other )
      : size_( other.size_ ),
        objects_( size_ > 0 ? new Object[ size_ ] : 0 )
    {}

    ~GeometryGridCache ()
    {
      if( objects_ != 0 )
        delete[] objects_;
    }

  private:
    GeometryGridCache &operator= ( const GeometryGridCache & );

  public:
    const Object &operator[] ( unsigned int i ) const
    {
      assert( i < size_ );
      return objects_[ i ];
    }

    Object &operator[] ( unsigned int i )
    {
      assert( i < size_ );
      return objects_[ i ];
    }

    void reserve ( unsigned int size )
    {
      if( size <= size_ )
        return;

      size_ = size;
      if( objects_ != 0 )
        delete[] objects_;
      objects_ = new Object[ size_ ];
    }

    unsigned int size () const
    {
      return size_;
    }
  };

}

#endif
