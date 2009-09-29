// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_STORAGE_HH
#define DUNE_GEOGRID_STORAGE_HH

namespace Dune
{

  template< class Object >
  class GeometryGridStorage
  {
    struct MemObject
    {
      Object object;
      MemObject *next;
    };

    MemObject *top_;

    GeometryGridStorage ()
      : top_( 0 )
    {}

    GeometryGridStorage ( const GeometryGridStorage & );

    ~GeometryGridStorage ()
    {
      while( top_ != 0 )
      {
        MemObject *ptr = top_;
        top_ = top_->next;
        delete ptr;
      }
    }

    static GeometryGridStorage &instance ()
    {
      static GeometryGridStorage singleton;
      return singleton;
    }

  public:
    static Object *alloc ()
    {
      MemObject *&top = instance().top_;

      MemObject *ptr = top;
      if( ptr != 0 )
        top = ptr->next;
      else
        ptr = new MemObject;
      return &(ptr->object);
    }

    static void free ( Object *object )
    {
      MemObject *&top = instance().top_;

      MemObject *ptr = reinterpret_cast< MemObject * >( object );
      if( ptr != 0 )
      {
        ptr->next = top;
        top = ptr;
      }
    }
  };

}

#endif
