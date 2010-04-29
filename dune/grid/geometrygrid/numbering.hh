// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_NUMBERING_HH
#define DUNE_GEOGRID_NUMBERING_HH

namespace Dune
{

  // IdenticaldNumbering
  // -------------------

  struct IdenticalNumbering
  {
    template< int codim >
    struct EntityNumbering;

    template< class Entity >
    EntityNumbering< Entity::codimension >
    operator[] ( const Entity &entity ) const
    {
      return EntityNumbering< Entity::codimension >();
    }
  };



  // IdenticalNumbering
  // ------------------

  template< int codim >
  struct IdenticalNumbering::EntityNumbering
  {
    enum Direction { Forward, Backward };

    template< Direction direction >
    unsigned int map ( const int subcodim, const unsigned int i ) const
    {
      return i;
    }
  };

}

#endif // #ifndef DUNE_GEOGRID_NUMBERING_HH
