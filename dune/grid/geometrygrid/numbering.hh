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

    struct IntersectionNumbering;

    template< int codim, int dim, class G, template< int, int, class > class I >
    EntityNumbering< codim > operator[] ( const Entity< codim, dim, G, I > &entity ) const;

    template< class G, template< class > class I >
    IntersectionNumbering operator[] ( const Intersection< G, I > &intersection ) const;
  };



  // IdenticalNumbering::EntityNumbering
  // -----------------------------------

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



  // IdenticalNumbering::IntersectionNumbering
  // -----------------------------------------

  struct IdenticalNumbering::IntersectionNumbering
  {
    enum Side { Inside, Outside };

    template< Side side >
    unsigned int map ( const unsigned int face ) const
    {
      return face;
    }
  };



  // Implementation of IdenticalNumbering
  // ------------------------------------

  template< int codim, int dim, class G, template< int, int, class > class I >
  inline IdenticalNumbering::EntityNumbering< codim >
  IdenticalNumbering::operator[] ( const Entity< codim, dim, G, I > &entity ) const
  {
    return EntityNumbering< codim >();
  }


  template< class G, template< class > class I >
  inline IdenticalNumbering::IntersectionNumbering
  IdenticalNumbering::operator[] ( const Intersection< G, I > &intersection ) const
  {
    return IntersectionNumbering();
  }

}

#endif // #ifndef DUNE_GEOGRID_NUMBERING_HH
