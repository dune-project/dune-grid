// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_REFERENCEELEMENT_HH
#define DUNE_GENERICGEOMETRY_REFERENCEELEMENT_HH

#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/geometrytypes.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Geometry >
    class ReferenceElement;


    template<>
    class ReferenceElement< Point >
    {
      typedef Point Geometry;

    public:
      enum { numCorners = Geometry :: numCorners };
      enum { dimension = Geometry :: dimension };

      template< class ctype >
      inline static void corner ( unsigned int i,
                                  FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        corner_( i, x );
      }

      template< class ctype >
      inline static ctype volume ()
      {
        return ctype( 1 );
      }

    protected:
      template< class ctype, int dim >
      inline static void corner_ ( unsigned int i,
                                   FieldVector< ctype, dim > &x )
      {}
    };


    template< class BaseGeometry >
    class ReferenceElement< Prism< BaseGeometry > >
    {
      typedef Prism< BaseGeometry > Geometry;

    public:
      enum { numCorners = Geometry :: numCorners };
      enum { dimension = Geometry :: dimension };

      template< class ctype >
      inline static void corner ( unsigned int i,
                                  FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        corner_( i, x );
      }

      template< class ctype >
      inline static ctype volume ()
      {
        return ReferenceElement< BaseGeometry > :: volume();
      }

    protected:
      template< class ctype, int dim >
      inline static void corner_ ( unsigned int i,
                                   FieldVector< ctype, dim > &x )
      {
        const unsigned int j = i % ReferenceElement< BaseGeometry > :: numCorners;
        ReferenceElement< BaseGeometry > :: corner_( j, x );
        if( j >= ReferenceElement< BaseGeometry > :: numCorners )
          x[ dimension - 1 ] = ctype( 1 );
      }
    };


    template< class BaseGeometry >
    class ReferenceElement< Pyramid< BaseGeometry > >
    {
      typedef Prism< BaseGeometry > Geometry;

    public:
      enum { numCorners = Geometry :: numCorners };
      enum { dimension = Geometry :: dimension };

      template< class ctype >
      inline static void corner ( unsigned int i,
                                  FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        corner_( i, x );
      }

      template< class ctype >
      inline static ctype volume ()
      {
        const ctype baseVol = ReferenceElement< BaseGeometry > :: volume();

        return baseVol / ctype( Faculty< dimension > :: value );
      }

    protected:
      template< class ctype, int dim >
      inline static void corner_ ( unsigned int i,
                                   FieldVector< ctype, dim > &x )
      {
        if( i < ReferenceElement< BaseGeometry > :: numCorners )
          ReferenceElement< BaseGeometry > :: corner( i, x );
        else
          x[ dimension - 1 ] = ctype( 1 );
      }
    };

  }

}

#endif
