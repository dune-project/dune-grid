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
    class ReferenceElementBase
    {
    public:
      enum { dimension = Geometry :: dimension };

      inline static GeometryType geometryType ()
      {
        return GeometryType( Geometry :: duneType, dimension );
      }
    };



    template< class Geometry >
    class ReferenceElement;


    template< GeometryType :: BasicType basicType >
    class ReferenceElement< Point< basicType > >
      : public ReferenceElementBase< Point< basicType > >
    {
      typedef ReferenceElementBase< Point< basicType > > Base;

    public:
      enum { numCorners = 1 };
      enum { dimension = Base :: dimension };

      template< int codim >
      struct Codim
      {
        enum { numSubEntities = (codim == 0) ? 1 : 0 };
      };

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
      : public ReferenceElementBase< Prism< BaseGeometry > >
    {
      typedef ReferenceElementBase< Prism< BaseGeometry > > Base;

    public:
      enum { numCorners = 2 * ReferenceElement< BaseGeometry > :: numCorners };
      enum { dimension = Base :: dimension };

      template< int codim >
      struct Codim
      {
        enum { numSubEntities = 2 * ReferenceElement< BaseGeometry > :: template Codim< codim - 1 > :: numSubEntities
                                + ReferenceElement< BaseGeometry > :: template Codim< codim > :: numSubEntities };
      };

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
      : public ReferenceElementBase< Pyramid< BaseGeometry > >
    {
      typedef ReferenceElementBase< Pyramid< BaseGeometry > > Base;

    public:
      enum { numCorners = ReferenceElement< BaseGeometry > :: numCorners + 1 };
      enum { dimension = Base :: dimension };

      template< int codim >
      struct Codim
      {
        enum { numSubEntities = ReferenceElement< BaseGeometry > :: template Codim< codim - 1 > :: numSubEntities
                                + ReferenceElement< BaseGeometry > :: template Codim< codim > :: numSubEntities };
      };

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
