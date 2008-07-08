// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_REFERENCEELEMENT_HH
#define DUNE_GENERICGEOMETRY_REFERENCEELEMENT_HH

#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Geometry >
    class ReferenceElement;



    // IntegrationOuterNormal
    // ----------------------

    template< class Geometry >
    class IntegrationOuterNormal;

    template< class BaseGeometry >
    class IntegrationOuterNormal< Prism< BaseGeometry > >
    {
      typedef Prism< BaseGeometry > Geometry;

      enum { dimension = Geometry :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_( unsigned int i,
                             FieldVector< ctype, dim > &n )
      {
        if( i >= NumSubEntities< BaseGeometry, 1 > :: value )
        {
          const unsigned int j = i - NumSubEntities< BaseGeometry, 1 > :: value;
          n[ dimension - 1 ] = (j == 0 ? ctype( -1 ) : ctype( 1 ));
        }
        else
          IntegrationOuterNormal< BaseGeometry > :: evaluate_( i, n );
      }

    public:
      template< class ctype >
      static void evaluate ( unsigned int i,
                             FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };

    template<>
    class IntegrationOuterNormal< Prism< Point > >
    {
      typedef Prism< Point > Geometry;

      enum { dimension = Geometry :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_( unsigned int i,
                             FieldVector< ctype, dim > &n )
      {
        n[ dimension - 1 ] = (i == 0 ? ctype( -1 ) : ctype( 1 ));
      }

    public:
      template< class ctype >
      static void evaluate ( unsigned int i,
                             FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };

    template< class BaseGeometry >
    class IntegrationOuterNormal< Pyramid< BaseGeometry > >
    {
      typedef Pyramid< BaseGeometry > Geometry;

      enum { dimension = Geometry :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_( unsigned int i,
                             FieldVector< ctype, dim > &n )
      {
        typedef SubEntityNumbering< BaseGeometry > Numbering;

        if( i < NumSubEntities< BaseGeometry, 1 > :: value )
        {
          const unsigned int j
            = Numbering :: template subEntity< 1, dimension-2 >( i, 0 );
          FieldVector< ctype, dim > x( ctype( 0 ) );
          ReferenceElement< BaseGeometry > :: corner_( j, x );

          IntegrationOuterNormal< BaseGeometry > :: evaluate_( i, n );
          n[ dimension - 1 ] = (x * n);
        }
        else
          n[ dimension - 1 ] = ctype( -1 );
      }

    public:
      template< class ctype >
      static void evaluate ( unsigned int i,
                             FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };

    template<>
    class IntegrationOuterNormal< Pyramid< Point > >
    {
      typedef Pyramid< Point > Geometry;

      enum { dimension = Geometry :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_( unsigned int i,
                             FieldVector< ctype, dim > &n )
      {
        n[ dimension - 1 ] = (i == 0 ? ctype( -1 ) : ctype( 1 ));
      }

    public:
      template< class ctype >
      static void evaluate ( unsigned int i,
                             FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };



    // ReferenceElement
    // ----------------

    template<>
    class ReferenceElement< Point >
    {
      typedef Point Geometry;

      template< class > friend class ReferenceElement;

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

    private:
      template< class ctype, int dim >
      inline static void corner_ ( unsigned int i,
                                   FieldVector< ctype, dim > &n )
      {}
    };


    template< class BaseGeometry >
    class ReferenceElement< Prism< BaseGeometry > >
    {
      typedef Prism< BaseGeometry > Geometry;

      template< class > friend class ReferenceElement;

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
      inline static void
      integrationOuterNormal ( unsigned int i,
                               FieldVector< ctype, dimension > &n )
      {
        IntegrationOuterNormal< Geometry > :: evaluate( i, n );
      }

      template< class ctype >
      inline static ctype volume ()
      {
        return ReferenceElement< BaseGeometry > :: volume();
      }

    private:
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

      template< class > friend class ReferenceElement;

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
      inline static void
      integrationOuterNormal ( unsigned int i,
                               FieldVector< ctype, dimension > &n )
      {
        IntegrationOuterNormal< Geometry > :: evaluate( i, n );
      }

      template< class ctype >
      inline static ctype volume ()
      {
        const ctype baseVol = ReferenceElement< BaseGeometry > :: volume();

        return baseVol / ctype( Faculty< dimension > :: value );
      }

    private:
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
