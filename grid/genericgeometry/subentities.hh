// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_SUBENTITIES_HH
#define DUNE_GENERICGEOMETRY_SUBENTITIES_HH

#include <vector>

#include <dune/common/static_assert.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/geometrytypes.hh>
#include <dune/grid/genericgeometry/codimtable.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Geometry, unsigned int codim >
    struct NumSubEntities;

    template< class Geometry, unsigned int codim, unsigned int i >
    struct SubGeometry;

    template< class Geometry, unsigned int codim, unsigned int i,
        unsigned int subcodim, unsigned int j >
    class SubEntityNumber;

    template< class Geometry >
    class SubEntityNumbering;



    // NumSubEntities
    // --------------

    template< class Geometry, unsigned int dim, unsigned int codim >
    class NumSubEntitiesImpl;

    template< unsigned int dim, unsigned int codim >
    class NumSubEntitiesImpl< Point, dim, codim >
    {
      typedef Point Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );

    public:
      enum { value = 1 };
    };

    template< class BaseGeometry, unsigned int dim, unsigned int codim >
    class NumSubEntitiesImpl< Prism< BaseGeometry >, dim, codim >
    {
      typedef Prism< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );

      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };
      enum { n = NumSubEntities< BaseGeometry, codim > :: value };

    public:
      enum { value = n + 2*m };
    };

    template< class BaseGeometry, unsigned int dim >
    class NumSubEntitiesImpl< Prism< BaseGeometry >, dim, 0 >
    {
      typedef Prism< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );

    public:
      enum { value = 1 };
    };

    template< class BaseGeometry, unsigned int dim >
    class NumSubEntitiesImpl< Prism< BaseGeometry >, dim, dim >
    {
      typedef Prism< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );

      enum { m = NumSubEntities< BaseGeometry, dim-1 > :: value };

    public:
      enum { value = 2*m };
    };

    template< class BaseGeometry, unsigned int dim, unsigned int codim >
    struct NumSubEntitiesImpl< Pyramid< BaseGeometry >, dim, codim >
    {
      typedef Pyramid< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );

      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };
      enum { n = NumSubEntities< BaseGeometry, codim > :: value };

    public:
      enum { value = m+n };
    };

    template< class BaseGeometry, unsigned int dim >
    class NumSubEntitiesImpl< Pyramid< BaseGeometry >, dim, 0 >
    {
      typedef Pyramid< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );

    public:
      enum { value = 1 };
    };

    template< class BaseGeometry, unsigned int dim >
    class NumSubEntitiesImpl< Pyramid< BaseGeometry >, dim, dim >
    {
      typedef Pyramid< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );

      enum { m = NumSubEntities< BaseGeometry, dim-1 > :: value };

    public:
      enum { value = m+1 };
    };


    template< class Geometry, unsigned int codim >
    struct NumSubEntities
    {
      enum { value = NumSubEntitiesImpl< Geometry, Geometry :: dimension, codim > :: value };
    };



    // SubGeometry
    // -----------

    template< class Geometry, unsigned int dim, unsigned int codim, unsigned int i >
    class SubGeometryImpl;

    template< unsigned int dim, unsigned int codim, unsigned int i >
    class SubGeometryImpl< Point, dim, codim, i >
    {
      typedef Point Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );
      dune_static_assert( (i < NumSubEntities< Geometry, codim > :: value),
                          "Invalid subentity index" );

    public:
      typedef Geometry type;
    };

    template< class BaseGeometry, unsigned int dim, unsigned int codim, unsigned int i >
    class SubGeometryImpl< Prism< BaseGeometry >, dim, codim, i >
    {
      typedef Prism< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );
      dune_static_assert( (i < NumSubEntities< Geometry, codim > :: value),
                          "Invalid subentity index" );

      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };
      enum { n = NumSubEntities< BaseGeometry, codim > :: value };

      enum { s = (i < n+m ? 0 : 1) };

      template< bool >
      struct PrismSub
      {
        typedef Prism< typename SubGeometry< BaseGeometry, codim, i > :: type > type;
      };

      template< bool >
      struct BaseSub
      {
        typedef typename SubGeometry< BaseGeometry, codim-1, i-(n+s*m) > :: type type;
      };

    public:
      typedef typename ProtectedIf< (i < n), PrismSub, BaseSub > :: type type;
    };

    template< class BaseGeometry, unsigned int dim, unsigned int i >
    class SubGeometryImpl< Prism< BaseGeometry >, dim, 0, i >
    {
      typedef Prism< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );
      dune_static_assert( (i < NumSubEntities< Geometry, 0 > :: value),
                          "Invalid subentity index" );
    public:
      typedef Geometry type;
    };

    template< class BaseGeometry, unsigned int dim, unsigned int i >
    class SubGeometryImpl< Prism< BaseGeometry >, dim, dim, i >
    {
      typedef Prism< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );
      dune_static_assert( (i < NumSubEntities< Geometry, dim > :: value),
                          "Invalid subentity index" );
    public:
      typedef Point type;
    };

    template< class BaseGeometry, unsigned int dim, unsigned int codim, unsigned int i >
    class SubGeometryImpl< Pyramid< BaseGeometry >, dim, codim, i >
    {
      typedef Pyramid< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );
      dune_static_assert( (i < NumSubEntities< Geometry, codim > :: value),
                          "Invalid subentity index" );

      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };

      template< bool >
      struct BaseSub
      {
        typedef typename SubGeometry< BaseGeometry, codim-1, i > :: type type;
      };

      template< bool >
      struct PyramidSub
      {
        typedef Pyramid< typename SubGeometry< BaseGeometry, codim, i-m > :: type > type;
      };

    public:
      typedef typename ProtectedIf< (i < m), BaseSub, PyramidSub > :: type type;
    };

    template< class BaseGeometry, unsigned int dim, unsigned int i >
    class SubGeometryImpl< Pyramid< BaseGeometry >, dim, 0, i >
    {
      typedef Pyramid< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );
      dune_static_assert( (i < NumSubEntities< Geometry, 0 > :: value),
                          "Invalid subentity index" );

    public:
      typedef Geometry type;
    };

    template< class BaseGeometry, unsigned int dim, unsigned int i >
    class SubGeometryImpl< Pyramid< BaseGeometry >, dim, dim, i >
    {
      typedef Pyramid< BaseGeometry > Geometry;
      dune_static_assert( (dim == Geometry :: dimension), "Wrong dimension" );
      dune_static_assert( (i < NumSubEntities< Geometry, dim > :: value),
                          "Invalid subentity index" );

    public:
      typedef Point type;
    };

    template< class Geometry, unsigned int codim, unsigned int i >
    struct SubGeometry
    {
      typedef typename SubGeometryImpl< Geometry, Geometry :: dimension, codim, i > :: type type;
    };



    // SubEntityNumber
    // ---------------

    template< class Geometry, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int subcodim, unsigned int j >
    struct SubEntityNumberHelper;

    template< class BaseGeometry, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int subcodim, unsigned int j >
    struct SubEntityNumberHelper< Prism< BaseGeometry >, codim, i, subdim, subcodim, j >
    {
      typedef Prism< BaseGeometry > Geometry;
      typedef typename GenericGeometry :: SubGeometry< Geometry, codim, i > :: type
      SubGeometry;

      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };
      enum { n = NumSubEntities< BaseGeometry, codim > :: value };

      enum { mb = NumSubEntities< BaseGeometry, codim+subcodim-1 > :: value };
      enum { nb = NumSubEntities< BaseGeometry, codim+subcodim > :: value };

      enum { s = (i < n+m ? 0 : 1) };

      template< bool >
      struct PrismSub
      {
        typedef typename GenericGeometry :: BaseGeometry< SubGeometry > :: type
        SubBaseGeometry;

        enum { ms = NumSubEntities< SubBaseGeometry, subcodim-1 > :: value };
        enum { ns = NumSubEntities< SubBaseGeometry, subcodim > :: value };

        enum { s = (j < ns+ms ? 0 : 1) };

        template< bool >
        struct PrismSubSub
        {
          enum { value = SubEntityNumber< BaseGeometry, codim, i, subcodim, j > :: value };
        };

        template< bool >
        struct BaseSubSub
        {
          enum
          {
            value = SubEntityNumber< BaseGeometry, codim, i, subcodim-1, j-(ns+s*ms) > :: value
                    + nb + s*mb
          };
        };

        enum { value = ProtectedIf< (j < ns), PrismSubSub, BaseSubSub > :: value };
      };

      template< bool >
      struct BaseSub
      {
        enum
        {
          value = SubEntityNumber< BaseGeometry, codim-1, i-(n+s*m), subcodim, j > :: value
                  + nb + s*mb
        };
      };

      enum { value = ProtectedIf< (i < n), PrismSub, BaseSub > :: value };
    };

    template< class BaseGeometry, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int j >
    struct SubEntityNumberHelper< Prism< BaseGeometry >, codim, i, subdim, 0, j >
    {
      enum { value = i };
    };

    template< class BaseGeometry, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int j >
    struct SubEntityNumberHelper< Prism< BaseGeometry >, codim, i, subdim, subdim, j >
    {
      typedef Prism< BaseGeometry > Geometry;

      typedef typename GenericGeometry :: SubGeometry< Geometry, codim, i > :: type
      SubGeometry;

      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };
      enum { n = NumSubEntities< BaseGeometry, codim > :: value };

      enum { mb = NumSubEntities< BaseGeometry, codim+subdim-1 > :: value };

      enum { s = (i < n+m ? 0 : 1) };

      template< bool >
      struct PrismSub
      {
        typedef typename GenericGeometry :: BaseGeometry< SubGeometry > :: type
        SubBaseGeometry;

        enum { ms = NumSubEntities< SubBaseGeometry, subdim-1 > :: value };

        enum { s = (j < ms ? 0 : 1) };

        enum
        {
          value = SubEntityNumber< BaseGeometry, codim, i, subdim-1, j-s*ms > :: value
                  + s*mb
        };
      };

      template< bool >
      struct BaseSub
      {
        enum
        {
          value = SubEntityNumber< BaseGeometry, codim-1, i-(n+s*m), subdim, j > :: value
                  + s*mb
        };
      };

      enum { value = ProtectedIf< (i < n), PrismSub, BaseSub > :: value };
    };

    template< class BaseGeometry, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int subcodim, unsigned int j >
    struct SubEntityNumberHelper< Pyramid< BaseGeometry >, codim, i, subdim, subcodim, j >
    {
      typedef Pyramid< BaseGeometry > Geometry;

      typedef typename GenericGeometry :: SubGeometry< Geometry, codim, i > :: type
      SubGeometry;

      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };

      enum { mb = NumSubEntities< BaseGeometry, codim+subcodim-1 > :: value };

      template< bool >
      struct BaseSub
      {
        enum
        {
          value = SubEntityNumber< BaseGeometry, codim-1, i, subcodim, j > :: value
        };
      };

      template< bool >
      struct PyramidSub
      {
        typedef typename GenericGeometry :: BaseGeometry< SubGeometry > :: type
        SubBaseGeometry;

        enum { ms = NumSubEntities< SubBaseGeometry, subcodim-1 > :: value };

        template< bool >
        struct BaseSubSub
        {
          enum { value = SubEntityNumber< BaseGeometry, codim, i-m, subcodim-1, j > :: value };
        };

        template< bool >
        struct PyramidSubSub
        {
          enum
          {
            value = SubEntityNumber< BaseGeometry, codim, i-m, subcodim, j-ms > :: value
                    + mb
          };
        };

        enum { value = ProtectedIf< (j < ms), BaseSubSub, PyramidSubSub > :: value };
      };

      enum { value = ProtectedIf< (i < m), BaseSub, PyramidSub > :: value };
    };

    template< class BaseGeometry, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int j >
    struct SubEntityNumberHelper< Pyramid< BaseGeometry >, codim, i, subdim, 0, j >
    {
      enum { value = i };
    };

    template< class BaseGeometry, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int j >
    struct SubEntityNumberHelper< Pyramid< BaseGeometry >, codim, i, subdim, subdim, j >
    {
      typedef Pyramid< BaseGeometry > Geometry;

      typedef typename GenericGeometry :: SubGeometry< Geometry, codim, i > :: type
      SubGeometry;

      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };

      enum { mb = NumSubEntities< BaseGeometry, codim+subdim-1 > :: value };

      template< bool >
      struct BaseSub
      {
        enum
        {
          value = SubEntityNumber< BaseGeometry, codim-1, i, subdim, j > :: value
        };
      };

      template< bool >
      struct PyramidSub
      {
        typedef typename GenericGeometry :: BaseGeometry< SubGeometry > :: type
        SubBaseGeometry;

        enum { ms = NumSubEntities< SubBaseGeometry, subdim-1 > :: value };

        template< bool >
        struct BaseSubSub
        {
          enum { value = SubEntityNumber< BaseGeometry, codim, i-m, subdim-1, j > :: value };
        };

        template< bool >
        struct PyramidSubSub
        {
          enum { value = mb };
        };

        enum { value = ProtectedIf< (j < ms), BaseSubSub, PyramidSubSub > :: value };
      };

      enum { value = ProtectedIf< (i < m), BaseSub, PyramidSub > :: value };
    };


    template< class Geometry, unsigned int codim, unsigned int i,
        unsigned int subcodim, unsigned int j >
    class SubEntityNumber
    {
      dune_static_assert( (codim <= Geometry :: dimension), "Invald codimension" );
      dune_static_assert( (i < NumSubEntities< Geometry, codim > :: value),
                          "Invalid subentity index" );

      typedef typename GenericGeometry :: SubGeometry< Geometry, codim, i > :: type
      SubGeometry;
      dune_static_assert( (subcodim <= SubGeometry :: dimension), "Invald subcodimension" );
      dune_static_assert( (j < NumSubEntities< SubGeometry, subcodim > :: value),
                          "Invalid subsubentity index" );

      template< bool >
      struct Border
      {
        enum { value = (codim == 0 ? j : i) };
      };

      template< bool >
      struct Interior
      {
        enum
        {
          value = SubEntityNumberHelper< Geometry, codim, i, SubGeometry :: dimension,
              subcodim, j > :: value
        };
      };

    public:
      enum
      {
        value = ProtectedIf< (codim == 0) || (codim == Geometry :: dimension),
            Border, Interior > :: value
      };
    };



    // SubEntityNumbering
    // ------------------

    template< class Geometry >
    class SubEntityNumbering
    {
      template< int codim >
      class Numbering;

      template< int codim >
      class CodimNumbering
        : public std :: vector< Numbering< codim > >
      {};

      template< int codim >
      struct Builder;

      CodimTable< CodimNumbering, Geometry :: dimension > codimNumbering_;

    public:
      template< unsigned int codim, unsigned int subcodim >
      static unsigned int subEntity ( unsigned int i, unsigned int j )
      {
        Int2Type< codim > codimVariable;
        const CodimNumbering< codim > &numbering = instance().codimNumbering_[ codimVariable ];
        return numbering[ i ].template number< subcodim >( j );
      }

    private:
      SubEntityNumbering ()
      {
        ForLoop< Builder, 0, Geometry :: dimension > :: apply( *this );
      }

      static const SubEntityNumbering &instance ()
      {
        static SubEntityNumbering inst;
        return inst;
      }
    };



    template< class Geometry >
    template< int codim >
    class SubEntityNumbering< Geometry > :: Numbering
    {
      std :: vector< unsigned int > numbering_[ Geometry :: dimension - codim + 1 ];

    public:
      template< int subcodim >
      const unsigned int &number ( unsigned int j ) const
      {
        return numbering_[ subcodim ][ j ];
      }

      template< int subcodim >
      unsigned int &number ( unsigned int j )
      {
        return numbering_[ subcodim ][ j ];
      }

      template< int subcodim >
      void resize ( unsigned int size )
      {
        numbering_[ subcodim ].resize( size );
      }
    };



    template< class Geometry >
    template< int codim >
    struct SubEntityNumbering< Geometry > :: Builder
    {
      typedef GenericGeometry :: SubEntityNumbering< Geometry > SubEntityNumbering;
      typedef typename SubEntityNumbering :: template CodimNumbering< codim > CodimNumbering;
      typedef typename SubEntityNumbering :: template Numbering< codim > Numbering;

      template< int i >
      struct Sub
      {
        typedef typename GenericGeometry :: SubGeometry< Geometry, codim, i > :: type
        SubGeometry;

        template< int subcodim >
        struct SubCodim
        {
          template< int j >
          struct SubSub
          {
            static void apply( Numbering &numbering )
            {
              numbering.template number< subcodim >( j )
                = SubEntityNumber< Geometry, codim, i, subcodim, j > :: value;
            }
          };

          static void apply ( Numbering &numbering )
          {
            enum { numSubSubs = NumSubEntities< SubGeometry, subcodim > :: value };
            numbering.template resize< subcodim >( numSubSubs );
            ForLoop< SubSub, 0, numSubSubs-1 > :: apply( numbering );
          }
        };

        static void apply ( CodimNumbering &numbering )
        {
          ForLoop< SubCodim, 0, SubGeometry :: dimension > :: apply( numbering[ i ] );
        }
      };

      static void apply ( SubEntityNumbering &numbering )
      {
        enum { numSubs = NumSubEntities< Geometry, codim > :: value };
        Int2Type< codim > codimVariable;
        CodimNumbering &codimNumbering = numbering.codimNumbering_[ codimVariable ];
        codimNumbering.resize( numSubs );
        ForLoop< Sub, 0, numSubs-1 > :: apply( codimNumbering );
      }
    };

  }

}

#endif
