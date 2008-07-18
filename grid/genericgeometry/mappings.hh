// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAPPINGS_HH
#define DUNE_GENERICGEOMETRY_MAPPINGS_HH

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/matrix.hh>

namespace Dune
{

  namespace GenericGeometry
  {


    /* Mapping have two template arguments:
     * Topology:    the generic geometry describing
     *              the domain.
     * CoordTraits: a traits class describing the
     *              vector and derivative types needed
     *              to describe the range of the mapping.
     *              This class fixes the local (domain) vector type
     *              and the range (world) vector type. Note
     *              that the dimension of the local coordinates (dimG)
     *              must be greater or equal to Topology::dimension;
     *              if dimG > Topology::dimension then only the
     *              first components are used.
     *              The dimension of the global coordinates (dimW)
     *              must be greater or equal to dimG.
     *              Matrix types both for the jacobian (dimW x dimG)
     *              its transpose (dimG x dimW) and
     *              a square matrix of dimension (dimG x dimG).
     * struct CoordTraits {
     *   enum {dimW = };              // world dimension
     *   enum {dimG = };              // grid dimension
     *   typedef ... FieldType;
     *   // general vector and matrix types
     *   template <int dim> struct Vector { typedef ... Type; };
     *   template <int dimR,dimC> struct Matrix { typedef ... Type; };
     *
     *   // Vector of global vectors denoting the edges of the range
     *   // domain, used to construct a mapping together with an offset.
     *   // Corners used are
     *   // p[offset],...,p[offset+Topology::numCorners]
     *   // Assumption: coord_vector[i] is of type const Vector<dimW>&
     *   typedef ... coord_vector;
     *
     *   // mapping is of the form Ax+b (used untested)
     *   enum {affine = 0|1};
     * };
     */



    // MappingTraits
    // -------------

    template< int DimG, class CoordTraits >
    struct MappingTraits
    {
      enum { dimG = DimG };
      enum { dimW = CoordTraits :: dimCoord };
      enum { affine = CoordTraits :: affine };

      typedef typename CoordTraits :: FieldType FieldType;
      typedef typename CoordTraits :: template Vector< dimG > :: Type LocalCoordType;
      typedef typename CoordTraits :: template Vector< dimW > :: Type GlobalCoordType;

      typedef typename CoordTraits :: template Matrix< dimW, dimG > :: Type JacobianType;
      typedef typename CoordTraits :: template Matrix< dimG, dimW > :: Type
      JacobianTransposedType;
    };



    // Main mapping class.
    // Has the same method as the GenericMapping class
    // defined below which it just forwards.
    // Some additional methods add also implemented like
    // integrationElement, JacobianTransposeInverse
    // It seems not to be possible to make this part
    // of the computation more efficient through
    // the Prism/Pyramid construction !?
    template< class Topology, class CoordTraits >
    class Mapping;

    /*
       template< class Topology,class CoordTraits,
              unsigned int codim>
       struct SubGeometryCoordVector {
       typedef MappingTraits<Topology::dimension,CoordTraits> Traits;
       typedef ReferenceElement<Topology,typename Traits :: FieldType> ReferenceElementType;
       typedef typename Traits :: GlobalCoordType Vector;
       int i_;
       const typename Traits :: CoordVector& coord_;
       SubGeometryCoordVector(const typename Traits :: CoordVector& coord,int i) :
        i_(i), coord_(coord)
       {}
       const Vector& operator[](int k) {
        const int l = ReferenceElementType :: template subNumbering<codim,CoordTraits::dimGrid>(i_,k);
        return coord_[l];
       }
       };
       template< class Topology,class CoordTraits>
       struct SubGeometryCoordVector<Topology,CoordTraits,0> {
       typedef MappingTraits<Topology::dimension,CoordTraits> Traits;
       typedef ReferenceElement<Topology,typename Traits :: FieldType> ReferenceElementType;
       typedef typename Traits :: GlobalCoordType Vector;
       int i_;
       const typename Traits :: CoordVector& coord_;
       SubGeometryCoordVector(const typename Traits :: CoordVector& coord,int i) :
        i_(i), coord_(coord)
       {}
       const Vector& operator[](int k) {
        return coord_[k];
       }
       };
     */



    // GenericMapping
    // --------------

    template< class Topology, class Traits, unsigned int offset = 0 >
    class GenericMapping;

    template< class Traits, unsigned int offset >
    class GenericMapping < Point, Traits, offset >
    {
      typedef Point Topology;

#if 0
      template< class, class, unsigned int >
      friend class GenericMapping;
#endif

    public:
      enum { dim = Topology :: dimension };

      enum { dimW = Traits :: dimW };
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      enum { alwaysAffine = true };

      static bool isZero ( const FieldType &a )
      {
        return std::abs(a)<1e-12;
      }

#if 0
    private:
      GlobalCoordType p_;
      bool zero_;

      void setProperties ()
      {
        zero_ = isZero(p_.two_norm2());
      }

    public:
      template< class CoordVector >
      explicit GenericMapping ( const CoordVector &coords )
        : p_( coords[ offset ] ),
          zero_(false)
      {
        setProperties();
      }
#endif

      template< unsigned int numCorners >
      static const GlobalCoordType &
      origin ( const GlobalCoordType *const (&coords)[ numCorners ] )
      {
        dune_static_assert( (offset < numCorners), "Invalid offset." );
        return *(coords[ offset ]);
      }

      template< unsigned int numCorners >
      static void phi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const GlobalCoordType &y = origin( coords );
        for( unsigned int i = 0; i < dimW; ++i )
          p[ i ] = factor * y[ i ];
      }

      template< unsigned int numCorners >
      static void phi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const GlobalCoordType &y = origin( coords );
        for( unsigned int i = 0; i < dimW; ++i )
          p[ i ] += factor * y[ i ];
      }

      template< unsigned int numCorners >
      static bool Dphi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        return true;
      }

      template< unsigned int numCorners >
      static bool Dphi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        return true;
      }

#if 0
      // returns Phi : G -> D, Phi_j(x)
      void phi_set(const LocalCoordType&,
                   GlobalCoordType& p) const {
        p = p_;
      }
      void phi_add(const LocalCoordType&,
                   const FieldType& fac,
                   GlobalCoordType& p) const {
        p.axpy(fac,p_);
      }
      // returns (d[i])_j = (d/dx_i Phi_j)
      // e.g. this gives the transpose of the jacobian
      void deriv_set(const LocalCoordType&,
                     JacobianTransposedType& d) const {
        d = 0;
      }
      void deriv_add(const LocalCoordType&,
                     const FieldType& fac,
                     JacobianTransposedType& d) const {}
      void origin_add(const FieldType& fac,
                      GlobalCoordType& p) const {
        p.axpy(fac,p_);
      }
#endif

      static bool inDomain ( const LocalCoordType &x, FieldType factor )
      {
        return isZero( x[ 0 ] );
      }

#if 0
      // returns if mapping is affine
      bool affine() const {
        return true;
      }
      // returns if mapping is constant
      bool constant() const {
        return true;
      }
      // returns if mapping is the zero mapping
      bool zero() const {
        return zero_;
      }

      template< unsigned int ofs >
      GenericMapping &operator-= ( const GenericMapping< Point, Traits, ofs > &other )
      {
        p_ -= other.p_;
        setProperties();
        return *this;
      }

      GenericMapping &operator-= ( const GenericMapping &other )
      {
        p_ -= other.p_;
        setProperties();
        return *this;
      }

      template< unsigned int ofs >
      GenericMapping &operator+= ( const GenericMapping< Point, Traits, ofs > &other )
      {
        p_ += other.p_;
        setProperties();
        return *this;
      }

      GenericMapping &operator+= ( const GenericMapping &other )
      {
        p_ += other.p_;
        setProperties();
        return *this;
      }
#endif
    };


    template< class BaseTopology, class Traits, unsigned int offset >
    class GenericMapping< Prism< BaseTopology >, Traits, offset >
    {
      typedef Prism< BaseTopology > Topology;

      typedef GenericMapping< BaseTopology, Traits, offset > BottomMapping;
      typedef GenericMapping< BaseTopology, Traits, offset + BaseTopology :: numCorners >
      TopMapping;

#if 0
      template< class, class, unsigned int >
      friend class GenericMapping;
#endif

    public:
      enum { dim  = Topology :: dimension };

      enum { dimW = Traits :: dimW };
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      enum { alwaysAffine = ((dim < 2) || Traits :: affine) };

#if 0
    private:
      BottomMapping bottom_;
      TopMapping top_;

      bool affine_, constant_, zero_;

      void setProperties ()
      {
        affine_ = (alwaysAffine || (top_.constant() && bottom_.affine()));
        constant_ = (top_.zero() && bottom_.constant());
        zero_ = (top_.zero() && bottom_.zero());
      }

    public:
      template< class CoordVector >
      explicit GenericMapping ( const CoordVector &coords )
        : bottom_( coords ),
          top_( coords ),
          affine_( alwaysAffine ),
          constant_( false ),
          zero_( false )
      {
        top_ -= bottom_;
        setProperties();
      }
#endif

      template< unsigned int numCorners >
      static const GlobalCoordType &
      origin ( const GlobalCoordType *const (&coords)[ numCorners ] )
      {
        return BottomMapping :: origin( coords );
      }

      template< unsigned int numCorners >
      static void phi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = FieldType( 1 ) - xn;
        BottomMapping :: phi_set( coords, x, factor * cxn, p );
        TopMapping :: phi_add( coords, x, factor * xn, p );
      }

      template< unsigned int numCorners >
      static void phi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = FieldType( 1 ) - xn;
        BottomMapping :: phi_add( coords, x, factor * cxn, p );
        TopMapping :: phi_add( coords, x, factor * xn, p );
      }

      template< unsigned int numCorners >
      static bool Dphi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        const FieldType xn = x[ dim-1 ];
        bool affine = true;
        if( alwaysAffine )
        {
          const FieldType cxn = FieldType( 1 ) - xn;
          BottomMapping :: Dphi_set( coords, x, factor * cxn, J );
          TopMapping :: Dphi_add( coords, x, factor * xn, J );
        }
        else
        {
          JacobianTransposedType Jtop;
          affine &= BottomMapping :: Dphi_set( coords, x, factor, J );
          affine &= TopMapping :: Dphi_set( coords, x, factor, Jtop );

          FieldType norm = FieldType( 0 );
          for( unsigned int i = 0; i < dim-1; ++i )
          {
            Jtop[ i ] -= J[ i ];
            norm += Jtop[ i ].two_norm2();
            J[ i ].axpy( xn, Jtop[ i ] );
          }
          affine &= (norm < 1e-12);
        }
        BottomMapping :: phi_set( coords, x, -factor, J[ dim-1 ] );
        TopMapping :: phi_add( coords, x, factor, J[ dim-1 ] );
        return affine;
      }

      template< unsigned int numCorners >
      static bool Dphi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        const FieldType xn = x[ dim-1 ];
        bool affine = true;
        if( alwaysAffine )
        {
          const FieldType cxn = FieldType( 1 ) - xn;
          BottomMapping :: Dphi_add( coords, x, factor * cxn, J );
          TopMapping :: Dphi_add( coords, x, factor * xn, J );
        }
        else
        {
          JacobianTransposedType Jbottom, Jtop;
          affine &= BottomMapping :: Dphi_set( coords, x, FieldType( 1 ), Jbottom );
          affine &= TopMapping :: Dphi_set( coords, x, FieldType( 1 ), Jtop );

          FieldType norm = FieldType( 0 );
          for( unsigned int i = 0; i < dim-1; ++i )
          {
            Jtop[ i ] -= Jbottom[ i ];
            norm += Jtop[ i ].two_norm2();
            J[ i ].axpy( factor, Jbottom[ i ] );
            J[ i ].axpy( factor*xn, Jtop[ i ] );
          }
          affine &= (norm < 1e-12);
        }
        BottomMapping :: phi_add( coords, x, -factor, J[ dim-1 ] );
        TopMapping :: phi_add( coords, x, factor, J[ dim-1 ] );
        return affine;
      }

#if 0
      // p = b(x)+t(x)*x_n
      void phi_set(const LocalCoordType& x,
                   GlobalCoordType& p) const {
        bottom_.phi_set(x,p);
        top_.phi_add(x,x[dim-1],p);
      }
      void phi_add(const LocalCoordType& x,
                   const FieldType& fac,
                   GlobalCoordType& p) const {
        bottom_.phi_add(x,fac,p);
        top_.phi_add(x,fac*x[dim-1],p);
      }
      // d[i]_j = db[i]_j + dt[i]_j * x_n (i=1,..,n-1, j=1,..,n)
      // d[n]_j = t(x)
      void deriv_set(const LocalCoordType& x,
                     JacobianTransposedType& d) const {
        bottom_.deriv_set(x,d);
        top_.deriv_add(x,x[dim-1],d);
        top_.phi_set(x,d[dim-1]);
      }
      void deriv_add(const LocalCoordType& x,
                     const FieldType& fac,
                     JacobianTransposedType& d) const {
        bottom_.deriv_add(x,fac,d);
        top_.deriv_add(x,fac*x[dim-1],d);
        top_.phi_add(x,fac,d[dim-1]);
      }
      void origin_add(const FieldType& fac,
                      GlobalCoordType& p) const {
        bottom_.origin_add(fac,p);
      }
#endif

      // check if x/fac is in domain of phi
      static bool inDomain ( const LocalCoordType &x, FieldType factor )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = factor - xn;
        return (xn > -1e-12) && (cxn > -1e-12)
               && BottomMapping :: inDomain( x, factor );
      }

#if 0
      bool affine() const {
        return Traits::affine || affine_;
      }
      bool constant() const {
        return constant_;
      }
      bool zero() const {
        return zero_;
      }

      template< unsigned int ofs >
      GenericMapping &
      operator-= ( const GenericMapping< Prism< BaseTopology >, Traits, ofs > &other )
      {
        top_ -= other.top_;
        bottom_ -= other.bottom_;
        setProperties();
        return *this;
      }

      GenericMapping &operator-= ( const GenericMapping &other )
      {
        top_ -= other.top_;
        bottom_ -= other.bottom_;
        setProperties();
        return *this;
      }

      template< unsigned int ofs >
      GenericMapping &
      operator+= ( const GenericMapping< Prism< BaseTopology >, Traits, ofs > &other )
      {
        top_ += other.top_;
        bottom_ += other.bottom_;
        setProperties();
        return *this;
      }

      GenericMapping &operator+= ( const GenericMapping &other )
      {
        top_ += other.top_;
        bottom_ += other.bottom_;
        setProperties();
        return *this;
      }
#endif
    };

    template< class BaseTopology, class Traits, unsigned int offset >
    class GenericMapping < Pyramid< BaseTopology >, Traits, offset >
    {
      typedef Pyramid< BaseTopology > Topology;

      typedef GenericMapping< BaseTopology, Traits, offset > BottomMapping;
      typedef GenericMapping< Point, Traits, offset + BaseTopology :: numCorners >
      TopMapping;

#if 0
      template< class, class, unsigned int >
      friend class GenericMapping;
#endif

    public:
      enum {dim  = Topology::dimension};

      enum { dimW = Traits :: dimW };
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      enum { alwaysAffine = (BottomMapping :: alwaysAffine || Traits :: affine) };

#if 0
    private:
      BottomMapping bottom_;
      GlobalCoordType top_;
      bool affine_,constant_,zero_;
      void setProperties() {
        affine_ = Traits::affine==1 ||
                  ( bottom_.affine() );
        constant_ = ( (top_.two_norm2()<1e-12) && bottom_.constant());
        zero_ = ( (top_.two_norm2()<1e-12) && bottom_.zero());
      }

    public:
      template< class CoordVector >
      explicit GenericMapping ( const CoordVector &coords )
        : bottom_( coords ),
          top_( coords[ offset + BaseTopology :: numCorners ] ),
          affine_(Traits::affine==1),
          constant_(false),
          zero_(false)
      {
        top_ -= coords[ offset ];
        setProperties();
      }
#endif

      template< unsigned int numCorners >
      static const GlobalCoordType &
      origin ( const GlobalCoordType *const (&coords)[ numCorners ] )
      {
        return BottomMapping :: origin( coords );
      }

      template< unsigned int numCorners >
      static void phi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const FieldType xn = x[ dim-1 ];
        if( alwaysAffine )
        {
          const GlobalCoordType &top = TopMapping :: origin( coords );
          const GlobalCoordType &bottom = BottomMapping :: origin( coords );

          BottomMapping :: phi_set( coords, x, factor, p );
          for( unsigned int i = 0; i < dimW; ++i )
            p[ i ] += (factor * xn) * (top[ i ] - bottom[ i ]);
        }
        else
        {
          TopMapping :: phi_set( coords, x, factor * xn, p );
          const FieldType cxn = FieldType( 1 ) - xn;
          if( cxn > 1e-12 )
          {
            const FieldType icxn = FieldType( 1 ) / cxn;
            LocalCoordType xb;
            for( unsigned int i = 0; i < dim-1; ++i )
              xb[ i ] = icxn * x[ i ];

            BottomMapping :: phi_add( coords, xb, factor * cxn, p );
          }
        }
      }

      template< unsigned int numCorners >
      static void phi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const FieldType xn = x[ dim-1 ];
        if( alwaysAffine )
        {
          const GlobalCoordType &top = TopMapping :: origin( coords );
          const GlobalCoordType &bottom = BottomMapping :: origin( coords );

          BottomMapping :: phi_add( coords, x, factor, p );
          for( unsigned int i = 0; i < dimW; ++i )
            p[ i ] += (factor * xn) * (top[ i ] - bottom[ i ]);
        }
        else
        {
          TopMapping :: phi_add( coords, x, factor * xn, p );
          const FieldType cxn = FieldType( 1 ) - xn;
          if( cxn > 1e-12 )
          {
            const FieldType icxn = FieldType( 1 ) / cxn;
            LocalCoordType xb;
            for( unsigned int i = 0; i < dim-1; ++i )
              xb[ i ] = icxn * x[ i ];

            BottomMapping :: phi_add( coords, xb, factor * cxn, p );
          }
        }
      }

      template< unsigned int numCorners >
      static bool Dphi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        GlobalCoordType &q = J[ dim-1 ];
        bool affine;
        if( alwaysAffine )
        {
          const GlobalCoordType &top = TopMapping :: origin( coords );
          const GlobalCoordType &bottom = BottomMapping :: origin( coords );

          affine = BottomMapping :: Dphi_set( coords, x, factor, J );
          for( unsigned int i = 0; i < dimW; ++i )
            q[ i ] = factor * (top[ i ] - bottom[ i ]);
        }
        else
        {
          const FieldType xn = x[ dim-1 ];
          const FieldType cxn = FieldType( 1 ) - xn;
          const FieldType icxn = FieldType( 1 ) / cxn;
          LocalCoordType xb;
          for( unsigned int i = 0; i < dim-1; ++i )
            xb[ i ] = icxn * x[ i ];
          affine = BottomMapping :: Dphi_set( coords, xb, factor, J );

          TopMapping :: phi_set( coords, x, factor, q );
          BottomMapping :: phi_add( coords, xb, -factor, q );
          xb *= factor;
          for( unsigned int j = 0; j < dim-1; ++j )
          {
            for( unsigned int i = 0; i < dimW; ++i )
              q[ i ] += J[ j ][ i ] * xb[ j ];
          }
        }
        return affine;
      }

      template< unsigned int numCorners >
      static bool Dphi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        GlobalCoordType &q = J[ dim-1 ];
        bool affine;
        if( alwaysAffine )
        {
          const GlobalCoordType &top = TopMapping :: origin( coords );
          const GlobalCoordType &bottom = BottomMapping :: origin( coords );

          affine = BottomMapping :: Dphi_add( coords, x, factor, J );
          for( unsigned int i = 0; i < dimW; ++i )
            q[ i ] = factor * (top[ i ] - bottom[ i ]);
        }
        else
        {
          const FieldType xn = x[ dim-1 ];
          const FieldType cxn = FieldType( 1 ) - xn;
          const FieldType icxn = FieldType( 1 ) / cxn;
          LocalCoordType xb;
          for( unsigned int i = 0; i < dim-1; ++i )
            xb[ i ] = icxn * x[ i ];
          affine = BottomMapping :: Dphi_add( coords, xb, factor, J );

          TopMapping :: phiaddd( coords, x, factor, q );
          BottomMapping :: phi_add( coords, xb, -factor, q );
          xb *= factor;
          for( unsigned int j = 0; j < dim-1; ++j )
          {
            for( unsigned int i = 0; i < dimW; ++i )
              q[ i ] += J[ j ][ i ] * xb[ j ];
          }
        }
        return affine;
      }

#if 0
      // if affine:
      //   p = b(x/(1-y))*(1-y) + pn * y
      //     = db * x + b0*(1-y) + pn * y
      //     = db * x + b0 + y * (pn - b0)
      //     = b(x) + t * y
      // else
      //   p = b(x/(1-y))*(1-y) + pn * y
      //     = b(x/(1-y))*(1-y) + t * y + b0 * y
      void phi_set(const LocalCoordType& x,
                   GlobalCoordType& p) const {
        if( affine() )
        {
          bottom_.phi_set(x,p);
          p.axpy(x[dim-1],top_);
        }
        else
        {
          double h=1.-x[dim-1];
          double hinv = 1./h;
          LocalCoordType xx(x);
          xx *= hinv;
          bottom_.phi_set(xx,p);
          p *= h;
          p.axpy(x[dim-1],top_);
          bottom_.origin_add(x[dim-1],p);
        }
      }
      void phi_add(const LocalCoordType& x,
                   const FieldType& fac,
                   GlobalCoordType& p) const {
        if (affine()) {
          bottom_.phi_add(x,fac,p);
          p.axpy(fac*x[dim-1],top_);
        }
        else {
          double h=1.-x[dim-1];
          double hinv = 1./h;
          LocalCoordType xx(x);
          xx *= hinv;
          bottom_.phi_add(xx,fac,p);
          p *= h;
          p.axpy(fac*x[dim-1],top_);
          bottom_.origin_add(fac*x[dim-1],p);
        }
      }
      // if affine:
      //   d[i]_j = db[i]_j (i=1,..,n-1, j=1,..,n)
      //   d[n]_j = t_j
      // else (h=1-y)
      //   d[i]_j = db[i]_j*h/h (i=1,..,n-1, j=1,..,n)
      //   d[n]_j = x db(x/h)*h-b(x/h) + pn
      //          = x db(x/h)*h-b(x/h) + t + b0
      void deriv_set(const LocalCoordType& x,
                     JacobianTransposedType& d) const {
        if (affine()) {
          bottom_.deriv_set(x,d);
          d[dim-1] = top_;
        }
        else {
          //   d[i]_j = db[i]_j*s/s (i=1,..,n-1, j=1,..,n)
          //   d[n]_j = x db(x/s)*s-b(x/s) + pn
          double h=1.-x[dim-1];
          double hinv = 1./h;
          LocalCoordType X(x);
          X[dim-1] = 0;
          X *= hinv;

          bottom_.deriv_set(X,d);
          d[dim-1] = top_;
          bottom_.origin_add(FieldType(1.),d[dim-1]);
          bottom_.phi_add(X,-1.,d[dim-1]);
          d.umtv(X,d[dim-1]);
        }
      }
      void deriv_add(const LocalCoordType& x,
                     const FieldType& fac,
                     JacobianTransposedType& d) const {
        if (affine()) {
          bottom_.deriv_add(x,d);
          d[dim-1] += top_;
        }
        else {
          //   d[i]_j = db[i]_j*s/s (i=1,..,n-1, j=1,..,n)
          //   d[n]_j = x db(x/s)*s-b(x/s) + pn
          double h=1.-x[dim-1];
          double hinv = 1./h;
          LocalCoordType X(x);
          X[dim-1] = 0;
          X *= hinv;

          bottom_.deriv_add(X,d);
          d[dim-1] += top_;
          bottom_.origin_add(FieldType(1.),d[dim-1]);
          bottom_.phi_add(X,-1.,d[dim-1]);
          d.umtv(X,d[dim-1]);
        }
      }
      void origin_add(const FieldType& fac,
                      GlobalCoordType& p) const {
        bottom_.origin_add(fac,p);
      }
#endif

      static bool inDomain ( const LocalCoordType &x, FieldType factor )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = factor - xn;
        return (xn > -1e-12) && (cxn > -1e-12)
               && BottomMapping :: inDomain( x, factor * cxn );
      }

#if 0
      bool affine() const {
        return Traits::affine || affine_;
      }
      bool constant() const {
        return constant_;
      }
      bool zero() const {
        return zero_;
      }

      template< unsigned int ofs >
      GenericMapping &
      operator-= ( const GenericMapping< Pyramid< BaseTopology >, Traits, ofs > &other )
      {
        top_ -= other.top_;
        bottom_ -= other.bottom_;
        setProperties();
        return *this;
      }

      GenericMapping &operator-= ( const GenericMapping &other )
      {
        top_ -= other.top_;
        bottom_ -= other.bottom_;
        setProperties();
        return *this;
      }

      template< unsigned int ofs >
      GenericMapping &
      operator+= ( const GenericMapping< Pyramid< BaseTopology >, Traits, ofs > &other )
      {
        top_ += other.top_;
        bottom_ += other.bottom_;
        setProperties();
        return *this;
      }

      GenericMapping &operator+= ( const GenericMapping &other )
      {
        top_ += other.top_;
        bottom_ += other.bottom_;
        setProperties();
        return *this;
      }
#endif
    };



    // Mapping
    // -------


    template< class Topology, class CoordTraits >
    class Mapping
    {
      typedef Mapping< Topology, CoordTraits > ThisType;

    public:
      typedef MappingTraits< Topology :: dimension, CoordTraits > Traits;

      enum {dimG = Traits :: dimG};
      enum {dimW = Traits :: dimW};
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      enum { numCorners = Topology :: numCorners };

      typedef GenericGeometry :: GenericMapping< Topology, Traits > GenericMapping;
      typedef ReferenceElement< Topology, FieldType > ReferenceElementType;

    protected:
      const GlobalCoordType *coords_[ numCorners ];
      //GenericMapping map_;

      mutable JacobianTransposedType jT_;
      mutable JacobianType jTInv_;
      mutable FieldType intEl_;
      mutable GlobalCoordType faceNormal_[Size<Topology,1>::value];
      mutable bool jTComputed, jTInvComputed, intElComputed;
      mutable FieldVector<bool,Size<Topology,1>::value> normalComputed;

    public:
      template< class CoordVector >
      explicit Mapping ( const CoordVector &coords )
      //: coords_( coords, 0 ),
      //: map_( coords ),
        : jTComputed( false ),
          jTInvComputed( false ),
          intElComputed( false ),
          normalComputed( false )
      {
        for( int i = 0; i < numCorners; ++i )
          coords_[ i ] = &(coords[ i ]);
      }

      bool affine () const
      {
        if( GenericMapping :: alwaysAffine )
          return true;
        jacobianT( baryCenter() );
        return jTComputed;
        //return map_.affine();
      }

      const GlobalCoordType &operator[] ( int i ) const
      {
        return *(coords_[ i ]);
      }

      GlobalCoordType global( const LocalCoordType &x ) const
      {
        GlobalCoordType p;
        //map_.phi_set(x,p);
        if( jTComputed )
        {
          MatrixHelper< CoordTraits > :: template ATx< dimG, dimW >( jT_, x, p );
          p += (*this)[ 0 ];
        }
        else
          GenericMapping :: phi_set( coords_, x, FieldType( 1 ), p );
        return p;
      }

      // additional methods
      LocalCoordType local ( const GlobalCoordType &p ) const
      {
        LocalCoordType x;
        GlobalCoordType y = p - (*this)[ 0 ];
        if( jTComputed )
          MatrixHelper< CoordTraits > :: template ATx< dimW, dimG >( jTInv_, y, x );
        else if( affine() )
          local_affine( baryCenter(), y, x );
        else
        {
          x = FieldType( baryCenter() );
          LocalCoordType dx;
          do
          { // DF^n dx^n = -F^n, x^{n+1} += dx^n
            y = p - global( x );
            local_affine( x, y, dx );
            x += dx;
          } while( dx.two_norm2() > 1e-12 );
        }
        return x;
      }

      bool checkInside ( const LocalCoordType &x )
      {
        return GenericMapping :: inDomain( x, FieldType( 1 ) );
      }

      const JacobianTransposedType &jacobianT ( const LocalCoordType &x ) const
      {
        if( !jTComputed )
        {
          //map_.deriv_set(x,jT_);
          jTComputed = GenericMapping :: Dphi_set( coords_, x, FieldType( 1 ), jT_ );
          //jTComputed = affine();
        }
        return jT_;
      }

      const JacobianType &jacobianInverseTransposed ( const LocalCoordType &x ) const
      {
        if( !jTInvComputed )
        {
          const JacobianTransposedType &JT = jacobianT( x );
          intEl_ = MatrixHelper< CoordTraits >
                   :: template rightInvA< dimG, dimW >( JT, jTInv_ );
          jTInvComputed = affine();
          intElComputed = affine();
        }
        return jTInv_;
      }

      FieldType integrationElement ( const LocalCoordType &x ) const
      {
        if( !intElComputed )
        {
          const JacobianTransposedType &JT = jacobianT( x );
          intEl_ = MatrixHelper< CoordTraits > :: template detAAT< dimG, dimW >( JT );
          intElComputed = affine();
        }
        return intEl_;
      }

      const GlobalCoordType &normal ( int face, const LocalCoordType &x ) const
      {
        if( !normalComputed[ face ] )
        {
          const JacobianType &JT = jacobianInverseTransposed( x );
          const LocalCoordType &refNormal
            =  ReferenceElementType :: integrationOuterNormal( face );
          MatrixHelper< CoordTraits >
          :: template Ax< dimW, dimG >( JT, refNormal, faceNormal_[ face ] );
          faceNormal_[ face ] *= intEl_;
          normalComputed[ face ] = affine();
        }
        return faceNormal_[ face ];
      }

      FieldType volume () const
      {
        const FieldType refVolume = ReferenceElementType :: volume();
        return refVolume * integrationElement( baryCenter() );
      }

    protected:
      static const LocalCoordType &baryCenter ()
      {
        return ReferenceElementType :: template baryCenter< 0 >( 0 );
      }

      // affine local to global mapping
      void local_affine ( const LocalCoordType &x,
                          const GlobalCoordType &p,
                          LocalCoordType &y ) const
      {
        const JacobianTransposedType &JT = jacobianT( x );
        MatrixHelper< CoordTraits > :: template xTRightInvA< dimG, dimW >( JT, p, y );
      }
    };


    // If not affine only volume is cached (based on intElCompute)
    // otherwise all quantities can be cached using:
    //   geoCompute:    assign if method called using barycenter
    //   geoPreCompute: assign in constructor using barycenter
    //   geoIsComputed: assign in constructor using barycenter using callback
    enum {geoCompute=0,geoPreCompute=1,geoIsComputed=2};

    template <class Traits>
    struct ComputeAll {
      enum {jTCompute = geoCompute,
            jTInvCompute = geoCompute,
            intElCompute = geoCompute,
            normalCompute = geoCompute};
      void jacobianT(typename Traits::JacobianTransposedType& d) const {}
      void integrationElement(typename Traits::FieldType& intEl) const {}
      void jacobianInverseTransposed(typename Traits::JacobianType& dInv) const {}
      void normal(int face, typename Traits::GlobalCoordType& n) const {}
    };



    // CachedMapping
    // -------------

    template< class Topology, class CoordTraits,
        template< class > class Caching = ComputeAll >
    class CachedMapping
      : public Mapping< Topology, CoordTraits >
    {
      typedef Mapping< Topology, CoordTraits > BaseType;
      typedef CachedMapping< Topology, CoordTraits, Caching > ThisType;

    public:
      enum { dimG = BaseType :: dimG };
      enum { dimW = BaseType :: dimW };
      typedef typename BaseType :: FieldType FieldType;
      typedef typename BaseType :: LocalCoordType LocalCoordType;
      typedef typename BaseType :: GlobalCoordType GlobalCoordType;
      typedef typename BaseType :: JacobianType JacobianType;
      typedef typename BaseType :: JacobianTransposedType JacobianTransposedType;

      typedef Caching< typename BaseType :: Traits > CachingType;
      typedef typename BaseType :: ReferenceElementType ReferenceElementType;

    protected:
      using BaseType :: baryCenter;

    public:
      template< class CoordVector >
      explicit CachedMapping ( const CoordVector &coords,
                               const CachingType &cache = CachingType() )
        : BaseType( coords )
      {
        assert(dim==dimG);
        if (affine()) {
          if (int(CachingType::jTCompute)==geoIsComputed) {
            cache.jacobianT(this->jT_);
            this->jTComputed = true;
          } else if (int(CachingType::jTCompute)==geoPreCompute) {
            BaseType :: jacobianT( baryCenter() );
          }
          if (int(CachingType::jTInvCompute)==geoIsComputed) {
            cache.jacobianInverseTransposed(this->jTInv_);
            this->jTInvComputed = true;
          } else if (int(CachingType::jTInvCompute)==geoPreCompute) {
            BaseType :: jacobianInverseTransposed( baryCenter() );
          }
          if (int(CachingType::intElCompute)==geoIsComputed) {
            cache.integrationElement(this->intEl_);
            this->intElComputed = true;
          } else if (int(CachingType::intElCompute)==geoPreCompute) {
            integrationElement( baryCenter() );
          }
        }
      }

      using BaseType :: affine;
      using BaseType :: operator[];
      using BaseType :: global;
      using BaseType :: local;
      using BaseType :: volume;
      using BaseType :: normal

            const JacobianTransposedType &jacobianT ( const LocalCoordType &x ) const
            {
              if (int(CachingType::jTCompute) == geoCompute || !affine()) {
                BaseType::jacobianT(x);
              }
              return this->jT_;
            }

            // additional methods
            FieldType integrationElement ( const LocalCoordType &x ) const
            {
              if (int(CachingType::intElCompute) == geoCompute || !affine()) {
                BaseType::integrationElement(x);
              }
              return this->intEl_;
            }

            const JacobianType &jacobianInverseTransposed ( const LocalCoordType &x ) const
            {
              if (int(CachingType::jTInvCompute) == geoCompute || !affine()) {
                BaseType::jacobianInverseTransposed(x);
              }
              return this->jTInv_;
            }

            /*
               private:
               template< unsigned int codim>
               struct SubGeometryCoordTraits : public CoordTraits {
               typedef SubGeometryCoordVector<codim> CoordVector;
               };
               public:
               template< unsigned int codim,
                      template<class> class SubCaching = ComputeAll>
               struct SubGeometryType {
               typedef SubGeometryCoordTraits<codim> SubCoordTraits;
               typedef Geometry< SubCoordTraits , SubCachingType > GeometryType;
               GeometryType subGeometry(int i,
                           SubCaching subCache = SubCachingType() ) {
                return GeometryType(SubGeometryCoordVector<codim,subcodim>(this->coords_,i,ii),subCache);
               }
               };
             */
    };

  }

}

#endif
