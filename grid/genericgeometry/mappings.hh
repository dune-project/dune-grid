// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAPPINGS_HH
#define DUNE_GENERICGEOMETRY_MAPPINGS_HH

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/geometrytypes.hh>

namespace Dune
{

  namespace GenericGeometry
  {


    /* Mapping have two template arguments:
     * Geometry:    the generic geometry describing
     *              the domain.
     * CoordTraits: a traits class describing the
     *              vector and derivative types needed
     *              to describe the range of the mapping.
     *              This class fixes the local (domain) vector type
     *              and the range (world) vector type. Note
     *              that the dimension of the local coordinates (dimG)
     *              must be greater or equal to Geometry::dimension;
     *              if dimG > Geometry::dimension then only the
     *              first components are used.
     *              The dimension of the global coordinates (dimW)
     *              must be greater or equal to dimG.
     *              Matrix types both for the jacobian (dimW x dimG)
     *              its transpose (dimG x dimW) and
     *              a square matrix of dimension (dimG x dimG).
     * struct CoordTraits {
     *   enum {dimW = };              // world dimension
     *   enum {dimG = };              // grid dimension
     *   typedef ... field_type;
     *   // Vector type of dimension dimG
     *   typedef ... local_type;
     *   // Vector type of dimension dimW
     *   typedef ... global_type;
     *   // Matrix type of dimension = w x g (for jacobian, w>g)
     *   typedef ... derivative_type;
     *   // Matrix type of dimension = g x g
     *   typedef ... localmapping_type;
     *   // Matrix type of dimension = g x w (for jacobian transposed)
     *   typedef ... derivativeT_type;
     *
     *   // Vector of global vectors denoting the edges of the range
     *   // domain, used to construct a mapping together with an offset.
     *   // Corners used are
     *   // p[offset],...,p[offset+Geometry::numCorners]
     *   typedef ... coord_vector;
     *
     *   // mapping is of the form Ax+b (used untested)
     *   enum {affine = 0|1};
     * };
     */

    // Main mapping class.
    // Has the same method as the GenericMapping class
    // defined below which it just forwards.
    // Some additional methods add also implemented like
    // integrationElement, JacobianTransposeInverse
    // It seems not to be possible to make this part
    // of the computation more efficient through
    // the Prism/Pyramid construction !?
    template< class Geometry, class CoordTraits >
    class Mapping;
    // compute the normal to the codim 1 subentity
    // with number face.
    // Use the ansatz:
    //
    template< class Mapping, int face >
    struct Normal;

    // *******************************************
    // *******************************************
    // *******************************************

    template< class Geometry, class CoordTraits >
    struct GenericMapping;

    template<class CoordTraits>
    class GenericMapping < Point, CoordTraits >
    {
      typedef Point Geometry;
    public:
      typedef CoordTraits Traits;
      enum {dimW = CoordTraits :: dimW};
      enum {dimG = CoordTraits :: dimG};
      typedef typename CoordTraits :: field_type FieldType;
      typedef typename CoordTraits :: local_type LocalCoordType;
      typedef typename CoordTraits :: global_type GlobalCoordType;
      typedef typename CoordTraits :: derivative_type JacobianType;
      typedef typename CoordTraits :: localmapping_type SquareMappingType;
      typedef typename CoordTraits :: derivativeT_type JacobianTransposeType;
      typedef typename CoordTraits :: coord_vector CoordVector;
    private:
      GlobalCoordType p_;
      bool zero_;
      void setProperties() {
        zero_ = (p_.two_norm2()<1e-12);
      }
    public:
      explicit GenericMapping(const CoordVector& coords,
                              int offset) :
        p_(coords[offset]),
        zero_(false)
      {
        setProperties();
      }
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
                     JacobianTransposeType& d) const {
        d = 0;
      }
      void deriv_add(const LocalCoordType&,
                     const FieldType& fac,
                     JacobianTransposeType& d) const {}

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

      GenericMapping& operator-=(const GenericMapping& other) {
        p_ -= other.p_;
        setProperties();
        return *this;
      }
      GenericMapping& operator+=(const GenericMapping& other) {
        p_ += other.p_;
        setProperties();
        return *this;
      }
    };


    template< class BaseGeometry, class CoordTraits >
    class GenericMapping < Prism< BaseGeometry >, CoordTraits >
    {
      typedef Prism< BaseGeometry > Geometry;
    public:
      typedef CoordTraits Traits;
      enum {dimW = CoordTraits :: dimW};
      enum {dimG = CoordTraits :: dimG};
      typedef typename CoordTraits :: field_type FieldType;
      typedef typename CoordTraits :: local_type LocalCoordType;
      typedef typename CoordTraits :: global_type GlobalCoordType;
      typedef typename CoordTraits :: derivative_type JacobianType;
      typedef typename CoordTraits :: localmapping_type SquareMappingType;
      typedef typename CoordTraits :: derivativeT_type JacobianTransposeType;
      typedef typename CoordTraits :: coord_vector CoordVector;
    private:
      GenericMapping<BaseGeometry,CoordTraits> bottom_;
      GenericMapping<BaseGeometry,CoordTraits> top_;
      bool affine_,constant_,zero_;
      void setProperties() {
        affine_ = CoordVector::affine==1 ||
                  (top_.constant() && bottom_.affine() );
        constant_ = (top_.zero() && bottom_.constant());
        zero_ = (top_.zero() && bottom_.zero());
      }
    public:
      explicit GenericMapping(const CoordVector& coords,
                              int offset) :
        bottom_(coords,offset),
        top_(coords,offset+BaseGeometry::numCorners),
        affine_(CoordVector::affine==1),
        constant_(false),
        zero_(false)
      {
        top_ -= bottom_;
        setProperties();
      }
      // p = b(x)+t(x)*x_n
      void phi_set(const LocalCoordType& x,
                   GlobalCoordType& p) const {
        bottom_.phi_set(x,p);
        top_.phi_add(x,x[dimG-1],p);
      }
      void phi_add(const LocalCoordType& x,
                   const FieldType& fac,
                   GlobalCoordType& p) const {
        bottom_.phi_add(x,fac,p);
        top_.phi_add(x,fac*x[dimG-1],p);
      }
      // d[i]_j = db[i]_j + dt[i]_j * x_n (i=1,..,n-1, j=1,..,n)
      // d[n]_j = t(x)
      void deriv_set(const LocalCoordType& x,
                     JacobianTransposeType& d) const {
        bottom_.deriv_set(x,d);
        top_.deriv_add(x,x[dimG-1],d);
        top_.phi_set(x,d[dimG-1]);
      }
      void deriv_add(const LocalCoordType& x,
                     const FieldType& fac,
                     JacobianTransposeType& d) const {
        bottom_.deriv_add(x,fac,d);
        top_.deriv_add(x,fac*x[dimG-1],d);
        top_.phi_add(x,fac,d[dimG-1]);
      }

      bool affine() const {
        return CoordVector::affine || affine_;
      }
      bool constant() const {
        return constant_;
      }
      bool zero() const {
        return zero_;
      }

      GenericMapping& operator-=(const GenericMapping& other) {
        top_ -= other.top_;
        bottom_ -= other.bottom_;
        setProperties();
        return *this;
      }
      GenericMapping& operator+=(const GenericMapping& other) {
        top_ += other.top_;
        bottom_ += other.bottom_;
        setProperties();
        return *this;
      }
    };

    template< class BaseGeometry, class CoordTraits >
    class GenericMapping < Pyramid< BaseGeometry >, CoordTraits >
    {
      typedef Pyramid< BaseGeometry > Geometry;
    public:
      typedef CoordTraits Traits;
      enum {dimW = CoordTraits :: dimW};
      enum {dimG = CoordTraits :: dimG};
      typedef typename CoordTraits :: field_type FieldType;
      typedef typename CoordTraits :: local_type LocalCoordType;
      typedef typename CoordTraits :: global_type GlobalCoordType;
      typedef typename CoordTraits :: derivative_type JacobianType;
      typedef typename CoordTraits :: localmapping_type SquareMappingType;
      typedef typename CoordTraits :: derivativeT_type JacobianTransposeType;
      typedef typename CoordTraits :: coord_vector CoordVector;
    private:
      GenericMapping<BaseGeometry,CoordTraits> bottom_;
      GlobalCoordType top_,pn_;
      bool affine_,constant_,zero_;
      void setProperties() {
        affine_ = CoordVector::affine==1 ||
                  ( bottom_.affine() );
        constant_ = (top_.zero() && bottom_.constant());
        zero_ = (top_.zero() && bottom_.zero());
      }
    public:
      explicit GenericMapping(const CoordVector& coords,
                              int offset) :
        bottom_(coords,offset),
        pn_(coords[offset+BaseGeometry::numCorners]),
        top_(coords[offset+BaseGeometry::numCorners]),
        affine_(CoordVector::affine==1),
        constant_(false),
        zero_(false)
      {
        top_ -= coords[offset];
        setProperties();
      }
      // if affine:
      //   p = b(x/(1-y))*(1-y) + pn * y
      //     = db * x + b0*(1-y) + pn * y
      //     = db * x + b0 + y * (pn - p0)
      //     = b(x) + t * y
      // else
      //   p = b(x/(1-y))*(1-y) + pn * y
      void phi_set(const LocalCoordType& x,
                   GlobalCoordType& p) const {
        if (affine()) {
          bottom_.phi_set(x,p);
          p.axpy(x[dimG-1],top_);
        }
        else {
          double h=1.-x[dimG-1];
          double hinv = 1./h;
          bottom_.phi_set(x*hinv,p);
          p *= h;
          p.axpy(x[dimG-1],pn_);
        }
      }
      void phi_add(const LocalCoordType& x,
                   const FieldType& fac,
                   GlobalCoordType& p) const {
        if (affine()) {
          bottom_.phi_add(x,fac,p);
          p.axpy(fac*x[dimG-1],top_);
        }
        else {
          double h=1.-x[dimG-1];
          double hinv = 1./h;
          bottom_.phi_add(x*hinv,fac,p);
          p *= h;
          p.axpy(fac*x[dimG-1],pn_);
        }
      }
      // if affine:
      //   d[i]_j = db[i]_j (i=1,..,n-1, j=1,..,n)
      //   d[n]_j = t_j
      // else (h=1-y)
      //   d[i]_j = db[i]_j*h/h (i=1,..,n-1, j=1,..,n)
      //   d[n]_j = x db(x/h)*s-b(x/h) + pn
      void deriv_set(const LocalCoordType& x,
                     JacobianTransposeType& d) const {
        if (affine()) {
          bottom_.deriv_set(x,d);
          d[dimG-1] = top_;
        }
        else
        {
          //   d[i]_j = db[i]_j*s/s (i=1,..,n-1, j=1,..,n)
          //   d[n]_j = x db(x/s)*s-b(x/s) + pn
          double h=1.-x[dimG-1];
          double hinv = 1./h;
          GlobalCoordType X(x);
          X[dimG-1] = 0;
          X *= hinv;

          bottom_.deriv_set(X,d);
          d[dimG-1] = pn_;
          bottom_.phi_add(X,-1.,d[dimG-1]);
          d.umtv(X,d[dimG-1]);
        }
      }
      void deriv_add(const LocalCoordType& x,
                     const FieldType& fac,
                     JacobianTransposeType& d) const {
        if (affine()) {
          bottom_.deriv_add(x,d);
          d[dimG-1] += top_;
        }
        else
        {
          //   d[i]_j = db[i]_j*s/s (i=1,..,n-1, j=1,..,n)
          //   d[n]_j = x db(x/s)*s-b(x/s) + pn
          double h=1.-x[dimG-1];
          double hinv = 1./h;
          GlobalCoordType X(x);
          X[dimG-1] = 0;
          X *= hinv;

          bottom_.deriv_add(X,d);
          d[dimG-1] += pn_;
          bottom_.phi_add(X,-1.,d[dimG-1]);
          d.umtv(X,d[dimG-1]);
        }
      }

      bool affine() const {
        return CoordVector::affine || affine_;
      }
      bool constant() const {
        return constant_;
      }
      bool zero() const {
        return zero_;
      }

      GenericMapping& operator-=(const GenericMapping& other) {
        top_ -= other.top_;
        bottom_ -= other.bottom_;
        setProperties();
        return *this;
      }
      GenericMapping& operator+=(const GenericMapping& other) {
        top_ += other.top_;
        bottom_ += other.bottom_;
        setProperties();
        return *this;
      }
    };

    template< class Geometry, class CoordTraits >
    class Mapping {
      typedef GenericMapping<Geometry,CoordTraits> GenericMappingType;
    public:
      typedef CoordTraits Traits;
      enum {dimW = CoordTraits :: dimW};
      enum {dimG = CoordTraits :: dimG};
      typedef typename CoordTraits :: field_type FieldType;
      typedef typename CoordTraits :: local_type LocalCoordType;
      typedef typename CoordTraits :: global_type GlobalCoordType;
      typedef typename CoordTraits :: derivative_type JacobianType;
      typedef typename CoordTraits :: localmapping_type SquareMappingType;
      typedef typename CoordTraits :: derivativeT_type JacobianTransposeType;
      typedef typename CoordTraits :: coord_vector CoordVector;
    private:
      GenericMappingType& map_;
    public:
      explicit Mapping(const CoordVector& coords) :
        map_(coords,0) {}
      void phi_set(const LocalCoordType& x,
                   GlobalCoordType& p) const {
        map_.phi_set(x,p);
      }
      void phi_add(const LocalCoordType& x,
                   const FieldType& fac,
                   GlobalCoordType& p) const {
        map_.phi_add(x,fac,p);
      }
      void deriv_set(const LocalCoordType& x,
                     JacobianTransposeType& d) const {
        map_.deriv_set(x,d);
      }
      void deriv_add(const LocalCoordType& x,
                     const FieldType& fac,
                     JacobianTransposeType& d) const {
        map_.deriv_add(x,fac,d);
      }
      bool affine() const {
        return map_.affine();
      }
      // additional methods
      FieldType integrationElement(const LocalCoordType& x)
      {
        JacobianTransposeType d;
        deriv_set(x,d);
        if (dimW==dimG)
          return std::abs(d.det());
        else {
          SquareMappingType R;
          JacobianType dInv;
          return QRDecompose<CoordTraits>::compute(d,R,dInv);
        }
      }
      FieldType jacobianInverseTransposed(const LocalCoordType& x,
                                          JacobianType& dInv)
      {
        JacobianTransposeType d;
        deriv_set(x,d);
        if (dimW==dimG) {
          dInv = d.inverse();
          return std::abs(d.det());
        }
        else {
          SquareMappingType R;
          return QRDecompose<CoordTraits>::compute(d,R,dInv);
        }
      }
    };

  }
}
#endif
