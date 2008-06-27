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

    /* struct CoordVectorInterface {
     *   typedef ... field_type;           // e.g. double
     *   typedef ... vector_type;          // e.g. FieldVector
     *   typedef ... derivative_type;      // e.g. FieldMatrix (like vector_type*)
     *   enum {affine = 0|1};              // mapping is of the form Ax+b (used untested)
     *   enum {dim};                       // size of vector_type
     * };
     */
    template< class Geometry, class CoordVector >
    class Mapping;
    template< class Mapping, int face >
    class Normal;

    template<class CoordVector>
    class Mapping < Point, CoordVector >
    {
      typedef Point Geometry;
      enum {dimG = Geometry::dimension};
      enum {dimW = CoordVector::dim};
    public:
      typedef typename CoordVector::field_type FieldType;
      typedef typename CoordVector::vector_type GlobalCoordType;
      typedef typename CoordVector::derivative_type GlobalDerivativeType;
    private:
      GlobalCoordType p_;
      bool zero_;
      void setProperties() {
        zero_ = (p_.two_norm2()<1e-12);
      }
    public:
      explicit Mapping(const CoordVector& coords,int offset=0) :
        p_(coords[offset]),
        zero_(false)
      {
        setProperties();
      }
      // returns Phi : G -> D, Phi_j(x)
      template <class LocalCoordType>
      void phi_set(const LocalCoordType&, GlobalCoordType& p) {
        p = p_;
      }
      template <class LocalCoordType>
      void phi_add(const LocalCoordType&, const FieldType& fac,
                   GlobalCoordType& p) {
        p.axpy(fac,p_);
      }
      // returns (d[i])_j = (d/dx_i Phi_j)
      // e.g. this gives the transpose of the jacobian
      template <class LocalCoordType>
      void deriv_set(const LocalCoordType&, GlobalDerivativeType& d) {
        d = 0;
      }
      template <class LocalCoordType>
      void deriv_add(const LocalCoordType&, const FieldType& fac,
                     GlobalDerivativeType& d) {}
      template <class LocalCoordType>
      double integrationElement(const LocalCoordType&) {
        return 0;
      }
      Mapping& operator-=(const Mapping& other) {
        p_ -= other.p_;
        setProperties();
        return *this;
      }
      Mapping& operator+=(const Mapping& other) {
        p_ += other.p_;
        setProperties();
        return *this;
      }
      // returns if mapping is affine
      bool affine() {
        return true;
      }
      // returns if mapping is constant
      bool constant() {
        return true;
      }
      // returns if mapping is the zero mapping
      bool zero() {
        return zero_;
      }
    };


    template< class BaseGeometry, class CoordVector >
    class Mapping < Prism< BaseGeometry >, CoordVector >
    {
      typedef Prism< BaseGeometry > Geometry;
      enum {dimG = Geometry::dimension};
      enum {dimW = CoordVector::dim};
    public:
      typedef typename CoordVector::field_type FieldType;
      typedef typename CoordVector::vector_type GlobalCoordType;
      typedef typename CoordVector::derivative_type GlobalDerivativeType;
    private:
      Mapping<BaseGeometry,CoordVector> bottom_;
      Mapping<BaseGeometry,CoordVector> top_;
      bool affine_,constant_,zero_;
      void setProperties() {
        affine_ = CoordVector::affine==1 ||
                  (top_.constant() && bottom_.affine() );
        constant_ = (top_.zero() && bottom_.constant());
        zero_ = (top_.zero() && bottom_.zero());
      }
    public:
      explicit Mapping(const CoordVector& coords,int offset=0) :
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
      template <class LocalCoordType>
      void phi_set(const LocalCoordType& x, GlobalCoordType& p) {
        bottom_.phi_set(x,p);
        top_.phi_add(x,x[dimG-1],p);
      }
      template <class LocalCoordType>
      void phi_add(const LocalCoordType& x, const FieldType& fac,
                   GlobalCoordType& p) {
        bottom_.phi_add(x,fac,p);
        top_.phi_add(x,fac*x[dimG-1],p);
      }
      // d[i]_j = db[i]_j + dt[i]_j * x_n (i=1,..,n-1, j=1,..,n)
      // d[n]_j = t(x)
      template <class LocalCoordType>
      void deriv_set(const LocalCoordType& x, GlobalDerivativeType& d) {
        bottom_.deriv_set(x,d);
        top_.deriv_add(x,x[dimG-1],d);
        top_.phi_set(x,d[dimG-1]);
      }
      template <class LocalCoordType>
      void deriv_add(const LocalCoordType& x, const FieldType& fac,
                     GlobalDerivativeType& d) {
        bottom_.deriv_add(x,fac,d);
        top_.deriv_add(x,fac*x[dimG-1],d);
        top_.phi_add(x,fac,d[dimG-1]);
      }
      template <class LocalCoordType>
      double integrationElement(const LocalCoordType& x) {
        GlobalDerivativeType d;
        deriv_set(x,d);
        if (dimW == dimG)
          return std::abs(d.det());
        else {
          return 0.;
          // return sqrt(IntegrationElement::compute(d));
        }
      }
      Mapping& operator-=(const Mapping& other) {
        top_ -= other.top_;
        bottom_ -= other.bottom_;
        setProperties();
        return *this;
      }
      Mapping& operator+=(const Mapping& other) {
        top_ += other.top_;
        bottom_ += other.bottom_;
        setProperties();
        return *this;
      }
      bool affine() {
        return affine_;
      }
      bool constant() {
        return constant_;
      }
      bool zero() {
        return zero_;
      }
    };

    template< class BaseGeometry, class CoordVector >
    class Mapping < Pyramid< BaseGeometry >, CoordVector >
    {
      typedef Prism< BaseGeometry > Geometry;
      enum {dimG = Geometry::dimension};
      enum {dimW = CoordVector::dim};
    public:
      typedef typename CoordVector::field_type FieldType;
      typedef typename CoordVector::vector_type GlobalCoordType;
      typedef typename CoordVector::derivative_type GlobalDerivativeType;
    private:
      Mapping<BaseGeometry,CoordVector> bottom_;
      GlobalCoordType top_,pn_;
      bool affine_,constant_,zero_;
      void setProperties() {
        affine_ = CoordVector::affine==1 ||
                  ( bottom_.affine() );
        constant_ = (top_.zero() && bottom_.constant());
        zero_ = (top_.zero() && bottom_.zero());
      }
    public:
      explicit Mapping(const CoordVector& coords,int offset=0) :
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
      template <class LocalCoordType>
      void phi_set(const LocalCoordType& x, GlobalCoordType& p) {
        if (CoordVector::affine) {       // compile time
          bottom_.phi_set(x,p);
          p.axpy(x[dimG-1],top_);
        }
        else if (affine()) {             // runtime
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
      template <class LocalCoordType>
      void phi_add(const LocalCoordType& x, const FieldType& fac,
                   GlobalCoordType& p) {
        if (CoordVector::affine) {
          bottom_.phi_add(x,fac,p);
          p.axpy(fac*x[dimG-1],top_);
        }
        else if (affine()) {
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
      template <class LocalCoordType>
      void deriv_set(const LocalCoordType& x, GlobalDerivativeType& d) {
        if (CoordVector::affine) {
          bottom_.deriv_set(x,d);
          d[dimG-1] = top_;
        }
        else if (affine()) {
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
      template <class LocalCoordType>
      void deriv_add(const LocalCoordType& x, const FieldType& fac,
                     GlobalDerivativeType& d) {
        // to be implemented
      }
      template <class LocalCoordType>
      double integrationElement(const LocalCoordType& x) {
        GlobalDerivativeType d;
        deriv_set(x,d);
        if (dimW == dimG)
          return std::abs(d.det());
        else {
          return 0.;
          // return sqrt(IntegrationElement::compute(d));
        }
      }
      Mapping& operator-=(const Mapping& other) {
        top_ -= other.top_;
        bottom_ -= other.bottom_;
        setProperties();
        return *this;
      }
      Mapping& operator+=(const Mapping& other) {
        top_ += other.top_;
        bottom_ += other.bottom_;
        setProperties();
        return *this;
      }
      bool affine() {
        return affine_;
      }
      bool constant() {
        return constant_;
      }
      bool zero() {
        return zero_;
      }
    };


  }
}
#endif
