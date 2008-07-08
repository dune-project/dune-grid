// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAPPINGS_HH
#define DUNE_GENERICGEOMETRY_MAPPINGS_HH

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>

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
     *   typedef ... field_type;
     *   // general vector and matrix types
     *   template <int dim> struct Vector { typedef ... Type; };
     *   template <int dimR,dimC> struct Matrix { typedef ... Type; };
     *
     *   // Vector of global vectors denoting the edges of the range
     *   // domain, used to construct a mapping together with an offset.
     *   // Corners used are
     *   // p[offset],...,p[offset+Topology::numCorners]
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
    template< class Topology, class CoordTraits >
    class Mapping;

    // *******************************************
    // *******************************************
    // *******************************************

    template< class Topology, class CoordTraits >
    struct GenericMapping;

    template<class CoordTraits>
    class GenericMapping < Point, CoordTraits >
    {
      typedef Point Topology;
    public:
      typedef CoordTraits Traits;
      enum {dimG = Topology :: dimension};
      enum {dimW = CoordTraits :: dimW};
      typedef typename CoordTraits :: field_type
      FieldType;
      typedef typename CoordTraits :: template Vector <dimG> :: Type
      LocalCoordType;
      typedef typename CoordTraits :: template Vector<dimW> :: Type
      GlobalCoordType;
      typedef typename CoordTraits :: template Matrix<dimW,dimG> :: Type
      JacobianType;
      typedef typename CoordTraits :: template Matrix<dimG,dimG> :: Type
      SquareMappingType;
      typedef typename CoordTraits :: template Matrix<dimG,dimW> :: Type
      JacobianTransposeType;
      typedef typename CoordTraits :: coord_vector
      CoordVector;
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


    template< class BaseTopology, class CoordTraits >
    class GenericMapping < Prism< BaseTopology >, CoordTraits >
    {
      typedef Prism< BaseTopology > Topology;
    public:
      typedef CoordTraits Traits;
      enum {dimG = Topology :: dimension};
      enum {dimW = CoordTraits :: dimW};
      typedef typename CoordTraits :: field_type
      FieldType;
      typedef typename CoordTraits :: template Vector<dimG> :: Type
      LocalCoordType;
      typedef typename CoordTraits :: template Vector<dimW> :: Type
      GlobalCoordType;
      typedef typename CoordTraits :: template Matrix<dimW,dimG> :: Type
      JacobianType;
      typedef typename CoordTraits :: template Matrix<dimG,dimG> :: Type
      SquareMappingType;
      typedef typename CoordTraits :: template Matrix<dimG,dimW> :: Type
      JacobianTransposeType;
      typedef typename CoordTraits :: coord_vector
      CoordVector;
    private:
      GenericMapping<BaseTopology,CoordTraits> bottom_;
      GenericMapping<BaseTopology,CoordTraits> top_;
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
        top_(coords,offset+BaseTopology::numCorners),
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

    template< class BaseTopology, class CoordTraits >
    class GenericMapping < Pyramid< BaseTopology >, CoordTraits >
    {
      typedef Pyramid< BaseTopology > Topology;
    public:
      typedef CoordTraits Traits;
      enum {dimG = Topology :: dimension};
      enum {dimW = CoordTraits :: dimW};
      typedef typename CoordTraits :: field_type
      FieldType;
      typedef typename CoordTraits :: template Vector<dimG> :: Type
      LocalCoordType;
      typedef typename CoordTraits :: template Vector<dimW> :: Type
      GlobalCoordType;
      typedef typename CoordTraits :: template Matrix<dimW,dimG> :: Type
      JacobianType;
      typedef typename CoordTraits :: template Matrix<dimG,dimG> :: Type
      SquareMappingType;
      typedef typename CoordTraits :: template Matrix<dimG,dimW> :: Type
      JacobianTransposeType;
      typedef typename CoordTraits :: coord_vector
      CoordVector;
    private:
      GenericMapping<BaseTopology,CoordTraits> bottom_;
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
        pn_(coords[offset+BaseTopology::numCorners]),
        top_(coords[offset+BaseTopology::numCorners]),
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

    template< class Topology, class CoordTraits >
    class Mapping {
      typedef Mapping<Topology,CoordTraits> ThisType;
      typedef GenericMapping<Topology,CoordTraits> GenericMappingType;
    public:
      typedef CoordTraits Traits;
      enum {dimG = Topology :: dimension};
      enum {dimW = CoordTraits :: dimW};
      typedef typename CoordTraits :: field_type
      FieldType;
      typedef typename CoordTraits :: template Vector<dimG> :: Type
      LocalCoordType;
      typedef typename CoordTraits :: template Vector<dimW> :: Type
      GlobalCoordType;
      typedef typename CoordTraits :: template Matrix<dimW,dimG> :: Type
      JacobianType;
      typedef typename CoordTraits :: template Matrix<dimG,dimG> :: Type
      SquareMappingType;
      typedef typename CoordTraits :: template Matrix<dimG,dimW> :: Type
      JacobianTransposeType;
      typedef typename CoordTraits :: coord_vector
      CoordVector;
      template< int c, int cc >
      struct SubEntityCoordVector {
        typedef GlobalCoordType Vector;
        int i_,ii_;
        const CoordVector& coord_;
        SubEntityCoordVector(const CoordVector& coord,
                             int i,int ii) :
          i_(i), ii_(ii), coord_(coord)
        {}
        const Vector& operator[](int k) {
          return coord_[k];
        }
      };
      template< int c, int cc >
      struct SubEntityTraits : public Traits {
        typedef SubEntityCoordVector<c,cc> CoordVector;
      };
    private:
      GenericMappingType map_;
      const CoordVector coords_;
      JacobianTransposeType d;
      SquareMappingType L;
      JacobianType Q;
    public:
      explicit Mapping(const CoordVector& coords) :
        map_(coords,0),
        coords_(coords)
      {}
      void phi(const LocalCoordType& x,
               GlobalCoordType& p) const {
        map_.phi_set(x,p);
      }
      void jacobianT(const LocalCoordType& x,
                     JacobianTransposeType& d) const {
        map_.deriv_set(x,d);
      }
      bool affine() const {
        return map_.affine();
      }
      // additional methods
      /*
         FieldType integrationElement(const LocalCoordType& x)
         {
         jacobianT(x,d);
         if (dimW==dimG)
          return std::abs(d.det());
         else {
          return QRDecompose<CoordTraits>::compute(d,L,Q);
         }
         }
         FieldType jacobianInverseTransposed(const LocalCoordType& x,
                                          JacobianType& dInv)
         {
         if (dimW==dimG) {
          deriv(x,dInv);
          dInv.invert();
          return std::abs(dInv.det());
         }
         else {
          double intEl = QRDecompose<CoordTraits>::compute(d,L,Q);
          QRDecompose<CoordTraits>::invert(L,Q,dInv);
          return intEl;
         }
         }
         FieldType volume() {
         return 0.;
         }
         void noraml(int face,GlobalCoordType& n) {
         }
       */
    };
  }
}
#endif
