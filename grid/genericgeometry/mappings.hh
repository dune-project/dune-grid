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
    template< int DimG, class CoordTraits >
    struct MappingTraits {
      enum {dimG = DimG};
      enum {dimW = CoordTraits :: dimCoord};
      enum {affine = CoordTraits :: affine};
      typedef typename CoordTraits :: FieldType
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
      typedef typename CoordTraits :: CoordVector CoordVector;
    };


    template< class Topology, class Traits >
    struct GenericMapping;

    template<class Traits>
    class GenericMapping < Point, Traits >
    {
      typedef Point Topology;
    public:
      enum {dim  = Topology::dimension};
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposeType JacobianTransposeType;
      typedef typename Traits :: CoordVector CoordVector;
      static bool isZero(const FieldType& a) {
        return std::abs(a)<1e-12;
      }
    private:
      GlobalCoordType p_;
      bool zero_;
      void setProperties() {
        zero_ = isZero(p_.two_norm2());
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
      // check if x/fac is in domain of phi
      bool inDomain(const LocalCoordType& x,FieldType fac) {
        return (isZero(x[0]));
      }

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


    template< class BaseTopology, class Traits >
    class GenericMapping < Prism< BaseTopology >, Traits >
    {
      typedef Prism< BaseTopology > Topology;
    public:
      enum {dim  = Topology::dimension};
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposeType JacobianTransposeType;
      typedef typename Traits :: CoordVector CoordVector;
    private:
      GenericMapping<BaseTopology,Traits> bottom_;
      GenericMapping<BaseTopology,Traits> top_;
      bool affine_,constant_,zero_;
      void setProperties() {
        affine_ = Traits::affine==1 ||
                  (top_.constant() && bottom_.affine() );
        constant_ = (top_.zero() && bottom_.constant());
        zero_ = (top_.zero() && bottom_.zero());
      }
    public:
      explicit GenericMapping(const CoordVector& coords,
                              int offset) :
        bottom_(coords,offset),
        top_(coords,offset+BaseTopology::numCorners),
        affine_(Traits::affine==1),
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
                     JacobianTransposeType& d) const {
        bottom_.deriv_set(x,d);
        top_.deriv_add(x,x[dim-1],d);
        top_.phi_set(x,d[dim-1]);
      }
      void deriv_add(const LocalCoordType& x,
                     const FieldType& fac,
                     JacobianTransposeType& d) const {
        bottom_.deriv_add(x,fac,d);
        top_.deriv_add(x,fac*x[dim-1],d);
        top_.phi_add(x,fac,d[dim-1]);
      }
      // check if x/fac is in domain of phi
      bool inDomain(const LocalCoordType& x,FieldType fac) {
        return (std::abs(x[dim-1])<1e-12 &&
                std::abs(fac-x[dim-1])<1e-12 &&
                bottom_.inDomain(x,fac));
      }

      bool affine() const {
        return Traits::affine || affine_;
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

    template< class BaseTopology, class Traits >
    class GenericMapping < Pyramid< BaseTopology >, Traits >
    {
      typedef Pyramid< BaseTopology > Topology;
    public:
      enum {dim  = Topology::dimension};
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposeType JacobianTransposeType;
      typedef typename Traits :: CoordVector CoordVector;
    private:
      GenericMapping<BaseTopology,Traits> bottom_;
      GlobalCoordType top_,pn_;
      bool affine_,constant_,zero_;
      void setProperties() {
        affine_ = Traits::affine==1 ||
                  ( bottom_.affine() );
        constant_ = ( (top_.two_norm2()<1e-12) && bottom_.constant());
        zero_ = ( (top_.two_norm2()<1e-12) && bottom_.zero());
      }
    public:
      explicit GenericMapping(const CoordVector& coords,
                              int offset) :
        bottom_(coords,offset),
        top_(coords[offset+BaseTopology::numCorners]),
        pn_(coords[offset+BaseTopology::numCorners]),
        affine_(Traits::affine==1),
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
          p.axpy(x[dim-1],top_);
        }
        else {
          double h=1.-x[dim-1];
          double hinv = 1./h;
          LocalCoordType xx(x);
          xx *= hinv;
          bottom_.phi_set(xx,p);
          p *= h;
          p.axpy(x[dim-1],pn_);
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
          p.axpy(fac*x[dim-1],pn_);
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
          d[dim-1] = pn_;
          bottom_.phi_add(X,-1.,d[dim-1]);
          d.umtv(X,d[dim-1]);
        }
      }
      void deriv_add(const LocalCoordType& x,
                     const FieldType& fac,
                     JacobianTransposeType& d) const {
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
          d[dim-1] += pn_;
          bottom_.phi_add(X,-1.,d[dim-1]);
          d.umtv(X,d[dim-1]);
        }
      }
      bool inDomain(const LocalCoordType& x,FieldType fac) {
        return (std::abs(x[dim-1])<1e-12 &&
                std::abs(fac-x[dim-1])<1e-12 &&
                bottom_.inDomain(x,fac*(1.-x[dim-1])));
      }

      bool affine() const {
        return Traits::affine || affine_;
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
      typedef MappingTraits<Topology::dimension,CoordTraits> Traits;
      typedef GenericMapping<Topology,Traits> GenericMappingType;
    public:
      enum {dimG = Traits :: dimG};
      enum {dimW = Traits :: dimW};
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;
      typedef typename Traits :: JacobianTransposeType JacobianTransposeType;
      typedef typename Traits :: CoordVector CoordVector;

      typedef ReferenceElement<Topology,FieldType> ReferenceElementType;

      template< unsigned int codim, unsigned int subcodim >
      struct SubGeometryCoordVector {
        typedef GlobalCoordType Vector;
        int i_,ii_;
        const CoordVector& coord_;
        SubGeometryCoordVector(const CoordVector& coord,
                               int i,int ii) :
          i_(i), ii_(ii), coord_(coord)
        {}
        const Vector& operator[](int k) {
          const int l = ReferenceElementType :: template subNumbering<codim,subcodim>(i_,ii_);
          return coord_[l];
        }
      };
      template< unsigned int codim, unsigned int subcodim >
      struct SubGeometryTraits : public Traits {
        typedef SubGeometryCoordVector<codim,subcodim> CoordVector;
      };
    protected:
      GenericMappingType map_;
      const CoordVector& coords_;
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
      void phiInvert(const GlobalCoordType& p,
                     LocalCoordType& x) {}
      FieldType volume() const {
        const LocalCoordType& bary = ReferenceElementType::template baryCenter<0>(0);
        return ReferenceElementType::volume() * integrationElement(bary);
      }
      FieldType integrationElement(const LocalCoordType& x) const {
        JacobianTransposeType d;
        jacobianT(x,d);
        return MatrixHelper<CoordTraits>::template detAAT<dimG,dimW>(d);
      }
      FieldType jacobianInverseTransposed(const LocalCoordType& x,
                                          JacobianType& dInv) const {
        JacobianTransposeType d;
        jacobianT(x,d);
        return MatrixHelper<CoordTraits>::template rightInvA<dimG,dimW>(d,dInv);
      }
      void normal(int face,const LocalCoordType& x, GlobalCoordType& n) const {
        JacobianType d;
        FieldType det = jacobianInverseTransposed(x,d);
        const LocalCoordType& refNormal =
          ReferenceElementType::integrationOuterNormal(face);
        MatrixHelper<CoordTraits>::template Ax<dimW,dimG>(d,refNormal,n);
        n *= det;
      }
    };
  }
}
#endif
