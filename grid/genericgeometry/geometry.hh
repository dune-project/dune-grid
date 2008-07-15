// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GEOMETRY_HH
#define DUNE_GENERICGEOMETRY_GEOMETRY_HH

#include <dune/common/fvector.hh>
#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/matrix.hh>
#include <dune/grid/genericgeometry/mappings.hh>
#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/common/geometry.hh>

namespace Dune
{
  namespace GenericGeometry
  {
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
      void jacobianT(typename Traits::JacobianTransposeType& d) const {}
      void integrationElement(typename Traits::FieldType& intEl) const {}
      void jacobianInverseTransposed(typename Traits::JacobianType& dInv) const {}
      void normal(int face, typename Traits::GlobalCoordType& n) const {}
    };

    template<class Topology,
        class CoordTraits,
        template <class> class Caching = ComputeAll>
    class Geometry : public Mapping <Topology,CoordTraits> {
      typedef Mapping<Topology,CoordTraits> BaseType;
      typedef Geometry<Topology,CoordTraits> ThisType;
    public:
      enum {dimG = BaseType :: dimG};
      enum {dimW = BaseType :: dimW};
      typedef typename BaseType :: FieldType FieldType;
      typedef typename BaseType :: LocalCoordType LocalCoordType;
      typedef typename BaseType :: GlobalCoordType GlobalCoordType;
      typedef typename BaseType :: JacobianType JacobianType;
      typedef typename BaseType :: JacobianTransposeType JacobianTransposeType;
      typedef typename BaseType :: CoordVector CoordVector;

      typedef Caching<typename BaseType::Traits> CachingType;
      typedef typename BaseType::ReferenceElementType ReferenceElementType;

      JacobianTransposeType jT_;
      JacobianType jTInv_;
      FieldType intEl_;
      GlobalCoordType faceNormal_[Size<Topology,1>::value];
      bool jTComputed, jTInvComputed, intElComputed;
      FieldVector<bool,Size<Topology,1>::value> normalComputed;
      const LocalCoordType& bary_;
      void jacobianT_() {
        if (!jTComputed) {
          BaseType::jacobianT(bary_,jT_);
          jTComputed = true;
        }
      }
      void integrationElement_()
      {
        if (!intElComputed) {
          const JacobianTransposeType& d = jacobianT(bary_);
          intEl_ = MatrixHelper<CoordTraits>::template detAAT<dimG,dimW>(d);
          intElComputed = true;
        }
      }
      void jacobianInverseTransposed_() {
        if (!jTInvComputed) {
          const JacobianTransposeType& d = jacobianT(bary_);
          intEl_ = MatrixHelper<CoordTraits>::template
                   rightInvA<dimG,dimW>(d,jTInv_);
          jTInvComputed = true;
          intElComputed = true;
        }
      }
      void normal_(int face) {
        if (!normalComputed[face]) {
          const JacobianType& d = jacobianInverseTransposed(bary_);
          const LocalCoordType& refNormal =
            ReferenceElementType::integrationOuterNormal(face);
          MatrixHelper<CoordTraits>::template Ax<dimW,dimG>(d,refNormal,faceNormal_[face]);
          faceNormal_[face] *= integrationElement(bary_);
          normalComputed[face] = true;
        }
      }
    public:
      explicit Geometry(const CoordVector& coords,
                        const CachingType& cache = CachingType()) :
        BaseType(coords),
        jTComputed(false),
        jTInvComputed(false),
        intElComputed(false),
        normalComputed(false),
        bary_(ReferenceElementType::template baryCenter<0>(0))
      {
        assert(dim==dimG);
        if (affine()) {
          if (int(CachingType::jTCompute)==geoIsComputed) {
            cache.jacobianT(jT_);
            jTComputed = true;
          } else if (int(CachingType::jTCompute)==geoPreCompute) {
            jacobianT_();
          }
          if (int(CachingType::jTInvCompute)==geoIsComputed) {
            cache.jacobianInverseTransposed(jTInv_);
            jTInvComputed = true;
          } else if (int(CachingType::jTInvCompute)==geoPreCompute) {
            jacobianInverseTransposed_();
          }
          if (int(CachingType::intElCompute)==geoIsComputed) {
            cache.integrationElement(intEl_);
            intElComputed = true;
          } else if (int(CachingType::intElCompute)==geoPreCompute) {
            integrationElement_();
          }
        }
      }
      bool affine() const {
        return BaseType::affine();
      }
      GlobalCoordType global(const LocalCoordType& x) const {
        GlobalCoordType p;
        BaseType::phi(x,p);
        return p;
      }
      const JacobianTransposeType& jacobianT(const LocalCoordType& x) {
        if (int(CachingType::jTCompute) == geoCompute || !affine()) {
          jacobianT_();
        }
        return jT_;
      }
      // additional methods
      FieldType integrationElement(const LocalCoordType& x)
      {
        if (int(CachingType::intElCompute) == geoCompute || !affine()) {
          integrationElement_();
        }
        return intEl_;
      }
      const JacobianType& jacobianInverseTransposed(const LocalCoordType& x) {
        if (int(CachingType::jTInvCompute) == geoCompute || !affine()) {
          jacobianInverseTransposed_();
        }
        return jTInv_;
      }
      const GlobalCoordType& normal(int face,const LocalCoordType& x) {
        if (!normalComputed[face] || !affine()) {
          normal_(face);
          normalComputed[face] = true;
        }
        return faceNormal_[face];
      }
      FieldType volume() {
        const LocalCoordType& bary = ReferenceElementType::template baryCenter<0>(0);
        const FieldType& refVol = ReferenceElementType::volume();
        return integrationElement(bary)*refVol;
      }
    private:
      // tut subcodim=0 das was man moechte?
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
      /*
         template< unsigned int codim, unsigned int subcodim>
         struct SubGeometryCoordTraits : public CoordTraits {
         enum {dimG = Traits::dimG - codim - subcodim};
         typedef SubGeometryCoordVector<codim,subcodim> CoordVector;
         };
       */
    public:
      /*
         template< unsigned int codim, unsigned int subcodim,
                class SubCachingType = ComputeAll<SubGeometryCoordTraits<codim,subcodim> >
         struct SubGeometryType {
         typedef SubGeometryCoordTraits<codim,subcodim> CoordTraits;
         typedef SimplexGeometry< CoordTraits , SubCachingType > GeometryType;
         GeometryType subGeometry(int i,int ii,
                                 SubCachingType subCache = SubCachingType() ) {
          return GeometryType(SubGeometryCoordVector<codim,subcodim>(this->coords_,i,ii),subCache);
         }
         };
       */
    };
  }
}
#endif
