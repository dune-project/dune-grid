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

      using BaseType :: baryCenter;

    public:
      explicit Geometry(const CoordVector& coords,
                        const CachingType& cache = CachingType()) :
        BaseType(coords)
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
      using BaseType::affine;
      using BaseType::operator[];
      using BaseType::global;
      using BaseType::local;
      using BaseType::volume;
      using BaseType::normal;
      const JacobianTransposeType& jacobianT(const LocalCoordType& x) const {
        if (int(CachingType::jTCompute) == geoCompute || !affine()) {
          BaseType::jacobianT(x);
        }
        return this->jT_;
      }
      // additional methods
      FieldType integrationElement(const LocalCoordType& x) const
      {
        if (int(CachingType::intElCompute) == geoCompute || !affine()) {
          BaseType::integrationElement(x);
        }
        return this->intEl_;
      }
      const JacobianType& jacobianInverseTransposed(const LocalCoordType& x) const {
        if (int(CachingType::jTInvCompute) == geoCompute || !affine()) {
          BaseType::jacobianInverseTransposed(x);
        }
        return this->jTInv_;
      }
    private:
      template< unsigned int codim>
      struct SubGeometryCoordVector {
        typedef GlobalCoordType Vector;
        int i_;
        const CoordVector& coord_;
        SubGeometryCoordVector(const CoordVector& coord,int i) :
          i_(i), coord_(coord)
        {}
        const Vector& operator[](int k) {
          const int l = ReferenceElementType :: template subNumbering<codim,dimG-codim>(i_,k);
          return coord_[l];
        }
      };
      template< unsigned int codim>
      struct SubGeometryCoordTraits : public CoordTraits {
        typedef SubGeometryCoordVector<codim> CoordVector;
      };
    public:
      /*
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
