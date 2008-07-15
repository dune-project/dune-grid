// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMTRY_HYBRIDGEOMETRIES_HH
#define DUNE_GENERICGEOMTRY_HYBRIDGEOMETRIES_HH

namespace Dune
{

  namespace GenericGeometry
  {

    template< unsigned int dimG, class Traits >
    struct HybridGeometry
    {
      enum { dimGrid = dimG };
      enum { dimWorld = Traits :: dimW };

      typedef typename Traits :: FieldType Field;
      typedef typename Traits :: template Vector< dimGrid > :: Type
      LocalCoordinate;
      typedef typename Traits :: template Vector< dimWorld > :: Type
      GlobalCoordinate;
      typedef typename Traits :: template Matrix< dimWorld, dimGrid > :: Type
      Jacobian;

      virtual ~HybridGeometry ()
      {}

      virtual GeometryType type () const = 0;

      virtual int corners () const = 0;

      virtual const GlobalCoordinate &operator[] ( int i ) const = 0;

      virtual GlobalCoordinate global ( const LocalCoordinate &local ) const = 0;

      virtual LocalCoordinate local ( const GlobalCoordinate &global ) const = 0;

      virtual bool checkInside ( const LocalCoordinate &local ) const = 0;

      virtual Field integrationElement ( const LocalCoordinate &local ) const = 0;

      virtual Field volume () const = 0;

      virtual const Jacobian &
      jacobianInverseTransposed ( const LocalCoordinate &local ) const = 0;
    };



    template< class Topology, class Traits >
    class VirtualGeometry
      : public HybridGeometry< Topology :: dimension, Traits >
    {
      Geometry< Topology, Traits > geometry_;

    public:
      explicit VirtualGeometry ( const typename Traits :: coord_vector &coordVector )
        : geometry_( coordVector )
      {}

      virtual GeometryType type () const
      {
        return geometry_.type();
      }

      virtual int corners () const
      {
        return geometry_.corners();
      }

      virtual const GlobalCoordinate &operator[] ( int i ) const
      {
        return geoetry_[ i ];
      }

      virtual GlobalCoordinate global ( const LocalCoordinate &local ) const
      {
        return geometry_.global( local );
      }

      virtual LocalCoordinate local ( const GlobalCoordinate &global ) const
      {
        return geometry_.local( global );
      }

      virtual bool checkInside ( const LocalCoordinate &local ) const
      {
        return geometry_.checkInside( local );
      }

      virtual Field integrationElement ( const LocalCoordinate &local ) const
      {
        return geometry_.integrationElement( local );
      }

      virtual Field volume () const
      {
        return geometry_.volume();
      }

      virtual const JacobianType &
      jacobianInverseTransposed ( const LocalCoordinate &local ) const
      {
        return geometry_.jacobianInverseTransposed( local );
      }
    };

  }

}

#endif
