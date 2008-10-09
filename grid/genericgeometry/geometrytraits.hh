// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GEOMETRYTRAITS_HH
#define DUNE_GENERICGEOMETRY_GEOMETRYTRAITS_HH

#include <dune/common/geometrytype.hh>

#include <dune/grid/genericgeometry/matrix.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class ctype, int dimG, int dimW >
    struct DefaultGeometryTraits
    {
      typedef DuneCoordTraits< ctype > CoordTraits;

      static const int dimGrid = dimG;
      static const int dimWorld = dimW;

      //   hybrid   [ true if Codim 0 is hybrid ]
      static const bool hybrid = true;
      //   dunetype [ for Codim 0, needed for (hybrid=false) ]
      // static const GeometryType :: BasicType dunetype = GeometryType :: simplex;

      // what basic geometry type shall the line be considered?
      static const GeometryType :: BasicType linetype = GeometryType :: simplex;

      template< class Topology >
      struct Mapping
      {
        typedef MappingTraits< CoordTraits, Topology :: dimension, dimWorld > Traits;
        typedef CoordPointerStorage< Topology, typename Traits :: GlobalCoordType >
        CornerStorage;
        typedef CornerMapping< Topology, Traits, CornerStorage > type;
      };

      template< class Traits >
      struct Caching
        : public ComputeAll< Traits >
      {};
    };



    /** \struct  GlobalGeometryTraits
     *  \ingroup GenericGeometry
     *  \brief   grid specific information required by GenericGeometry::Geometry
     *
     *  Every implementation of a DUNE Geometry is required to have the same
     *  template parameter list:
     *  \code
     *  template< int mydim, int cdim, class Grid >
     *  \endcode
     *  Consequently, there is no direct way to pass compile time static
     *  information to a unified implementation such as the generic geometries.
     *  The structure GeometryTraits realizes an indirect way to do this.
     *
     *  For every grid implementation using the generic geometries, this
     *  structure must be specialized. The following default implementation
     *  can be used (via derivation) to provide the necessary information. It
     *  contains exactly the fields that are necessary:
     *  \code
     *  template< class ctype, int dimG, int dimW >
     *  struct DefaultGeometryTraits
     *  {
     *    ...
     *  };
     *  \endcode
     */
    template< class Grid >
    struct GlobalGeometryTraits;

    template< class Grid >
    struct GlobalGeometryTraits< const Grid >
      : public GlobalGeometryTraits< Grid >
    {};



    /** \struct  LocalGeometryTraits
     *  \ingroup GenericGeometry
     *  \brief   grid specific information required by GenericGeometry::LocalGeometry
     *
     *  Every implementation of a DUNE Geometry is required to have the same
     *  template parameter list:
     *  \code
     *  template< int mydim, int cdim, class Grid >
     *  \endcode
     *  Consequently, there is no direct way to pass compile time static
     *  information to a unified implementation such as the generic geometries.
     *  The structure GeometryTraits realizes an indirect way to do this.
     *
     *  For every grid implementation using the generic geometries, this
     *  structure must be specialized. The following default implementation
     *  can be used (via derivation) to provide the necessary information. It
     *  contains exactly the fields that are necessary:
     *  \code
     *  template< class ctype, int dimG, int dimW >
     *  struct DefaultGeometryTraits
     *  {
     *    ...
     *  };
     *  \endcode
     */
    template< class Grid >
    struct LocalGeometryTraits;

    template< class Grid >
    struct LocalGeometryTraits< const Grid >
      : public LocalGeometryTraits< Grid >
    {};
  }

}

#endif
