// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_NODEPROJECTION_HH
#define DUNE_ALBERTA_NODEPROJECTION_HH

#include <dune/grid/common/boundaryprojection.hh>

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/elementinfo.hh>

namespace Dune
{

  namespace Alberta
  {

    // NoProjection
    // ------------

    template< int dim >
    struct NoProjection
    {
      static const int dimension = dim;

      typedef Alberta::ElementInfo< dimension > ElementInfo;

      // note: GlobalVector is an array type; global is the return value
      void operator() ( const ElementInfo &elementInfo, const LocalVector local,
                        GlobalVector global ) const
      {}
    };



    // DuneBoundaryProjection
    // ----------------------

    template< int dim >
    struct DuneBoundaryProjection
    {
      static const int dimension = dim;

      typedef Alberta::ElementInfo< dimension > ElementInfo;
      typedef FieldVector< Real, dimWorld > GlobalCoordinate;

      typedef Dune::DuneBoundaryProjection< dimWorld > Projection;

      DuneBoundaryProjection ( const Projection &projection )
        : projection_( &projection )
      {}

      // note: GlobalVector is an array type; global is the return value
      void operator() ( const ElementInfo &elementInfo, const LocalVector local,
                        GlobalVector global ) const
      {
        GlobalCoordinate x;
        for( int i = 0; i < dimWorld; ++i )
          x[ i ] = global[ i ];
        GlobalCoordinate y = projection() ( x );
        for( int i = 0; i < dimWorld; ++i )
          global[ i ] = y[ i ];
      }

      const Projection &projection () const
      {
        return *projection_;
      }

    private:
      const Projection *projection_;
    };



    // DuneGlobalBoundaryProjectionFactory
    // -----------------------------------

    template< int dim >
    struct DuneGlobalBoundaryProjectionFactory
    {
      static const int dimension = dim;

      typedef DuneBoundaryProjection< dimension > Projection;

      typedef Alberta::ElementInfo< dimension > ElementInfo;

      DuneGlobalBoundaryProjectionFactory ( const Projection &projection )
        : projection_( projection )
      {}

      Projection projection ( const ElementInfo &elementInfo, const int face ) const
      {
        return Projection( projection_ );
      };

    private:
      const Projection &projection_;
    };



    // BasicNodeProjection
    // -------------------

    struct BasicNodeProjection
      : public ALBERTA NODE_PROJECTION
    {
      unsigned int boundaryIndex_;

      explicit BasicNodeProjection ( unsigned int boundaryIndex )
        : boundaryIndex_( boundaryIndex )
      {
        func = 0;
      }
    };




    // NodeProjection
    // --------------

    template< int dim, class Projection = NoProjection< dim > >
    class NodeProjection
      : public BasicNodeProjection
    {
      typedef NodeProjection< dim, Projection > This;
      typedef BasicNodeProjection Base;

    public:
      static const int dimension = dim;

      typedef Alberta::ElementInfo< dimension > ElementInfo;

    private:
      Projection projection_;

    public:
      NodeProjection ( unsigned int boundaryIndex, const Projection &projection )
        : Base( boundaryIndex ),
          projection_( projection )
      {
        func = apply;
      }

    private:
      // note: global is the return type (it is an array type and hence no
      //       reference is needed)
      static void apply ( GlobalVector global, const EL_INFO *info, const LocalVector local )
      {
        const ElementInfo elementInfo = ElementInfo::createFake( *info );

        assert( (info->fill_flag & FillFlags< dimension >::projection) != 0 );
        const This *nodeProjection = static_cast< const This * >( info->active_projection );

        assert( nodeProjection != NULL );
        nodeProjection->projection_( elementInfo, local, global );
      }
    };

  }

}

#endif // #ifndef DUNE_ALBERTA_NODEPROJECTION_HH
