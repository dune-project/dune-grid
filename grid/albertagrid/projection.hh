// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_NODEPROJECTION_HH
#define DUNE_ALBERTA_NODEPROJECTION_HH

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/elementinfo.hh>

namespace Dune
{

  namespace Alberta
  {

    template< int dim >
    struct NoProjection
    {
      static const int dimension = dim;

      typedef Alberta::ElementInfo< dimension > ElementInfo;

      void operator() ( const ElementInfo &elementInfo, const LocalVector &local,
                        GlobalVector &global ) const
      {}
    };



    // NodeProjection
    // --------------

    template< int dim, class Projection = NoProjection< dim > >
    class NodeProjection
      : protected ALBERTA NODE_PROJECTION
    {
      typedef NodeProjection< dim, void > This;

    public:
      static const int dimension = dim;

      typedef Alberta::ElementInfo< dimension > ElementInfo;

    private:
      unsigned int boundaryIndex_;
      Projection projection_;

    public:
      NodeProjection ( unsigned int boundaryIndex, const Projection &projection )
        : boundaryIndex_( boundaryIndex ),
          projection_( projection )
      {
        func = apply;
      }

    private:
      // note: global is the return type (it is an array type and hence no
      //       reference is needed)
      static void apply ( GlobalVector global, const EL_INFO *info,
                          const LocalVector local )
      {
        const ElementInfo elementInfo = ElementInfo::createFake( info );

        assert( info->fill_flag & FillFlags< dimension >::projection != 0 );
        const This *nodeProjection = static_cast< const This * >( info->active_projection );

        assert( nodeProjection != NULL );
        nodeProjection->projection_( elementInfo, local, global );
      }
    };

  }

}

#endif // #ifndef DUNE_ALBERTA_NODEPROJECTION_HH
