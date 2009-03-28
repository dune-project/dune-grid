// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_NODEPROJECTION_HH
#define DUNE_ALBERTA_NODEPROJECTION_HH

#include <dune/grid/albertagrid/misc.hh>

namespace Dune
{

  namespace Alberta
  {

    // NodeProjection
    // --------------

#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    class NodeProjection
      : protected ALBERTA NODE_PROJECTION
    {
      typedef NodeProjection< dim > This;

    public:
      static const int dimension = dim;

      typedef Dune::BoundarySegment< dimension, dimWorld > BoundarySegment;

    private:
      const BoundarySegment *boundarySegment_;

    public:
      NodeProjection ( const BoundarySegment *const boundarySegment )
        : boundarySegment_( boundarySegment )
      {
        func = apply;
      }

      void operator() ( const ALBERTA EL_INFO *info, const LocalVector &local
                        const GlobalVector &global ) const
      {}

    private:
      // note: global is the return type (it is an array type and hence no
      //       reference is needed)
      static void apply ( GlobalVector global, const EL_INFO *info,
                          const LocalVector local )
      {
        assert( info->fill_flag & FillFlags< dimension >::projection != 0 );
        const This *projection = static_cast< const This * >( info->active_projection );
        assert( projection != NULL );
        (*projection)( global, info, local );
      }
    };
#endif // #if DUNE_ALBERTA_VERSION >= 0x200

  }

}

#endif // #ifndef DUNE_ALBERTA_NODEPROJECTION_HH
