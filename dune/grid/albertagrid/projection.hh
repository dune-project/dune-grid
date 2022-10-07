// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_NODEPROJECTION_HH
#define DUNE_ALBERTA_NODEPROJECTION_HH

#include <memory>

#include <dune/grid/common/boundaryprojection.hh>

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/elementinfo.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Proj, class Impl >
    class ProjectionFactory;



    // DuneBoundaryProjection
    // ----------------------

    template< int dim >
    class DuneBoundaryProjection
    {
      typedef DuneBoundaryProjection< dim > This;

    public:
      static const int dimension = dim;

      typedef Alberta::ElementInfo< dimension > ElementInfo;
      typedef FieldVector< Real, dimWorld > GlobalCoordinate;

      typedef Dune::DuneBoundaryProjection< dimWorld > Projection;
      typedef std::shared_ptr< const Projection > ProjectionPtr;

      explicit DuneBoundaryProjection ( const ProjectionPtr &projection )
        : projection_( projection )
      {}

      // note: GlobalVector is an array type; global is the return value
      void operator() ( const ElementInfo & /* elementInfo */, const LocalVector /* local */,
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
      ProjectionPtr projection_;
    };



    // ProjectionFactoryInterface
    // --------------------------

    template< class Proj, class Impl >
    class ProjectionFactoryInterface
    {
      typedef ProjectionFactoryInterface< Proj, Impl > This;

      friend class ProjectionFactory< Proj, Impl >;

    public:
      typedef Proj Projection;

      static const int dimension = Projection::dimension;

      typedef Alberta::ElementInfo< dimension > ElementInfo;

    private:
      ProjectionFactoryInterface ()
      {}

      ProjectionFactoryInterface ( const This &other );
      This &operator= ( const This &other );

    public:
      bool hasProjection ( const ElementInfo &elementInfo, const int face ) const
      {
        return asImpl().hasProjection( elementInfo, face );
      }

      bool hasProjection ( const ElementInfo &elementInfo ) const
      {
        return asImpl().hasProjection( elementInfo );
      }

      Projection projection ( const ElementInfo &elementInfo, const int face ) const
      {
        return asImpl().projection( elementInfo, face );
      };

      Projection projection ( const ElementInfo &elementInfo ) const
      {
        return asImpl().projection( elementInfo );
      };

    protected:
      const Impl &asImpl () const
      {
        return static_cast< const Impl & >( *this );
      }
    };



    // ProjectionFactory
    // -----------------

    template< class Proj, class Impl >
    class ProjectionFactory
      : public ProjectionFactoryInterface< Proj, Impl >
    {
      typedef ProjectionFactory< Proj, Impl > This;
      typedef ProjectionFactoryInterface< Proj, Impl > Base;

    public:
      typedef typename Base::Projection Projection;
      typedef typename Base::ElementInfo ElementInfo;

    protected:
      ProjectionFactory ()
      {}

    private:
      bool hasProjection ( const ElementInfo &elementInfo, const int face ) const;
      bool hasProjection ( const ElementInfo &elementInfo ) const;

      Projection projection ( const ElementInfo &elementInfo, const int face ) const;
      Projection projection ( const ElementInfo &elementInfo ) const;
    };



    // DuneGlobalBoundaryProjectionFactory
    // -----------------------------------

    template< int dim >
    class DuneGlobalBoundaryProjectionFactory
      : public ProjectionFactory< DuneBoundaryProjection< dim >, DuneGlobalBoundaryProjectionFactory< dim > >
    {
      typedef DuneGlobalBoundaryProjectionFactory< dim > This;
      typedef ProjectionFactory< DuneBoundaryProjection< dim >, This > Base;

    public:
      typedef typename Base::Projection Projection;
      typedef typename Base::ElementInfo ElementInfo;

      typedef typename Projection::ProjectionPtr DuneProjectionPtr;

      DuneGlobalBoundaryProjectionFactory ( const DuneProjectionPtr &projection )
        : projection_( projection )
      {}

      bool hasProjection ( const ElementInfo &elementInfo, const int face ) const
      {
        return true;
      }

      bool hasProjection ( const ElementInfo &elementInfo ) const
      {
        return true;
      }

      Projection projection ( const ElementInfo &elementInfo, const int face ) const
      {
        return projection_;
      };

      Projection projection ( const ElementInfo &elementInfo ) const
      {
        return projection_;
      };

    private:
      const Projection projection_;
    };



    // BasicNodeProjection
    // -------------------

    struct BasicNodeProjection
      : public ALBERTA NODE_PROJECTION
    {
      explicit BasicNodeProjection ( unsigned int boundaryIndex )
        : boundaryIndex_( boundaryIndex )
      {
        func = 0;
      }

      virtual ~BasicNodeProjection ()
      {}

      unsigned int boundaryIndex () const
      {
        return boundaryIndex_;
      }

    private:
      unsigned int boundaryIndex_;
    };



    // NodeProjection
    // --------------

    template< int dim, class Projection >
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

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_NODEPROJECTION_HH
