// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_MESHPOINTER_HH
#define DUNE_ALBERTA_MESHPOINTER_HH

/** \file
 *  \author Martin Nolte
 *  \brief  provides a wrapper for ALBERTA's mesh structure
 */

#include <limits>
#include <string>

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/albertagrid/macrodata.hh>
#include <dune/grid/albertagrid/projection.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    // External Forward Declarations
    // -----------------------------

    template< int dim >
    class HierarchyDofNumbering;

    // MeshPointer
    // -----------

    template< int dim >
    class MeshPointer
    {
      typedef Alberta::ElementInfo< dim > ElementInfo;
      typedef typename ElementInfo::MacroElement MacroElement;
      typedef typename ElementInfo::FillFlags FillFlags;

      class BoundaryProvider;

      template< int dimWorld >
      struct Library;

    public:
      class MacroIterator;

      MeshPointer ()
        : mesh_( 0 )
      {}

      explicit MeshPointer ( Mesh *mesh )
        : mesh_( mesh )
      {}

      operator Mesh * () const
      {
        return mesh_;
      }

      explicit operator bool () const
      {
        return (bool)mesh_;
      }

      MacroIterator begin () const
      {
        return MacroIterator( *this, false );
      }

      MacroIterator end () const
      {
        return MacroIterator( *this, true );
      }

      int numMacroElements () const;
      int size ( int codim ) const;

      // create a mesh from a macrodata structure
      // params:  macroData - macro data structure
      // returns: number of boundary segments
      unsigned int create ( const MacroData< dim > &macroData );

      // create a mesh from a macrodata structure, adding projections
      // params:  macroData         - macro data structure
      //          projectionFactory - factory for the projections
      // returns: number of boundary segments
      template< class Proj, class Impl >
      unsigned int create ( const MacroData< dim > &macroData,
                            const ProjectionFactoryInterface< Proj, Impl > &projectionFactory );

      // create a mesh from a file
      // params:  filename - file name of an Alberta macro triangulation
      //          binary   - read binary?
      // returns: number of boundary segments
      unsigned int create ( const std::string &filename, bool binary = false );

      // read back a mesh from a file
      // params:  filename - file name of an Alberta save file
      //          time     - variable to receive the time stored in the file
      // returns: number of boundary segments
      //
      // notes: - projections are not preserved
      //        - we assume that projections are added in the same order they
      //          inserted in when the grid was created (otherwise the boundary
      //          indices change)
      unsigned int read ( const std::string &filename, Real &time );

      bool write ( const std::string &filename, Real time ) const;

      void release ();

      template< class Functor >
      void hierarchicTraverse ( Functor &functor,
                                typename FillFlags::Flags fillFlags = FillFlags::standard ) const;

      template< class Functor >
      void leafTraverse ( Functor &functor,
                          typename FillFlags::Flags fillFlags = FillFlags::standard ) const;

      bool coarsen ( typename FillFlags::Flags fillFlags = FillFlags::nothing );

      bool refine ( typename FillFlags::Flags fillFlags = FillFlags::nothing );

    private:
      static ALBERTA NODE_PROJECTION *
      initNodeProjection ( [[maybe_unused]] Mesh *mesh, ALBERTA MACRO_EL *macroElement, int n );
      template< class ProjectionProvider >
      static ALBERTA NODE_PROJECTION *
      initNodeProjection ( Mesh *mesh, ALBERTA MACRO_EL *macroElement, int n );

      Mesh *mesh_;
    };



    // MeshPointer::Library
    // --------------------

    template< int dim >
    template< int dimWorld >
    struct MeshPointer< dim >::Library
    {
      typedef Alberta::MeshPointer< dim > MeshPointer;

      static inline unsigned int boundaryCount = 0;
      static inline const void *projectionFactory = nullptr;

      static void
      create ( MeshPointer &ptr, const MacroData< dim > &macroData,
               ALBERTA NODE_PROJECTION *(*initNodeProjection)( Mesh *, ALBERTA MACRO_EL *, int ) );
      static void release ( MeshPointer &ptr );
    };



    // MeshPointer::MacroIterator
    // --------------------------

    template< int dim >
    class MeshPointer< dim >::MacroIterator
    {
      typedef MacroIterator This;

      friend class MeshPointer< dim >;

    public:
      typedef Alberta::MeshPointer< dim > MeshPointer;
      typedef Alberta::ElementInfo< dim > ElementInfo;

      MacroIterator ()
        : mesh_(),
          index_( -1 )
      {}

    private:

      explicit MacroIterator ( const MeshPointer &mesh, bool end = false )
        : mesh_( mesh ),
          index_( end ? mesh.numMacroElements() : 0 )
      {}

    public:
      bool done () const
      {
        return (index_ >= mesh().numMacroElements());
      }

      bool equals ( const MacroIterator &other ) const
      {
        return (index_ == other.index_);
      }

      void increment ()
      {
        assert( !done() );
        ++index_;
      }

      const MacroElement &macroElement () const
      {
        assert( !done() );
        return static_cast< const MacroElement & >( mesh().mesh_->macro_els[ index_ ] );
      }

      const MeshPointer &mesh () const
      {
        return mesh_;
      }

      This &operator++ ()
      {
        increment();
        return *this;
      }

      ElementInfo operator* () const
      {
        return elementInfo();
      }

      bool operator== ( const MacroIterator &other ) const
      {
        return equals( other );
      }

      bool operator!= ( const MacroIterator &other ) const
      {
        return !equals( other );
      }

      ElementInfo
      elementInfo ( typename FillFlags::Flags fillFlags = FillFlags::standard ) const
      {
        if( done() )
          return ElementInfo();
        else
          return ElementInfo( mesh(), macroElement(), fillFlags );
      }

    private:
      MeshPointer mesh_;
      int index_;
    };



    // Implementation of MeshPointer
    // -----------------------------

    template< int dim >
    inline int MeshPointer< dim >::numMacroElements () const
    {
      return (mesh_ ? mesh_->n_macro_el : 0);
    }


    template<>
    inline int MeshPointer< 1 >::size( int codim ) const
    {
      assert( (codim >= 0) && (codim <= 1) );
      return (codim == 0 ? mesh_->n_elements : mesh_->n_vertices);
    }

    template<>
    inline int MeshPointer< 2 >::size( int codim ) const
    {
      assert( (codim >= 0) && (codim <= 2) );
      if( codim == 0 )
        return mesh_->n_elements;
      else if( codim == 2 )
        return mesh_->n_vertices;
      else
        return mesh_->n_edges;
    }

    template<>
    inline int MeshPointer< 3 >::size( int codim ) const
    {
      assert( (codim >= 0) && (codim <= 3) );
      if( codim == 0 )
        return mesh_->n_elements;
      else if( codim == 3 )
        return mesh_->n_vertices;
      else if( codim == 1 )
        return mesh_->n_faces;
      else
        return mesh_->n_edges;
    }


    template< int dim >
    inline unsigned int MeshPointer< dim >
    ::create ( const MacroData< dim > &macroData )
    {
      release();

      Library< dimWorld >::boundaryCount = 0;
      Library< dimWorld >::create( *this, macroData, &initNodeProjection );
      return Library< dimWorld >::boundaryCount;
    }


    template< int dim >
    template< class Proj, class Impl >
    inline unsigned int MeshPointer< dim >
    ::create ( const MacroData< dim > &macroData,
               const ProjectionFactoryInterface< Proj, Impl > &projectionFactory )
    {
      typedef ProjectionFactoryInterface< Proj, Impl > ProjectionFactory;

      release();

      Library< dimWorld >::boundaryCount = 0;
      Library< dimWorld >::projectionFactory = &projectionFactory;
      Library< dimWorld >::create( *this, macroData, &initNodeProjection< ProjectionFactory > );
      Library< dimWorld >::projectionFactory = nullptr;
      return Library< dimWorld >::boundaryCount;
    }




    template< int dim >
    inline unsigned int MeshPointer< dim >
    ::create ( const std::string &filename, bool binary )
    {
      MacroData< dim > macroData;
      macroData.read( filename, binary );
      const unsigned int boundaryCount = create( macroData );
      macroData.release();
      return boundaryCount;
    }


    template< int dim >
    inline unsigned int MeshPointer< dim >::read ( const std::string &filename, Real &time )
    {
      release();

      Library< dimWorld >::boundaryCount = 0;
      mesh_ = ALBERTA read_mesh_xdr( filename.c_str(), &time, NULL, NULL );
      return Library< dimWorld >::boundaryCount;
    }


    template< int dim >
    inline bool MeshPointer< dim >::write ( const std::string &filename, Real time ) const
    {
      int success = ALBERTA write_mesh_xdr( mesh_, filename.c_str(), time );
      return (success == 0);
    }


    template< int dim >
    inline void MeshPointer< dim >::release ()
    {
      Library< dimWorld >::release( *this );
    }


    template< int dim >
    template< class Functor >
    inline void MeshPointer< dim >
    ::hierarchicTraverse ( Functor &functor,
                           typename FillFlags::Flags fillFlags ) const
    {
      const MacroIterator eit = end();
      for( MacroIterator it = begin(); it != eit; ++it )
      {
        const ElementInfo info = it.elementInfo( fillFlags );
        info.hierarchicTraverse( functor );
      }
    }


    template< int dim >
    template< class Functor >
    inline void MeshPointer< dim >
    ::leafTraverse ( Functor &functor,
                     typename FillFlags::Flags fillFlags ) const
    {
      const MacroIterator eit = end();
      for( MacroIterator it = begin(); it != eit; ++it )
      {
        const ElementInfo info = it.elementInfo( fillFlags );
        info.leafTraverse( functor );
      }
    }


    template< int dim >
    inline bool MeshPointer< dim >::coarsen ( typename FillFlags::Flags fillFlags )
    {
      const bool coarsened = (ALBERTA coarsen( mesh_, fillFlags ) == meshCoarsened);
      if( coarsened )
        ALBERTA dof_compress( mesh_ );
      return coarsened;
    }

    template< int dim >
    inline bool MeshPointer< dim >::refine ( typename FillFlags::Flags fillFlags )
    {
      return (ALBERTA refine( mesh_, fillFlags ) == meshRefined);
    }


    template< int dim >
    inline ALBERTA NODE_PROJECTION *
    MeshPointer< dim >::initNodeProjection ( Mesh *mesh, ALBERTA MACRO_EL *macroEl, int n )
    {
      const MacroElement &macroElement = static_cast< const MacroElement & >( *macroEl );
      if( (n > 0) && macroElement.isBoundary( n-1 ) )
        return new BasicNodeProjection( Library< dimWorld >::boundaryCount++ );
      else
        return 0;
    }


    template< int dim >
    template< class ProjectionFactory >
    inline ALBERTA NODE_PROJECTION *
    MeshPointer< dim >::initNodeProjection ( Mesh *mesh, ALBERTA MACRO_EL *macroEl, int n )
    {
      typedef typename ProjectionFactory::Projection Projection;

      const MacroElement &macroElement = static_cast< const MacroElement & >( *macroEl );

      MeshPointer< dim > meshPointer( mesh );
      ElementInfo elementInfo( meshPointer, macroElement, FillFlags::standard );
      const ProjectionFactory &projectionFactory = *static_cast< const ProjectionFactory * >( Library< dimWorld >::projectionFactory );
      if( (n > 0) && macroElement.isBoundary( n-1 ) )
      {
        const unsigned int boundaryIndex = Library< dimWorld >::boundaryCount++;
        if( projectionFactory.hasProjection( elementInfo, n-1 ) )
        {
          Projection projection = projectionFactory.projection( elementInfo, n-1 );
          return new NodeProjection< dim, Projection >( boundaryIndex, projection );
        }
        else
          return new BasicNodeProjection( boundaryIndex );
      }
      else if( (dim < dimWorld) && (n == 0) )
      {
        const unsigned int boundaryIndex = std::numeric_limits< unsigned int >::max();
        if( projectionFactory.hasProjection( elementInfo ) )
        {
          Projection projection = projectionFactory.projection( elementInfo );
          return new NodeProjection< dim, Projection >( boundaryIndex, projection );
        }
        else
          return 0;
      }
      else
        return 0;
    }

  } // namespace Alberta

} // namespace Dune

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_MESHPOINTER_HH
