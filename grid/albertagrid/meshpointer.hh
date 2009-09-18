// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_MESHPOINTER_HH
#define DUNE_ALBERTA_MESHPOINTER_HH

/** \file
 *  \author Martin Nolte
 *  \brief  provides a wrapper for ALBERTA's mesh structure
 */

#include <cassert>
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
      Mesh *mesh_;

      typedef Alberta::ElementInfo< dim > ElementInfo;
      typedef typename ElementInfo::MacroElement MacroElement;
      typedef typename ElementInfo::FillFlags FillFlags;

      class BoundaryProvider;
      class MacroIteratorBase;

      template< int dimWorld >
      struct Library
      {
        static unsigned int boundaryCount;
        static const void *projectionFactory;
      };

    public:
      class MacroIterator;

      MeshPointer ()
        : mesh_( NULL )
      {}

      explicit MeshPointer ( Mesh *mesh )
        : mesh_( mesh )
      {}

      operator Mesh * () const
      {
        return mesh_;
      }

      bool operator! () const
      {
        return (mesh_ == NULL);
      }

      MacroIterator begin () const
      {
        return MacroIterator( *this, false );
      }

      MacroIterator end () const
      {
        return MacroIterator( *this, true );
      }

      std::string name () const
      {
        if( mesh_ != NULL )
          return mesh_->name;
        else
          return std::string();
      }

      int size ( int codim ) const;

      void create ( const MacroData< dim > &macroData, const std::string &name );
      template< class ProjectionFactory >
      void create ( const MacroData< dim > &macroData, const std::string &name, const ProjectionFactory &projectionFactory );
      void create ( const std::string &filename, const std::string &name, bool binary = false );

      void read ( const std::string &filename, Real &time );

      bool write ( const std::string &filename, Real time ) const;

      void release ()
      {
        if( mesh_ != NULL )
        {
          ALBERTA free_mesh( mesh_ );
          mesh_ = NULL;
        }
      }

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
      initNodeProjection ( Mesh *mesh, ALBERTA MACRO_EL *macroElement, int n );
      template< class ProjectionProvider >
      static ALBERTA NODE_PROJECTION *
      initNodeProjection ( Mesh *mesh, ALBERTA MACRO_EL *macroElement, int n );
    };



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
    inline void MeshPointer< dim >
    ::create ( const MacroData< dim > &macroData, const std::string &name )
    {
      release();

      Library< dimWorld >::boundaryCount = 0;
#if DUNE_ALBERTA_VERSION >= 0x300
      mesh_ = GET_MESH( dim, name.c_str(), macroData, &initNodeProjection, NULL );
#else
      mesh_ = GET_MESH( dim, name.c_str(), macroData, &initNodeProjection );
#endif
    }


    template< int dim >
    template< class ProjectionFactory >
    inline void MeshPointer< dim >
    ::create ( const MacroData< dim > &macroData, const std::string &name, const ProjectionFactory &projectionFactory )
    {
      release();

      Library< dimWorld >::boundaryCount = 0;
      Library< dimWorld >::projectionFactory = &projectionFactory;
#if DUNE_ALBERTA_VERSION >= 0x300
      mesh_ = GET_MESH( dim, name.c_str(), macroData, &initNodeProjection< ProjectionFactory >, NULL );
#else
      mesh_ = GET_MESH( dim, name.c_str(), macroData, &initNodeProjection< ProjectionFactory > );
#endif
      Library< dimWorld >::projectionFactory = 0;
    }




    template< int dim >
    inline void MeshPointer< dim >
    ::create ( const std::string &filename, const std::string &name, bool binary )
    {
      MacroData< dim > macroData;
      macroData.read( filename, binary );
      create( macroData, name );
      macroData.release();
    }


    template< int dim >
    inline void MeshPointer< dim >::read ( const std::string &filename, Real &time )
    {
      release();
#if DUNE_ALBERTA_VERSION >= 0x300
      mesh_ = ALBERTA read_mesh_xdr( filename.c_str(), &time, NULL, NULL );
#else
      mesh_ = ALBERTA read_mesh_xdr( filename.c_str(), &time, NULL );
#endif
    }


    template< int dim >
    inline bool MeshPointer< dim >::write ( const std::string &filename, Real time ) const
    {
      int success = ALBERTA write_mesh_xdr( mesh_, filename.c_str(), time );
      return (success == 0);
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


#if DUNE_ALBERTA_VERSION >= 0x300
    template< int dim >
    inline bool MeshPointer< dim >::coarsen ( typename FillFlags::Flags fillFlags )
    {
      const bool coarsened = (ALBERTA coarsen( mesh_, fillFlags ) == meshCoarsened);
      if( coarsened )
        ALBERTA dof_compress( mesh_ );
      return coarsened;
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x300

#if DUNE_ALBERTA_VERSION <= 0x200
    template< int dim >
    inline bool MeshPointer< dim >::coarsen ( typename FillFlags::Flags fillFlags )
    {
      assert( fillFlags == FillFlags::nothing );
      const bool coarsened = (ALBERTA coarsen( mesh_ ) == meshCoarsened);
      if( coarsened )
        ALBERTA dof_compress( mesh_ );
      return coarsened;
    }
#endif // #if DUNE_ALBERTA_VERSION <= 0x200


#if DUNE_ALBERTA_VERSION >= 0x300
    template< int dim >
    inline bool MeshPointer< dim >::refine ( typename FillFlags::Flags fillFlags )
    {
      return (ALBERTA refine( mesh_, fillFlags ) == meshRefined);
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x300

#if DUNE_ALBERTA_VERSION <= 0x200
    template< int dim >
    inline bool MeshPointer< dim >::refine ( typename FillFlags::Flags fillFlags )
    {
      assert( fillFlags == FillFlags::nothing );
      return (ALBERTA refine( mesh_ ) == meshRefined);
    }
#endif // #if DUNE_ALBERTA_VERSION <= 0x200


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
      if( (n > 0) && macroElement.isBoundary( n-1 ) )
      {
        MeshPointer< dim > meshPointer( mesh );
        ElementInfo elementInfo( meshPointer, macroElement, FillFlags::standard );
        const ProjectionFactory &projectionFactory = *static_cast< const ProjectionFactory * >( Library< dimWorld >::projectionFactory );
        Projection projection = projectionFactory.projection( elementInfo, n-1 );
        return new NodeProjection< dim, Projection >( Library< dimWorld >::boundaryCount++, projection );
      }
      return 0;
    }



    // MeshPointer::MacroIteratorBase
    // ------------------------------

    template< int dim >
    class MeshPointer< dim >::MacroIteratorBase
    {
      friend class MacroIterator;

    public:
      typedef Alberta::MeshPointer< dim > MeshPointer;
      typedef Alberta::ElementInfo< dim > ElementInfo;

    private:
      MeshPointer mesh_;
      int index_;

      explicit MacroIteratorBase ( const MeshPointer &mesh, bool end = false )
        : mesh_( mesh ),
          index_( end ? numMacroElements() : 0 )
      {}

    public:
      bool done () const
      {
        return (index_ >= numMacroElements());
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

    private:
      int numMacroElements () const
      {
        return mesh().mesh_->n_macro_el;
      }
    };



    // MeshPointer::MacroIterator
    // --------------------------

    template< int dim >
    class MeshPointer< dim >::MacroIterator
      : public MeshPointer< dim >::MacroIteratorBase
    {
      typedef MacroIterator This;
      typedef MacroIteratorBase Base;

      friend class MeshPointer< dim >;

    public:
      typedef Alberta::MeshPointer< dim > MeshPointer;
      typedef Alberta::ElementInfo< dim > ElementInfo;

    private:
      explicit MacroIterator ( const MeshPointer &mesh, bool end = false )
        : Base( mesh, end )
      {}

    public:
      using Base::done;
      using Base::equals;
      using Base::increment;
      using Base::macroElement;
      using Base::mesh;

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
    };

  }

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_MESHPOINTER_HH
