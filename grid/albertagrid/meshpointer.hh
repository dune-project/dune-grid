// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_MESHPOINTER_HH
#define DUNE_ALBERTA_MESHPOINTER_HH

#include <string>

#include <dune/grid/albertagrid/elementinfo.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    // MeshPointer
    // -----------

    template< int dim >
    class MeshPointer
    {
      Mesh *mesh_;

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

      void create ( const std::string &name, const std::string &filename );

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
    };


#if DUNE_ALBERTA_VERSION >= 0x200
    typedef NODE_PROJECTION *InitBoundary ( Mesh *, MacroElement *, int );
    static InitBoundary *initBoundary = 0;

    template< int dim >
    inline void MeshPointer< dim >
    ::create ( const std::string &name, const std::string &filename )
    {
      release();

      MACRO_DATA *macro = read_macro( filename.c_str() );

#if DUNE_ALBERTA_VERSION >= 0x201
      mesh_ = GET_MESH( dim, name.c_str(), macro, NULL, NULL );
#else
      mesh_ = GET_MESH( dim, name.c_str(), macro, NULL );
#endif

      free_macro_data( macro );

      if( mesh_ != NULL )
      {
        AlbertHelp::initDofAdmin< dim >( mesh_ );

        typedef AlbertHelp::AlbertLeafData< dim, dim+1 > LeafData;
        init_leaf_data( mesh_, sizeof( typename LeafData::Data ),
                        LeafData::AlbertLeafRefine,
                        LeafData::AlbertLeafCoarsen );
      }
    }
#endif // #if DUNE_ABLERTA_VERSION >= 0x200

#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    inline void MeshPointer< dim >
    ::create ( const std::string &name, const std::string &filename )
    {
      release();

      typedef AlbertHelp::AlbertLeafData< dim, dim+1 > LeafData;
      mesh_ = get_mesh( name.c_str(), AlbertHelp::initDofAdmin< dim >, LeafData::initLeafData );
      if( mesh_ != NULL )
        read_macro( mesh_, filename.c_str(), BoundaryProvider::initBoundary );
    }
#endif // #if DUNE_ABLERTA_VERSION < 0x200



#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    inline void MeshPointer< dim >::read ( const std::string &filename, Real &time )
    {
      release();

      typedef AlbertHelp::AlbertLeafData< dim, dim+1 > LeafData;
      mesh_ = read_mesh_xdr( filename.c_str(), &time, LeafData::initLeafData, BoundaryProvider::initBoundary );
    }
#endif // #if DUNE_ABLERTA_VERSION < 0x200

#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    inline void MeshPointer< dim >::read ( const std::string &filename, Real &time )
    {
      release();
#if DUNE_ALBERTA_VERSION >= 0x201
      mesh_ = ALBERTA read_mesh_xdr( filename.c_str(), &time, NULL, NULL );
#else
      mesh_ = ALBERTA read_mesh_xdr( filename.c_str(), &time, NULL );
#endif

      if( mesh_ != NULL )
      {
        typedef AlbertHelp::AlbertLeafData< dim, dim+1 > LeafData;
        init_leaf_data( mesh_, sizeof( typename LeafData::Data ),
                        LeafData::AlbertLeafRefine,
                        LeafData::AlbertLeafCoarsen );
      }
    }
#endif // #if DUNE_ABLERTA_VERSION >= 0x200



    template< int dim >
    inline bool MeshPointer< dim >::write ( const std::string &filename, Real time ) const
    {
      int success = ALBERTA write_mesh_xdr( mesh_, filename.c_str(), time );
      return (success == 0);
    }



    // MeshPointer::MacroIterator
    // --------------------------

#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    class MeshPointer< dim >::MacroIterator
    {
      friend class MeshPointer< dim >;

    public:
      typedef Alberta::MeshPointer< dim > MeshPointer;
      typedef Alberta::ElementInfo< dim > ElementInfo;

    private:
      MeshPointer mesh_;
      int index_;

      explicit MacroIterator ( const MeshPointer &mesh, bool end = false )
        : mesh_( mesh ),
          index_( end ? numMacroElements() : 0 )
      {}

    public:
      MacroIterator &operator++ ()
      {
        assert( !done() );
        ++index_;
        return *this;
      }

      ElementInfo operator* () const
      {
        return (done() ? ElementInfo() : ElementInfo( mesh(), macroElement() ));
      }

      bool operator== ( const MacroIterator &other ) const
      {
        return (index_ == other.index_);
      }

      bool operator!= ( const MacroIterator &other ) const
      {
        return (index_ != other.index_);
      }

      bool done () const
      {
        return (index_ >= numMacroElements());
      }

      MacroElement &macroElement () const
      {
        assert( !done() );
        return mesh().mesh_->macro_els[ index_ ];
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
#endif // #if DUNE_ABLERTA_VERSION >= 0x200


#if DUNE_ALBERTA_VERSION < 0x200
    template< int dim >
    class MeshPointer< dim >::MacroIterator
    {
      friend class MeshPointer< dim >;

    public:
      typedef Alberta::MeshPointer< dim > MeshPointer;
      typedef Alberta::ElementInfo< dim > ElementInfo;

    private:
      MeshPointer mesh_;
      MacroElement *element_;

      explicit MacroIterator ( const MeshPointer &mesh, bool end = false )
        : mesh_( mesh ),
          element_( end ? NULL : mesh.mesh_->first_macro_el )
      {}

    public:
      MacroIterator &operator++ ()
      {
        assert( !done() );
        element_ = element_->next;
        return *this;
      }

      ElementInfo operator* () const
      {
        return (done() ? ElementInfo() : ElementInfo( mesh(), macroElement() ));
      }

      bool operator== ( const MacroIterator &other ) const
      {
        return (element_ == other.element_);
      }

      bool operator!= ( const MacroIterator &other ) const
      {
        return (element_ != other.element_);
      }

      bool done () const
      {
        return (element_ == NULL);
      }

      MacroElement &macroElement () const
      {
        assert( !done() );
        return *element_;
      }

      const MeshPointer &mesh () const
      {
        return mesh_;
      }
    };
#endif // #if DUNE_ABLERTA_VERSION < 0x200

  }

}

#endif // #if HAVE_ALBERTA

#endif
