// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_NATIVEITERATOR_HH
#define DUNE_ALBERTA_NATIVEITERATOR_HH

#include <dune/grid/albertagrid/elementinfo.hh>

namespace Dune
{

  namespace Alberta
  {

    // MacroIterator
    // -------------

#if DUNE_ALBERTA_VERSION >= 0x200
    template< int dim >
    class MacroIterator
    {
      Mesh *mesh_;
      int index_;

    public:
      typedef Alberta::ElementInfo< dim > ElementInfo;

      explicit MacroIterator ( Mesh &mesh, bool end = false )
        : mesh_( &mesh ),
          index_( end ? mesh.n_macro_el : 0 )
      {}

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
        return (index_ >= mesh().n_macro_el);
      }

      MacroElement &macroElement () const
      {
        assert( !done() );
        return mesh_->macro_els[ index_ ];
      }

      Mesh &mesh () const
      {
        return *mesh_;
      }
    };
#else
    template< int dim >
    class MacroIterator
    {
      Mesh *mesh_;
      MacroElement *element_;

    public:
      typedef Alberta::ElementInfo< dim > ElementInfo;

      explicit MacroIterator ( Mesh &mesh, bool end = false )
        : mesh_( &mesh ),
          element_( end ? NULL : mesh.first_macro_el )
      {}

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

      Mesh &mesh () const
      {
        return *mesh_;
      }
    };
#endif

  }

}

#endif
