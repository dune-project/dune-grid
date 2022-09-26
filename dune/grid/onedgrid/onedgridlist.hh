// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONEDGRID_LIST_HH
#define DUNE_ONEDGRID_LIST_HH

#include <dune/common/iteratorfacades.hh>

namespace Dune {
  /** \file
      \brief A simple doubly-linked list needed in OneDGrid
      \todo I'd love to get rid of this and use std::list instead.
      Unfortunately, there are problems.  I need to store pointers/iterators
      within one element which point to another element (e.g. the element father).
   */
  template<class T>
  class OneDGridListIterator
    : public BidirectionalIteratorFacade<OneDGridListIterator<T>,T>
  {
  public:
    bool equals(const OneDGridListIterator& other) const {
      return pointer_ == other.pointer_;
    }

    T& dereference() {
      return *pointer_;
    }

    void increment() {
      pointer_ = pointer_->succ_;
    }

    void decrement() {
      pointer_ = pointer_->pred_;
    }

    OneDGridListIterator() {}

    OneDGridListIterator(T* pointer) {
      pointer_ = pointer;
    }

    OneDGridListIterator operator=(T* pointer) {
      pointer_ = pointer;
    }

    operator T*() {return pointer_;}

  private:
    T* pointer_;
  };

  template<class T>
  class OneDGridList
  {

  public:
    typedef T* iterator;
    typedef const T* const_iterator;

    OneDGridList() : numelements(0), begin_(0), rbegin_(0) {}

    int size() const {return numelements;}

    iterator push_back (const T& value) {

      T* i = rbegin();

      // New list element by copy construction
      T* t = new T(value);

      // einfuegen
      if (begin_==0) {
        // insert in empty list
        begin_ = t;
        rbegin_ = t;
      }
      else
      {
        // insert after element i
        t->pred_ = i;
        t->succ_ = i->succ_;
        i->succ_ = t;

        if (t->succ_!=0)
          t->succ_->pred_ = t;

        // new tail?
        if (rbegin_==i)
          rbegin_ = t;
      }

      // adjust size, return iterator
      numelements = numelements+1;
      return t;
    }

    iterator insert (iterator i, const T& value) {

      // Insert before 'one-after-the-end' --> append to the list
      if (i==end())
        return push_back(value);

      // New list element by copy construction
      T* t = new T(value);

      // insert
      if (begin_==0)
      {
        // insert in empty list
        begin_=t;
        rbegin_=t;
      }
      else
      {
        // insert before element i
        t->succ_ = i;
        t->pred_ = i->pred_;
        i->pred_ = t;

        if (t->pred_!=0)
          t->pred_->succ_ = t;
        // new head?
        if (begin_==i)
          begin_ = t;
      }

      // adjust size, return iterator
      numelements = numelements+1;
      return t;
    }

    void erase (iterator& i)
    {
      // test argument
      if (i==0)
        return;

      // handle pointers to predecessor and succesor
      if (i->succ_!=0)
        i->succ_->pred_ = i->pred_;
      if (i->pred_!=0)
        i->pred_->succ_ = i->succ_;

      // head & tail
      if (begin_==i)
        begin_=i->succ_;
      if (rbegin_==i)
        rbegin_ = i->pred_;

      // adjust size
      numelements = numelements-1;

      // Actually delete the object
      delete(i);
    }

    iterator begin() {
      return begin_;
    }

    const_iterator begin() const {
      return begin_;
    }

    iterator end() {
      return NULL;
    }

    const_iterator end() const {
      return NULL;
    }

    iterator rbegin() {
      return rbegin_;
    }

    const_iterator rbegin() const {
      return rbegin_;
    }

  private:

    int numelements;

    T* begin_;
    T* rbegin_;

  };   // end class OneDGridList

} // namespace Dune

#endif
