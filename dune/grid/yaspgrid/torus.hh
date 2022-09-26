// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_TORUS_HH
#define DUNE_GRID_YASPGRID_TORUS_HH

#include <array>
#include <bitset>
#include <cmath>
#include <deque>
#include <iostream>
#include <vector>

#if HAVE_MPI
#include <mpi.h>
#endif

#include <dune/common/binaryfunctions.hh>
#include <dune/common/streamoperators.hh>
#include <dune/grid/common/exceptions.hh>

#include "partitioning.hh"

/** \file
 *  \brief This file provides the infrastructure for toroidal communication in YaspGrid.
 */

namespace Dune
{

  /*! Torus provides all the functionality to handle a toroidal communication structure:

     - Map a set of processes (given by an MPI communicator) to a torus of dimension d. The "optimal"
     torus dimensions are determined by a coarse mesh. The maximum side length is minimized.

     - Provide lists of neighboring processes and a method for nearest neighbor exchange
     using asynchronous communication with MPI. The periodic case is handled where one process
     might have to exchange several messages with the same process. (Logically, a process has always
     \f$3^d-1\f$ neighbors, but several of these logical neighbors might be identical)

     - Provide means to partition a grid to the torus.

   */
  template<class Communication, int d>
  class Torus {
  public:
    //! type used to pass tupels in and out
    typedef std::array<int, d> iTupel;


  private:
    struct CommPartner {
      int rank;
      iTupel delta;
      int index;
    };

    struct CommTask {
      int rank;      // process to send to / receive from
      int size;      // size of buffer
      void *buffer;  // buffer to send / receive
    };

  public:
    //! constructor making uninitialized object
    Torus ()
    {}

    //! make partitioner from communicator and coarse mesh size
    Torus (Communication comm, int tag, iTupel size, const YLoadBalance<d>* lb)
      : _comm(comm), _tag(tag)
    {
      // determine dimensions
      lb->loadbalance(size, _comm.size(), _dims);

      // compute increments for lexicographic ordering
      int inc = 1;
      for (int i=0; i<d; i++)
      {
        _increment[i] = inc;
        inc *= _dims[i];
      }

      // check whether the load balancing matches the size of the communicator
      if (inc != _comm.size())
        DUNE_THROW(Dune::Exception, "Communicator size and result of the given load balancer do not match!");

      // make full schedule
      proclists();
    }

    //! return own rank
    int rank () const
    {
      return _comm.rank();
    }

    //! return own coordinates
    iTupel coord () const
    {
      return rank_to_coord(_comm.rank());
    }

    //! return number of processes
    int procs () const
    {
      return _comm.size();
    }

    //! return dimensions of torus
    const iTupel & dims () const
    {
      return _dims;
    }

    //! return dimensions of torus in direction i
    int dims (int i) const
    {
      return _dims[i];
    }

    //! return communicator
    Communication comm () const
    {
      return _comm;
    }

    //! return tag used by torus
    int tag () const
    {
      return _tag;
    }

    //! return true if coordinate is inside torus
    bool inside (iTupel c) const
    {
      for (int i=d-1; i>=0; i--)
        if (c[i]<0 || c[i]>=_dims[i]) return false;
      return true;
    }

    //! map rank to coordinate in torus using lexicographic ordering
    iTupel rank_to_coord (int rank) const
    {
      iTupel coord;
      rank = rank%_comm.size();
      for (int i=d-1; i>=0; i--)
      {
        coord[i] = rank/_increment[i];
        rank = rank%_increment[i];
      }
      return coord;
    }

    //! map coordinate in torus to rank using lexicographic ordering
    int coord_to_rank (iTupel coord) const
    {
      for (int i=0; i<d; i++) coord[i] = coord[i]%_dims[i];
      int rank = 0;
      for (int i=0; i<d; i++) rank += coord[i]*_increment[i];
      return rank;
    }

    //! return rank of process where its coordinate in direction dir has offset cnt (handles periodic case)
    int rank_relative (int rank, int dir, int cnt) const
    {
      iTupel coord = rank_to_coord(rank);
      coord[dir] = (coord[dir]+_dims[dir]+cnt)%_dims[dir];
      return coord_to_rank(coord);
    }

    //! assign color to given coordinate
    int color (const iTupel & coord) const
    {
      int c = 0;
      int power = 1;

      // interior coloring
      for (int i=0; i<d; i++)
      {
        if (coord[i]%2==1) c += power;
        power *= 2;
      }

      // extra colors for boundary processes
      for (int i=0; i<d; i++)
      {
        if (_dims[i]>1 && coord[i]==_dims[i]-1) c += power;
        power *= 2;
      }

      return c;
    }

    //! assign color to given rank
    int color (int rank) const
    {
      return color(rank_to_coord(rank));
    }

    //! return the number of neighbors, which is \f$3^d-1\f$
    int neighbors () const
    {
      int n=1;
      for (int i=0; i<d; ++i)
        n *= 3;
      return n-1;
    }

    //! return true if neighbor with given delta is a neighbor under the given periodicity
    bool is_neighbor (iTupel delta, std::bitset<d> periodic) const
    {
      iTupel coord = rank_to_coord(_comm.rank()); // my own coordinate with 0 <= c_i < dims_i

      for (int i=0; i<d; i++)
      {
        if (delta[i]<0)
        {
          // if I am on the boundary and domain is not periodic => no neighbor
          if (coord[i]==0 && periodic[i]==false) return false;
        }
        if (delta[i]>0)
        {
          // if I am on the boundary and domain is not periodic => no neighbor
          if (coord[i]==_dims[i]-1 && periodic[i]==false) return false;
        }
      }
      return true;
    }

    /** \brief partition the given grid onto the torus and return the piece of the process with given rank; returns load imbalance
     * @param rank rank of our processor
     * @param origin_in global origin
     * @param size_in global size
     * @param origin_out origin of this processors interior
     * @param size_out size of this processors interior
     */
    double partition (int rank, iTupel origin_in, iTupel size_in, iTupel& origin_out, iTupel& size_out) const
    {
      iTupel coord = rank_to_coord(rank);
      double maxsize = 1;
      double sz = 1;

      // make a tensor product partition
      for (int i=0; i<d; i++)
      {
        // determine
        int m = size_in[i]/_dims[i];
        int r = size_in[i]%_dims[i];

        sz *= size_in[i];

        if (coord[i]<_dims[i]-r)
        {
          origin_out[i] = origin_in[i] + coord[i]*m;
          size_out[i] = m;
          maxsize *= m;
        }
        else
        {
          origin_out[i] = origin_in[i] + (_dims[i]-r)*m + (coord[i]-(_dims[i]-r))*(m+1);
          size_out[i] = m+1;
          maxsize *= m+1;
        }
      }
      return maxsize/(sz/_comm.size());
    }

    /*!
       ProcListIterator provides access to a list of neighboring processes. There are always
       \f$ 3^d-1 \f$ entries in such a list. Two lists are maintained, one for sending and one for
       receiving. The lists are sorted in such a way that in sequence message delivery ensures that
       e.g. a message send to the left neighbor is received as a message from the right neighbor.
     */
    class ProcListIterator {
    public:
      //! make an iterator
      ProcListIterator (typename std::deque<CommPartner>::const_iterator iter)
      {
        i = iter;
      }

      //! return rank of neighboring process
      int rank () const
      {
        return i->rank;
      }

      //! return distance vector
      iTupel delta () const
      {
        return i->delta;
      }

      //! return index in proclist
      int index () const
      {
        return i->index;
      }

      //! return 1-norm of distance vector
      int distance () const
      {
        int dist = 0;
        iTupel delta=i->delta;
        for (int j=0; j<d; ++j)
          dist += std::abs(delta[j]);
        return dist;
      }

      //! Return true when two iterators point to same member
      bool operator== (const ProcListIterator& iter) const
      {
        return i == iter.i;
      }


      //! Return true when two iterators do not point to same member
      bool operator!= (const ProcListIterator& iter) const
      {
        return i != iter.i;
      }

      //! Increment iterator to next cell.
      ProcListIterator& operator++ ()
      {
        ++i;
        return *this;
      }

    private:
      typename std::deque<CommPartner>::const_iterator i;
    };

    //! first process in send list
    ProcListIterator sendbegin () const
    {
      return ProcListIterator(_sendlist.begin());
    }

    //! end of send list
    ProcListIterator sendend () const
    {
      return ProcListIterator(_sendlist.end());
    }

    //! first process in receive list
    ProcListIterator recvbegin () const
    {
      return ProcListIterator(_recvlist.begin());
    }

    //! last process in receive list
    ProcListIterator recvend () const
    {
      return ProcListIterator(_recvlist.end());
    }

    //! store a send request; buffers are sent in order; handles also local requests with memcpy
    void send (int rank, void* buffer, int size) const
    {
      CommTask task;
      task.rank = rank;
      task.buffer = buffer;
      task.size = size;
      if (rank!=_comm.rank())
        _sendrequests.push_back(task);
      else
        _localsendrequests.push_back(task);
    }

    //! store a receive request; buffers are received in order; handles also local requests with memcpy
    void recv (int rank, void* buffer, int size) const
    {
      CommTask task;
      task.rank = rank;
      task.buffer = buffer;
      task.size = size;
      if (rank!=_comm.rank())
        _recvrequests.push_back(task);
      else
        _localrecvrequests.push_back(task);
    }

    //! exchange messages stored in request buffers; clear request buffers afterwards
    void exchange () const
    {
      // handle local requests first
      if (_localsendrequests.size()!=_localrecvrequests.size())
      {
        std::cout << "[" << rank() << "]: ERROR: local sends/receives do not match in exchange!" << std::endl;
        return;
      }
      for (unsigned int i=0; i<_localsendrequests.size(); i++)
      {
        if (_localsendrequests[i].size!=_localrecvrequests[i].size)
        {
          std::cout << "[" << rank() << "]: ERROR: size in local sends/receive does not match in exchange!" << std::endl;
          return;
        }
        memcpy(_localrecvrequests[i].buffer,_localsendrequests[i].buffer,_localsendrequests[i].size);
      }
      _localsendrequests.clear();
      _localrecvrequests.clear();

#if HAVE_MPI
      // handle foreign requests

      std::vector<MPI_Request> requests(_sendrequests.size() + _recvrequests.size());
      MPI_Request* req = requests.data();

      // issue sends to foreign processes
      for (unsigned int i=0; i<_sendrequests.size(); i++)
        if (_sendrequests[i].rank!=rank())
        {
          //          std::cout << "[" << rank() << "]" << " send " << _sendrequests[i].size << " bytes "
          //                    << "to " << _sendrequests[i].rank << " p=" << _sendrequests[i].buffer << std::endl;
          MPI_Isend(_sendrequests[i].buffer, _sendrequests[i].size, MPI_BYTE,
                    _sendrequests[i].rank, _tag, _comm, req++);
        }

      // issue receives from foreign processes
      for (unsigned int i=0; i<_recvrequests.size(); i++)
        if (_recvrequests[i].rank!=rank())
        {
          //          std::cout << "[" << rank() << "]"  << " recv " << _recvrequests[i].size << " bytes "
          //                    << "fm " << _recvrequests[i].rank << " p=" << _recvrequests[i].buffer << std::endl;
          MPI_Irecv(_recvrequests[i].buffer, _recvrequests[i].size, MPI_BYTE,
                    _recvrequests[i].rank, _tag, _comm, req++);
        }

      // Wait for communication to complete
      MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

      // clear request buffers
      _sendrequests.clear();
      _recvrequests.clear();
#endif
    }

    //! global max
    double global_max (double x) const
    {
      double res = 0.0;
      _comm.template allreduce<Dune::Max<double>,double>(&x, &res, 1);
      return res;
    }

    //! print contents of torus object
    void print (std::ostream& s) const
    {
      s << "[" << rank() <<  "]: Torus " << procs() << " processor(s) arranged as " << dims() << std::endl;
      for (ProcListIterator i=sendbegin(); i!=sendend(); ++i)
      {
        s << "[" << rank() <<  "]: send to   "
          << "rank=" << i.rank()
          << " index=" << i.index()
          << " delta=" << i.delta() << " dist=" << i.distance() << std::endl;
      }
      for (ProcListIterator i=recvbegin(); i!=recvend(); ++i)
      {
        s << "[" << rank() <<  "]: recv from "
          << "rank=" << i.rank()
          << " index=" << i.index()
          << " delta=" << i.delta() << " dist=" << i.distance() << std::endl;
      }
    }

  private:

    void proclists ()
    {
      // compile the full neighbor list
      CommPartner cp;
      iTupel delta;

      std::fill(delta.begin(), delta.end(), -1);
      bool ready = false;
      iTupel me, nb;
      me = rank_to_coord(_comm.rank());
      int index = 0;
      int last = neighbors()-1;
      while (!ready)
      {
        // find neighbors coordinates
        for (int i=0; i<d; i++)
          nb[i] = ( me[i]+_dims[i]+delta[i] ) % _dims[i];

        // find neighbors rank
        int nbrank = coord_to_rank(nb);

        // check if delta is not zero
        for (int i=0; i<d; i++)
          if (delta[i]!=0)
          {
            cp.rank = nbrank;
            cp.delta = delta;
            cp.index = index;
            _recvlist.push_back(cp);
            cp.index = last-index;
            _sendlist.push_front(cp);
            index++;
            break;
          }

        // next neighbor
        ready = true;
        for (int i=0; i<d; i++)
          if (delta[i]<1)
          {
            (delta[i])++;
            ready=false;
            break;
          }
          else
          {
            delta[i] = -1;
          }
      }

    }

    Communication _comm;

    iTupel _dims;
    iTupel _increment;
    int _tag;
    std::deque<CommPartner> _sendlist;
    std::deque<CommPartner> _recvlist;

    mutable std::vector<CommTask> _sendrequests;
    mutable std::vector<CommTask> _recvrequests;
    mutable std::vector<CommTask> _localsendrequests;
    mutable std::vector<CommTask> _localrecvrequests;

  };

  //! Output operator for Torus
  template <class Communication, int d>
  inline std::ostream& operator<< (std::ostream& s, const Torus<Communication, d> & t)
  {
    t.print(s);
    return s;
  }
}

#endif
