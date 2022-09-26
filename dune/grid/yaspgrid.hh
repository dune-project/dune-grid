// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_HH
#define DUNE_GRID_YASPGRID_HH

#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <stack>
#include <type_traits>

#include <dune/grid/common/backuprestore.hh>
#include <dune/grid/common/grid.hh>     // the grid base classes
#include <dune/grid/common/capabilities.hh> // the capabilities
#include <dune/common/hybridutilities.hh>
#include <dune/common/bigunsignedint.hh>
#include <dune/common/math.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/indexidset.hh>
#include <dune/grid/common/datahandleif.hh>


#if HAVE_MPI
#include <dune/common/parallel/mpicommunication.hh>
#endif

/*! \file yaspgrid.hh
 * YaspGrid stands for yet another structured parallel grid.
 * It will implement the dune grid interface for structured grids
 * with arbitrary overlap, parallel features with two overlap
 * models, periodic boundaries and a fast implementation allowing on-the-fly computations.
 */

namespace Dune {

  /* some sizes for building global ids
   */
  const int yaspgrid_dim_bits = 24; // bits for encoding each dimension
  const int yaspgrid_level_bits = 5; // bits for encoding level number


  //************************************************************************
  // forward declaration of templates

  template<int dim, class Coordinates>                             class YaspGrid;
  template<int mydim, int cdim, class GridImp>  class YaspGeometry;
  template<int codim, int dim, class GridImp>   class YaspEntity;
  template<int codim, class GridImp>            class YaspEntitySeed;
  template<int codim, PartitionIteratorType pitype, class GridImp> class YaspLevelIterator;
  template<class GridImp>            class YaspIntersectionIterator;
  template<class GridImp>            class YaspIntersection;
  template<class GridImp>            class YaspHierarchicIterator;
  template<class GridImp, bool isLeafIndexSet>                     class YaspIndexSet;
  template<class GridImp>            class YaspGlobalIdSet;
  template<class GridImp>            class YaspPersistentContainerIndex;

} // namespace Dune

#include <dune/grid/yaspgrid/coordinates.hh>
#include <dune/grid/yaspgrid/torus.hh>
#include <dune/grid/yaspgrid/ygrid.hh>
#include <dune/grid/yaspgrid/yaspgridgeometry.hh>
#include <dune/grid/yaspgrid/yaspgridentity.hh>
#include <dune/grid/yaspgrid/yaspgridintersection.hh>
#include <dune/grid/yaspgrid/yaspgridintersectioniterator.hh>
#include <dune/grid/yaspgrid/yaspgridhierarchiciterator.hh>
#include <dune/grid/yaspgrid/yaspgridentityseed.hh>
#include <dune/grid/yaspgrid/yaspgridleveliterator.hh>
#include <dune/grid/yaspgrid/yaspgridindexsets.hh>
#include <dune/grid/yaspgrid/yaspgrididset.hh>
#include <dune/grid/yaspgrid/yaspgridpersistentcontainer.hh>

namespace Dune {

#if HAVE_MPI
  using YaspCommunication = Communication<MPI_Comm>;
#else
  using YaspCommunication = Communication<No_Comm>;
#endif

  template<int dim, class Coordinates>
  struct YaspGridFamily
  {
    typedef YaspCommunication CCType;

    typedef GridTraits<dim,                                     // dimension of the grid
        dim,                                                    // dimension of the world space
        Dune::YaspGrid<dim, Coordinates>,
        YaspGeometry,YaspEntity,
        YaspLevelIterator,                                      // type used for the level iterator
        YaspIntersection,              // leaf  intersection
        YaspIntersection,              // level intersection
        YaspIntersectionIterator,              // leaf  intersection iter
        YaspIntersectionIterator,              // level intersection iter
        YaspHierarchicIterator,
        YaspLevelIterator,                                      // type used for the leaf(!) iterator
        YaspIndexSet< const YaspGrid< dim, Coordinates >, false >,                  // level index set
        YaspIndexSet< const YaspGrid< dim, Coordinates >, true >,                  // leaf index set
        YaspGlobalIdSet<const YaspGrid<dim, Coordinates> >,
        bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+dim>,
        YaspGlobalIdSet<const YaspGrid<dim, Coordinates> >,
        bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+dim>,
        CCType,
        DefaultLevelGridViewTraits, DefaultLeafGridViewTraits,
        YaspEntitySeed>
    Traits;
  };

#ifndef DOXYGEN
  template<int dim, int codim>
  struct YaspCommunicateMeta {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,codim))
      {
        g.template communicateCodim<DataHandle,codim>(data,iftype,dir,level);
      }
      YaspCommunicateMeta<dim,codim-1>::comm(g,data,iftype,dir,level);
    }
  };

  template<int dim>
  struct YaspCommunicateMeta<dim,0> {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,0))
        g.template communicateCodim<DataHandle,0>(data,iftype,dir,level);
    }
  };
#endif

  //************************************************************************
  /*!
   * \brief [<em> provides \ref Dune::Grid </em>]
   * \brief Provides a distributed structured cube mesh.
   * \ingroup GridImplementations
   *
   * YaspGrid stands for yet another structured parallel grid.
   * It implements the dune grid interface for structured grids
   * with arbitrary overlap (including zero),
   * periodic boundaries, and a fast implementation allowing on-the-fly computations.
   *
   * YaspGrid supports three coordinate modes: \ref EquidistantCoordinates,
   * \ref EquidistantOffsetCoordinates, and \ref Dune::TensorProductCoordinates.
   *
   * \tparam dim The dimension of the grid and its surrounding world
   * \tparam Coordinates The coordinate mode of the grid.
   */
  template<int dim, class Coordinates = EquidistantCoordinates<double, dim> >
  class YaspGrid
    : public GridDefaultImplementation<dim,dim,typename Coordinates::ctype,YaspGridFamily<dim, Coordinates> >
  {

    template<int, PartitionIteratorType, typename>
    friend class YaspLevelIterator;

    template<typename>
    friend class YaspHierarchicIterator;

    using Base = GridDefaultImplementation<dim,dim,typename Coordinates::ctype,YaspGridFamily<dim, Coordinates> >;
  protected:

  public:
    //! Type used for coordinates
    typedef typename Coordinates::ctype ctype;

    using Communication = typename Base::Communication;
    using CollectiveCommunicationType [[deprecated("Use Communication. Will be removed after release 2.9")]] = typename Base::Communication;

#ifndef DOXYGEN
    typedef typename Dune::YGrid<Coordinates> YGrid;
    typedef typename Dune::YGridList<Coordinates>::Intersection Intersection;

    /** \brief A single grid level within a YaspGrid
     */
    struct YGridLevel {

      /** \brief Level number of this level grid */
      int level() const
      {
        return level_;
      }

      Coordinates coords;

      std::array<YGrid, dim+1> overlapfront;
      std::array<YGridComponent<Coordinates>, Dune::power(2,dim)> overlapfront_data;
      std::array<YGrid, dim+1> overlap;
      std::array<YGridComponent<Coordinates>, Dune::power(2,dim)> overlap_data;
      std::array<YGrid, dim+1> interiorborder;
      std::array<YGridComponent<Coordinates>, Dune::power(2,dim)> interiorborder_data;
      std::array<YGrid, dim+1> interior;
      std::array<YGridComponent<Coordinates>, Dune::power(2,dim)> interior_data;

      std::array<YGridList<Coordinates>,dim+1> send_overlapfront_overlapfront;
      std::array<std::deque<Intersection>, Dune::power(2,dim)>  send_overlapfront_overlapfront_data;
      std::array<YGridList<Coordinates>,dim+1> recv_overlapfront_overlapfront;
      std::array<std::deque<Intersection>, Dune::power(2,dim)>  recv_overlapfront_overlapfront_data;

      std::array<YGridList<Coordinates>,dim+1> send_overlap_overlapfront;
      std::array<std::deque<Intersection>, Dune::power(2,dim)>  send_overlap_overlapfront_data;
      std::array<YGridList<Coordinates>,dim+1> recv_overlapfront_overlap;
      std::array<std::deque<Intersection>, Dune::power(2,dim)>  recv_overlapfront_overlap_data;

      std::array<YGridList<Coordinates>,dim+1> send_interiorborder_interiorborder;
      std::array<std::deque<Intersection>, Dune::power(2,dim)>  send_interiorborder_interiorborder_data;
      std::array<YGridList<Coordinates>,dim+1> recv_interiorborder_interiorborder;
      std::array<std::deque<Intersection>, Dune::power(2,dim)>  recv_interiorborder_interiorborder_data;

      std::array<YGridList<Coordinates>,dim+1> send_interiorborder_overlapfront;
      std::array<std::deque<Intersection>, Dune::power(2,dim)>  send_interiorborder_overlapfront_data;
      std::array<YGridList<Coordinates>,dim+1> recv_overlapfront_interiorborder;
      std::array<std::deque<Intersection>, Dune::power(2,dim)>  recv_overlapfront_interiorborder_data;

      // general
      YaspGrid<dim,Coordinates>* mg;  // each grid level knows its multigrid
      int overlapSize;           // in mesh cells on this level
      bool keepOverlap;

      /** \brief The level number within the YaspGrid level hierarchy */
      int level_;
    };

    //! define types used for arguments
    typedef std::array<int, dim> iTupel;
    typedef FieldVector<ctype, dim> fTupel;

    // communication tag used by multigrid
    enum { tag = 17 };
#endif

    //! return reference to torus
    const Torus<Communication, dim>& torus () const
    {
      return _torus;
    }

    //! return number of cells on finest level in given direction on all processors
    int globalSize(int i) const
    {
      return levelSize(maxLevel(),i);
    }

    //! return number of cells on finest level on all processors
    iTupel globalSize() const
    {
      return levelSize(maxLevel());
    }

    //! return size of the grid (in cells) on level l in direction i
    int levelSize(int l, int i) const
    {
      return _coarseSize[i] * (1 << l);
    }

    //! return size vector of the grid (in cells) on level l
    iTupel levelSize(int l) const
    {
      iTupel s;
      for (int i=0; i<dim; ++i)
        s[i] = levelSize(l,i);
      return s;
    }

    //! return whether the grid is periodic in direction i
    bool isPeriodic(int i) const
    {
      return _periodic[i];
    }

    bool getRefineOption() const
    {
      return keep_ovlp;
    }

    //! Iterator over the grid levels
    typedef typename ReservedVector<YGridLevel,32>::const_iterator YGridLevelIterator;

    //! return iterator pointing to coarsest level
    YGridLevelIterator begin () const
    {
      return _levels.begin();
    }

    //! return iterator pointing to given level
    YGridLevelIterator begin (int i) const
    {
      if (i<0 || i>maxLevel())
        DUNE_THROW(GridError, "level not existing");
      return std::next(_levels.begin(), i);
    }

    //! return iterator pointing to one past the finest level
    YGridLevelIterator end () const
    {
      return _levels.end();
    }

    // static method to create the default load balance strategy
    static const YLoadBalanceDefault<dim>* defaultLoadbalancer()
    {
      static YLoadBalanceDefault<dim> lb;
      return & lb;
    }

  protected:
    /** \brief Make a new YGridLevel structure
     *
     * \param coords      the coordinate container
     * \param periodic    indicate periodicity for each direction
     * \param o_interior  origin of interior (non-overlapping) cell decomposition
     * \param overlap     to be used on this grid level
     */
    void makelevel (const Coordinates& coords, std::bitset<dim> periodic, iTupel o_interior, int overlap)
    {
      YGridLevel& g = _levels.back();
      g.overlapSize = overlap;
      g.mg = this;
      g.level_ = maxLevel();
      g.coords = coords;
      g.keepOverlap = keep_ovlp;

      using Dune::power;

      // set the inserting positions in the corresponding arrays of YGridLevelStructure
      typename std::array<YGridComponent<Coordinates>, power(2,dim)>::iterator overlapfront_it = g.overlapfront_data.begin();
      typename std::array<YGridComponent<Coordinates>, power(2,dim)>::iterator overlap_it = g.overlap_data.begin();
      typename std::array<YGridComponent<Coordinates>, power(2,dim)>::iterator interiorborder_it = g.interiorborder_data.begin();
      typename std::array<YGridComponent<Coordinates>, power(2,dim)>::iterator interior_it = g.interior_data.begin();

      typename std::array<std::deque<Intersection>, power(2,dim)>::iterator
        send_overlapfront_overlapfront_it = g.send_overlapfront_overlapfront_data.begin();
      typename std::array<std::deque<Intersection>, power(2,dim)>::iterator
        recv_overlapfront_overlapfront_it = g.recv_overlapfront_overlapfront_data.begin();

      typename std::array<std::deque<Intersection>, power(2,dim)>::iterator
        send_overlap_overlapfront_it = g.send_overlap_overlapfront_data.begin();
      typename std::array<std::deque<Intersection>, power(2,dim)>::iterator
        recv_overlapfront_overlap_it = g.recv_overlapfront_overlap_data.begin();

      typename std::array<std::deque<Intersection>, power(2,dim)>::iterator
        send_interiorborder_interiorborder_it = g.send_interiorborder_interiorborder_data.begin();
      typename std::array<std::deque<Intersection>, power(2,dim)>::iterator
        recv_interiorborder_interiorborder_it = g.recv_interiorborder_interiorborder_data.begin();

      typename std::array<std::deque<Intersection>, power(2,dim)>::iterator
        send_interiorborder_overlapfront_it = g.send_interiorborder_overlapfront_data.begin();
      typename std::array<std::deque<Intersection>, power(2,dim)>::iterator
        recv_overlapfront_interiorborder_it = g.recv_overlapfront_interiorborder_data.begin();

      // have a null array for constructor calls around
      std::array<int,dim> n;
      std::fill(n.begin(), n.end(), 0);

      // determine origin of the grid with overlap and store whether an overlap area exists in direction i.
      std::bitset<dim> ovlp_low(0ULL);
      std::bitset<dim> ovlp_up(0ULL);

      iTupel o_overlap;
      iTupel s_overlap;

      // determine at where we have overlap and how big the size of the overlap partition is
      for (int i=0; i<dim; i++)
      {
        // the coordinate container has been contructed to hold the entire grid on
        // this processor, including overlap. this is the element size.
        s_overlap[i] = g.coords.size(i);

        //in the periodic case there is always overlap
        if (periodic[i])
        {
          o_overlap[i] = o_interior[i]-overlap;
          ovlp_low[i] = true;
          ovlp_up[i] = true;
        }
        else
        {
          //check lower boundary
          if (o_interior[i] - overlap < 0)
            o_overlap[i] = 0;
          else
          {
            o_overlap[i] = o_interior[i] - overlap;
            ovlp_low[i] = true;
          }

          //check upper boundary
          if (o_overlap[i] + g.coords.size(i) < globalSize(i))
            ovlp_up[i] = true;
        }
      }

      for (unsigned int codim = 0; codim < dim + 1; codim++)
      {
        // set the begin iterator for the corresponding ygrids
        g.overlapfront[codim].setBegin(overlapfront_it);
        g.overlap[codim].setBegin(overlap_it);
        g.interiorborder[codim].setBegin(interiorborder_it);
        g.interior[codim].setBegin(interior_it);
        g.send_overlapfront_overlapfront[codim].setBegin(send_overlapfront_overlapfront_it);
        g.recv_overlapfront_overlapfront[codim].setBegin(recv_overlapfront_overlapfront_it);
        g.send_overlap_overlapfront[codim].setBegin(send_overlap_overlapfront_it);
        g.recv_overlapfront_overlap[codim].setBegin(recv_overlapfront_overlap_it);
        g.send_interiorborder_interiorborder[codim].setBegin(send_interiorborder_interiorborder_it);
        g.recv_interiorborder_interiorborder[codim].setBegin(recv_interiorborder_interiorborder_it);
        g.send_interiorborder_overlapfront[codim].setBegin(send_interiorborder_overlapfront_it);
        g.recv_overlapfront_interiorborder[codim].setBegin(recv_overlapfront_interiorborder_it);

        // find all combinations of unit vectors that span entities of the given codimension
        for (unsigned int index = 0; index < (1<<dim); index++)
        {
          // check whether the given shift is of our codimension
          std::bitset<dim> r(index);
          if (r.count() != dim-codim)
            continue;

          // get an origin and a size array for subsequent modification
          std::array<int,dim> origin(o_overlap);
          std::array<int,dim> size(s_overlap);

          // build overlapfront
          // we have to extend the element size by one in all directions without shift.
          for (int i=0; i<dim; i++)
            if (!r[i])
              size[i]++;
          *overlapfront_it = YGridComponent<Coordinates>(origin, r, &g.coords, size, n, size);

          // build overlap
          for (int i=0; i<dim; i++)
          {
            if (!r[i])
            {
              if (ovlp_low[i])
              {
                origin[i]++;
                size[i]--;
              }
              if (ovlp_up[i])
                size[i]--;
            }
          }
          *overlap_it = YGridComponent<Coordinates>(origin,size,*overlapfront_it);

          // build interiorborder
          for (int i=0; i<dim; i++)
          {
            if (ovlp_low[i])
            {
              origin[i] += overlap;
              size[i] -= overlap;
              if (!r[i])
              {
                origin[i]--;
                size[i]++;
              }
            }
            if (ovlp_up[i])
            {
              size[i] -= overlap;
              if (!r[i])
                size[i]++;
            }
          }
          *interiorborder_it = YGridComponent<Coordinates>(origin,size,*overlapfront_it);

          // build interior
          for (int i=0; i<dim; i++)
          {
            if (!r[i])
            {
              if (ovlp_low[i])
              {
                origin[i]++;
                size[i]--;
              }
              if (ovlp_up[i])
                size[i]--;
            }
          }
          *interior_it = YGridComponent<Coordinates>(origin, size, *overlapfront_it);

          intersections(*overlapfront_it,*overlapfront_it,*send_overlapfront_overlapfront_it, *recv_overlapfront_overlapfront_it);
          intersections(*overlap_it,*overlapfront_it,*send_overlap_overlapfront_it, *recv_overlapfront_overlap_it);
          intersections(*interiorborder_it,*interiorborder_it,*send_interiorborder_interiorborder_it,*recv_interiorborder_interiorborder_it);
          intersections(*interiorborder_it,*overlapfront_it,*send_interiorborder_overlapfront_it,*recv_overlapfront_interiorborder_it);

          // advance all iterators pointing to the next insertion point
          ++overlapfront_it;
          ++overlap_it;
          ++interiorborder_it;
          ++interior_it;
          ++send_overlapfront_overlapfront_it;
          ++recv_overlapfront_overlapfront_it;
          ++send_overlap_overlapfront_it;
          ++recv_overlapfront_overlap_it;
          ++send_interiorborder_interiorborder_it;
          ++recv_interiorborder_interiorborder_it;
          ++send_interiorborder_overlapfront_it;
          ++recv_overlapfront_interiorborder_it;
        }

        // set end iterators in the corresonding ygrids
        g.overlapfront[codim].finalize(overlapfront_it);
        g.overlap[codim].finalize(overlap_it);
        g.interiorborder[codim].finalize(interiorborder_it);
        g.interior[codim].finalize(interior_it);
        g.send_overlapfront_overlapfront[codim].finalize(send_overlapfront_overlapfront_it,g.overlapfront[codim]);
        g.recv_overlapfront_overlapfront[codim].finalize(recv_overlapfront_overlapfront_it,g.overlapfront[codim]);
        g.send_overlap_overlapfront[codim].finalize(send_overlap_overlapfront_it,g.overlapfront[codim]);
        g.recv_overlapfront_overlap[codim].finalize(recv_overlapfront_overlap_it,g.overlapfront[codim]);
        g.send_interiorborder_interiorborder[codim].finalize(send_interiorborder_interiorborder_it,g.overlapfront[codim]);
        g.recv_interiorborder_interiorborder[codim].finalize(recv_interiorborder_interiorborder_it,g.overlapfront[codim]);
        g.send_interiorborder_overlapfront[codim].finalize(send_interiorborder_overlapfront_it,g.overlapfront[codim]);
        g.recv_overlapfront_interiorborder[codim].finalize(recv_overlapfront_interiorborder_it,g.overlapfront[codim]);
      }
    }

#ifndef DOXYGEN
    /** \brief special data structure to communicate ygrids
     * Historically, this was needed because Ygrids had virtual functions and
     * a communicated virtual function table pointer introduced a bug. After the
     * change to tensorproductgrid, the dynamic polymorphism was removed, still this
     * is kept because it allows to communicate ygrids, that only have index, but no
     * coordinate information. This is sufficient, because all communicated YGrids are
     * intersected with a local grid, which has coordinate information.
     */
    struct mpifriendly_ygrid {
      mpifriendly_ygrid ()
      {
        std::fill(origin.begin(), origin.end(), 0);
        std::fill(size.begin(), size.end(), 0);
      }
      mpifriendly_ygrid (const YGridComponent<Coordinates>& grid)
        : origin(grid.origin()), size(grid.size())
      {}
      iTupel origin;
      iTupel size;
    };
#endif

    /** \brief Construct list of intersections with neighboring processors
     *
     * \param recvgrid the grid stored in this processor
     * \param sendgrid the subgrid to be sent to neighboring processors
     * \param sendlist the deque to fill with send intersections
     * \param recvlist the deque to fill with recv intersections
     * \returns two lists: Intersections to be sent and Intersections to be received
     */
    void intersections(const YGridComponent<Coordinates>& sendgrid, const YGridComponent<Coordinates>& recvgrid,
                        std::deque<Intersection>& sendlist, std::deque<Intersection>& recvlist)
    {
      iTupel size = globalSize();

      // the exchange buffers
      std::vector<YGridComponent<Coordinates> > send_recvgrid(_torus.neighbors());
      std::vector<YGridComponent<Coordinates> > recv_recvgrid(_torus.neighbors());
      std::vector<YGridComponent<Coordinates> > send_sendgrid(_torus.neighbors());
      std::vector<YGridComponent<Coordinates> > recv_sendgrid(_torus.neighbors());

      // new exchange buffers to send simple struct without virtual functions
      std::vector<mpifriendly_ygrid> mpifriendly_send_recvgrid(_torus.neighbors());
      std::vector<mpifriendly_ygrid> mpifriendly_recv_recvgrid(_torus.neighbors());
      std::vector<mpifriendly_ygrid> mpifriendly_send_sendgrid(_torus.neighbors());
      std::vector<mpifriendly_ygrid> mpifriendly_recv_sendgrid(_torus.neighbors());

      // fill send buffers; iterate over all neighboring processes
      // non-periodic case is handled automatically because intersection will be zero
      for (typename Torus<Communication,dim>::ProcListIterator i=_torus.sendbegin(); i!=_torus.sendend(); ++i)
      {
        // determine if we communicate with this neighbor (and what)
        bool skip = false;
        iTupel coord = _torus.coord();   // my coordinates
        iTupel delta = i.delta();        // delta to neighbor
        iTupel nb = coord;               // the neighbor
        for (int k=0; k<dim; k++) nb[k] += delta[k];
        iTupel v;                    // grid movement
        std::fill(v.begin(), v.end(), 0);

        for (int k=0; k<dim; k++)
        {
          if (nb[k]<0)
          {
            if (_periodic[k])
              v[k] += size[k];
            else
              skip = true;
          }
          if (nb[k]>=_torus.dims(k))
          {
            if (_periodic[k])
              v[k] -= size[k];
            else
              skip = true;
          }
          // neither might be true, then v=0
        }

        // store moved grids in send buffers
        if (!skip)
        {
          send_sendgrid[i.index()] = sendgrid.move(v);
          send_recvgrid[i.index()] = recvgrid.move(v);
        }
        else
        {
          send_sendgrid[i.index()] = YGridComponent<Coordinates>();
          send_recvgrid[i.index()] = YGridComponent<Coordinates>();
        }
      }

      // issue send requests for sendgrid being sent to all neighbors
      for (typename Torus<Communication,dim>::ProcListIterator i=_torus.sendbegin(); i!=_torus.sendend(); ++i)
      {
        mpifriendly_send_sendgrid[i.index()] = mpifriendly_ygrid(send_sendgrid[i.index()]);
        _torus.send(i.rank(), &mpifriendly_send_sendgrid[i.index()], sizeof(mpifriendly_ygrid));
      }

      // issue recv requests for sendgrids of neighbors
      for (typename Torus<Communication,dim>::ProcListIterator i=_torus.recvbegin(); i!=_torus.recvend(); ++i)
        _torus.recv(i.rank(), &mpifriendly_recv_sendgrid[i.index()], sizeof(mpifriendly_ygrid));

      // exchange the sendgrids
      _torus.exchange();

      // issue send requests for recvgrid being sent to all neighbors
      for (typename Torus<Communication,dim>::ProcListIterator i=_torus.sendbegin(); i!=_torus.sendend(); ++i)
      {
        mpifriendly_send_recvgrid[i.index()] = mpifriendly_ygrid(send_recvgrid[i.index()]);
        _torus.send(i.rank(), &mpifriendly_send_recvgrid[i.index()], sizeof(mpifriendly_ygrid));
      }

      // issue recv requests for recvgrid of neighbors
      for (typename Torus<Communication,dim>::ProcListIterator i=_torus.recvbegin(); i!=_torus.recvend(); ++i)
        _torus.recv(i.rank(), &mpifriendly_recv_recvgrid[i.index()], sizeof(mpifriendly_ygrid));

      // exchange the recvgrid
      _torus.exchange();

      // process receive buffers and compute intersections
      for (typename Torus<Communication,dim>::ProcListIterator i=_torus.recvbegin(); i!=_torus.recvend(); ++i)
      {
        // what must be sent to this neighbor
        Intersection send_intersection;
        mpifriendly_ygrid yg = mpifriendly_recv_recvgrid[i.index()];
        recv_recvgrid[i.index()] = YGridComponent<Coordinates>(yg.origin,yg.size);
        send_intersection.grid = sendgrid.intersection(recv_recvgrid[i.index()]);
        send_intersection.rank = i.rank();
        send_intersection.distance = i.distance();
        if (!send_intersection.grid.empty()) sendlist.push_front(send_intersection);

        Intersection recv_intersection;
        yg = mpifriendly_recv_sendgrid[i.index()];
        recv_sendgrid[i.index()] = YGridComponent<Coordinates>(yg.origin,yg.size);
        recv_intersection.grid = recvgrid.intersection(recv_sendgrid[i.index()]);
        recv_intersection.rank = i.rank();
        recv_intersection.distance = i.distance();
        if(!recv_intersection.grid.empty()) recvlist.push_back(recv_intersection);
      }
    }

  protected:

    typedef const YaspGrid<dim,Coordinates> GridImp;

    void init()
    {
      indexsets.push_back( std::make_shared< YaspIndexSet<const YaspGrid<dim, Coordinates>, false > >(*this,0) );
      boundarysegmentssize();
    }

    void boundarysegmentssize()
    {
      // sizes of local macro grid
      std::array<int, dim> sides;
      {
        for (int i=0; i<dim; i++)
        {
          sides[i] =
            ((begin()->overlap[0].dataBegin()->origin(i) == 0)+
             (begin()->overlap[0].dataBegin()->origin(i) + begin()->overlap[0].dataBegin()->size(i)
                    == levelSize(0,i)));
        }
      }
      nBSegments = 0;
      for (int k=0; k<dim; k++)
      {
        int offset = 1;
        for (int l=0; l<dim; l++)
        {
          if (l==k) continue;
          offset *= begin()->overlap[0].dataBegin()->size(l);
        }
        nBSegments += sides[k]*offset;
      }
    }

  public:

    // define the persistent index type
    typedef bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+dim> PersistentIndexType;

    //! the GridFamily of this grid
    typedef YaspGridFamily<dim, Coordinates> GridFamily;
    // the Traits
    typedef typename YaspGridFamily<dim, Coordinates>::Traits Traits;

    // need for friend declarations in entity
    typedef YaspIndexSet<YaspGrid<dim, Coordinates>, false > LevelIndexSetType;
    typedef YaspIndexSet<YaspGrid<dim, Coordinates>, true > LeafIndexSetType;
    typedef YaspGlobalIdSet<YaspGrid<dim, Coordinates> > GlobalIdSetType;

    /** Standard constructor for a YaspGrid with a given Coordinates object
     *  @param coordinates Object that stores or computes the vertex coordinates
     *  @param periodic tells if direction is periodic or not
     *  @param overlap size of overlap on coarsest grid (same in all directions)
     *  @param comm the communication object for this grid. An MPI communicator can be given here.
     *  @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (const Coordinates& coordinates,
              std::bitset<dim> periodic = std::bitset<dim>(0ULL),
              int overlap = 1,
              Communication comm = Communication(),
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
      : ccobj(comm)
      , leafIndexSet_(*this)
      , _periodic(periodic)
      , _overlap(overlap)
      , keep_ovlp(true)
      , adaptRefCount(0)
      , adaptActive(false)
    {
      _levels.resize(1);

      // Number of elements per coordinate direction on the coarsest level
      for (std::size_t i=0; i<dim; i++)
        _coarseSize[i] = coordinates.size(i);

      // Construct the communication torus
      _torus = decltype(_torus)(comm,tag,_coarseSize,lb);

      iTupel o;
      std::fill(o.begin(), o.end(), 0);
      iTupel o_interior(o);
      iTupel s_interior(_coarseSize);

      _torus.partition(_torus.rank(),o,_coarseSize,o_interior,s_interior);

      // Set domain size
      if (std::is_same<Coordinates,EquidistantCoordinates<ctype,dim> >::value
        || std::is_same<Coordinates,EquidistantOffsetCoordinates<ctype,dim> >::value)
      {
        for (std::size_t i=0; i<dim; i++)
          _L[i] = coordinates.size(i) * coordinates.meshsize(i,0);
      }
      if (std::is_same<Coordinates,TensorProductCoordinates<ctype,dim> >::value)
      {
        //determine sizes of vector to correctly construct torus structure and store for later size requests
        for (int i=0; i<dim; i++)
          _L[i] = coordinates.coordinate(i,_coarseSize[i]) - coordinates.coordinate(i,0);
      }

#if HAVE_MPI
      // TODO: Settle on a single value for all coordinate types
      int mysteryFactor = (std::is_same<Coordinates,EquidistantCoordinates<ctype,dim> >::value) ? 1 : 2;

      // check whether the grid is large enough to be overlapping
      for (int i=0; i<dim; i++)
      {
        // find out whether the grid is too small to
        int toosmall = (s_interior[i] <= mysteryFactor * overlap) &&    // interior is very small
            (periodic[i] || (s_interior[i] != _coarseSize[i]));    // there is an overlap in that direction
        // communicate the result to all those processes to have all processors error out if one process failed.
        int global = 0;
        MPI_Allreduce(&toosmall, &global, 1, MPI_INT, MPI_LOR, comm);
        if (global)
          DUNE_THROW(Dune::GridError,"YaspGrid does not support degrees of freedom shared by more than immediately neighboring subdomains."
                                     " Note that this also holds for DOFs on subdomain boundaries."
                                     " Increase grid elements or decrease overlap accordingly.");
      }
#endif // #if HAVE_MPI

      if (std::is_same<Coordinates,EquidistantCoordinates<ctype,dim> >::value
        || std::is_same<Coordinates,EquidistantOffsetCoordinates<ctype,dim> >::value)
      {
        iTupel s_overlap(s_interior);
        for (int i=0; i<dim; i++)
        {
          if ((o_interior[i] - overlap > 0) || (periodic[i]))
            s_overlap[i] += overlap;
          if ((o_interior[i] + s_interior[i] + overlap <= _coarseSize[i]) || (periodic[i]))
            s_overlap[i] += overlap;
        }

        FieldVector<ctype,dim> upperRightWithOverlap;
        for (int i=0; i<dim; i++)
          upperRightWithOverlap[i] = coordinates.coordinate(i,0) + coordinates.meshsize(i,0) * s_overlap[i];

        if constexpr (std::is_same_v<Coordinates,EquidistantCoordinates<ctype,dim>>)
        {
          // New coordinate object that additionally contains the overlap elements
          EquidistantCoordinates<ctype,dim> coordinatesWithOverlap(upperRightWithOverlap,s_overlap);

          // add level (the this-> is needed to make g++-6 happy)
          this->makelevel(coordinatesWithOverlap,periodic,o_interior,overlap);
        }

        if constexpr (std::is_same_v<Coordinates,EquidistantOffsetCoordinates<ctype,dim>>)
        {
          Dune::FieldVector<ctype,dim> lowerleft;
          for (int i=0; i<dim; i++)
            lowerleft[i] = coordinates.origin(i);

          // New coordinate object that additionally contains the overlap elements
          EquidistantOffsetCoordinates<ctype,dim> coordinatesWithOverlap(lowerleft,upperRightWithOverlap,s_overlap);

          // add level (the this-> is needed to make g++-6 happy)
          this->makelevel(coordinatesWithOverlap,periodic,o_interior,overlap);
        }
      }

      if constexpr (std::is_same_v<Coordinates,TensorProductCoordinates<ctype,dim>>)
      {
        std::array<std::vector<ctype>,dim> newCoords;
        std::array<int, dim> offset(o_interior);

        // find the relevant part of the coords vector for this processor and copy it to newCoords
        for (int i=0; i<dim; ++i)
        {
          //define the coordinate range to be used
          std::size_t begin = o_interior[i];
          std::size_t end   = begin + s_interior[i] + 1;

          // check whether we are not at the physical boundary. In that case overlap is a simple
          // extension of the coordinate range to be used
          if (o_interior[i] - overlap > 0)
          {
            begin = begin - overlap;
            offset[i] -= overlap;
          }
          if (o_interior[i] + s_interior[i] + overlap < _coarseSize[i])
            end = end + overlap;

          //copy the selected part in the new coord vector
          newCoords[i].resize(end-begin);
          auto newCoordsIt = newCoords[i].begin();
          for (std::size_t j=begin; j<end; j++)
          {
            *newCoordsIt = coordinates.coordinate(i, j);
            newCoordsIt++;
          }

          // Check whether we are at the physical boundary and have a periodic grid.
          // In this case the coordinate vector has to be tweaked manually.
          if ((periodic[i]) && (o_interior[i] + s_interior[i] + overlap >= _coarseSize[i]))
          {
            // we need to add the first <overlap> cells to the end of newcoords
            for (int j=0; j<overlap; ++j)
              newCoords[i].push_back(newCoords[i].back() - coordinates.coordinate(i,j) + coordinates.coordinate(i,j+1));
          }

          if ((periodic[i]) && (o_interior[i] - overlap <= 0))
          {
            offset[i] -= overlap;

            // we need to add the last <overlap> cells to the begin of newcoords
            std::size_t reverseCounter = coordinates.size(i);
            for (int j=0; j<overlap; ++j, --reverseCounter)
              newCoords[i].insert(newCoords[i].begin(), newCoords[i].front()
                                  - coordinates.coordinate(i,reverseCounter) + coordinates.coordinate(i,reverseCounter-1));
          }
        }

        TensorProductCoordinates<ctype,dim> coordinatesWithOverlap(newCoords, offset);

        // add level (the this-> is needed to make g++-6 happy)
        this->makelevel(coordinatesWithOverlap,periodic,o_interior,overlap);
      }

      init();
    }

    /** Standard constructor for an equidistant YaspGrid
     *  @param L extension of the domain
     *  @param s number of cells on coarse mesh in each direction
     *  @param periodic tells if direction is periodic or not
     *  @param overlap size of overlap on coarsest grid (same in all directions)
     *  @param comm the communication object for this grid. An MPI communicator can be given here.
     *  @param lb pointer to an overloaded YLoadBalance instance
     */
    template<class C = Coordinates,
             typename std::enable_if_t< std::is_same_v<C, EquidistantCoordinates<ctype,dim> >, int> = 0>
    YaspGrid (Dune::FieldVector<ctype, dim> L,
              std::array<int, std::size_t{dim}> s,
              std::bitset<std::size_t{dim}> periodic = std::bitset<std::size_t{dim}>{0ULL},
              int overlap = 1,
              Communication comm = Communication(),
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
      : ccobj(comm), _torus(comm,tag,s,lb), leafIndexSet_(*this),
        _L(L), _periodic(periodic), _coarseSize(s), _overlap(overlap),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
    {
      _levels.resize(1);

      iTupel o;
      std::fill(o.begin(), o.end(), 0);
      iTupel o_interior(o);
      iTupel s_interior(s);

      _torus.partition(_torus.rank(),o,s,o_interior,s_interior);

#if HAVE_MPI
      // check whether the grid is large enough to be overlapping
      for (int i=0; i<dim; i++)
      {
        // find out whether the grid is too small to
        int toosmall = (s_interior[i] / 2 <= overlap) &&    // interior is very small
            (periodic[i] || (s_interior[i] != s[i]));    // there is an overlap in that direction
        // communicate the result to all those processes to have all processors error out if one process failed.
        int global = 0;
        MPI_Allreduce(&toosmall, &global, 1, MPI_INT, MPI_LOR, comm);
        if (global)
          DUNE_THROW(Dune::GridError,"YaspGrid does not support degrees of freedom shared by more than immediately neighboring subdomains."
                                     " Note that this also holds for DOFs on subdomain boundaries."
                                     " Increase grid elements or decrease overlap accordingly.");
      }
#endif // #if HAVE_MPI

      iTupel s_overlap(s_interior);
      for (int i=0; i<dim; i++)
      {
        if ((o_interior[i] - overlap > 0) || (periodic[i]))
          s_overlap[i] += overlap;
        if ((o_interior[i] + s_interior[i] + overlap <= _coarseSize[i]) || (periodic[i]))
          s_overlap[i] += overlap;
      }

      FieldVector<ctype,dim> upperRightWithOverlap;

      for (int i=0; i<dim; i++)
        upperRightWithOverlap[i] = (L[i] / s[i]) * s_overlap[i];

      // New coordinate object that additionally contains the overlap elements
      EquidistantCoordinates<ctype,dim> cc(upperRightWithOverlap,s_overlap);

      // add level
      makelevel(cc,periodic,o_interior,overlap);

      init();
    }

    /** Constructor for an equidistant YaspGrid with non-trivial origin
     *  @param lowerleft Lower left corner of the domain
     *  @param upperright Upper right corner of the domain
     *  @param s number of cells on coarse mesh in each direction
     *  @param periodic tells if direction is periodic or not
     *  @param overlap size of overlap on coarsest grid (same in all directions)
     *  @param comm the communication object for this grid. An MPI communicator can be given here.
     *  @param lb pointer to an overloaded YLoadBalance instance
     */
    template<class C = Coordinates,
             typename std::enable_if_t< std::is_same_v<C, EquidistantOffsetCoordinates<ctype,dim> >, int> = 0>
    YaspGrid (Dune::FieldVector<ctype, dim> lowerleft,
              Dune::FieldVector<ctype, dim> upperright,
              std::array<int, std::size_t{dim}> s,
              std::bitset<std::size_t{dim}> periodic = std::bitset<std::size_t{dim}>(0ULL),
              int overlap = 1,
              Communication comm = Communication(),
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
      : ccobj(comm), _torus(comm,tag,s,lb), leafIndexSet_(*this),
        _L(upperright - lowerleft),
        _periodic(periodic), _coarseSize(s), _overlap(overlap),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
    {
      _levels.resize(1);

      iTupel o;
      std::fill(o.begin(), o.end(), 0);
      iTupel o_interior(o);
      iTupel s_interior(s);

      _torus.partition(_torus.rank(),o,s,o_interior,s_interior);

#if HAVE_MPI
      // check whether the grid is large enough to be overlapping
      for (int i=0; i<dim; i++)
      {
        // find out whether the grid is too small to
        int toosmall = (s_interior[i] / 2 <= overlap) &&    // interior is very small
            (periodic[i] || (s_interior[i] != s[i]));    // there is an overlap in that direction
        // communicate the result to all those processes to have all processors error out if one process failed.
        int global = 0;
        MPI_Allreduce(&toosmall, &global, 1, MPI_INT, MPI_LOR, comm);
        if (global)
          DUNE_THROW(Dune::GridError,"YaspGrid does not support degrees of freedom shared by more than immediately neighboring subdomains."
                                     " Note that this also holds for DOFs on subdomain boundaries."
                                     " Increase grid elements or decrease overlap accordingly.");
      }
#endif // #if HAVE_MPI

      iTupel s_overlap(s_interior);
      for (int i=0; i<dim; i++)
      {
        if ((o_interior[i] - overlap > 0) || (periodic[i]))
          s_overlap[i] += overlap;
        if ((o_interior[i] + s_interior[i] + overlap <= _coarseSize[i]) || (periodic[i]))
          s_overlap[i] += overlap;
      }

      FieldVector<ctype,dim> upperRightWithOverlap;
      for (int i=0; i<dim; i++)
        upperRightWithOverlap[i] = lowerleft[i]
                                 + s_overlap[i] * (upperright[i]-lowerleft[i]) / s[i];

      EquidistantOffsetCoordinates<ctype,dim> cc(lowerleft,upperRightWithOverlap,s_overlap);

      // add level
      makelevel(cc,periodic,o_interior,overlap);

      init();
    }

    /** @brief Standard constructor for a tensorproduct YaspGrid
     *  @param coords coordinate vectors to be used for coarse grid
     *  @param periodic tells if direction is periodic or not
     *  @param overlap size of overlap on coarsest grid (same in all directions)
     *  @param comm the communication object for this grid. An MPI communicator can be given here.
     *  @param lb pointer to an overloaded YLoadBalance instance
     */
    template<class C = Coordinates,
             typename std::enable_if_t< std::is_same_v<C, TensorProductCoordinates<ctype,dim> >, int> = 0>
    YaspGrid (std::array<std::vector<ctype>, std::size_t{dim}> coords,
              std::bitset<std::size_t{dim}> periodic = std::bitset<std::size_t{dim}>(0ULL),
              int overlap = 1,
              Communication comm = Communication(),
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
      : ccobj(comm), _torus(comm,tag,Dune::Yasp::sizeArray<dim>(coords),lb),
        leafIndexSet_(*this), _periodic(periodic), _overlap(overlap),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
    {
      if (!Dune::Yasp::checkIfMonotonous(coords))
        DUNE_THROW(Dune::GridError,"Setup of a tensorproduct grid requires monotonous sequences of coordinates.");

      _levels.resize(1);

      //determine sizes of vector to correctly construct torus structure and store for later size requests
      for (int i=0; i<dim; i++) {
        _coarseSize[i] = coords[i].size() - 1;
        _L[i] = coords[i][_coarseSize[i]] - coords[i][0];
      }

      iTupel o;
      std::fill(o.begin(), o.end(), 0);
      iTupel o_interior(o);
      iTupel s_interior(_coarseSize);

      _torus.partition(_torus.rank(),o,_coarseSize,o_interior,s_interior);

#if HAVE_MPI
      // check whether the grid is large enough to be overlapping
      for (int i=0; i<dim; i++)
      {
        // find out whether the grid is too small to
        int toosmall = (s_interior[i] / 2 <= overlap) &&               // interior is very small
             (periodic[i] || (s_interior[i] != _coarseSize[i]));    // there is an overlap in that direction
        // communicate the result to all those processes to have all processors error out if one process failed.
        int global = 0;
        MPI_Allreduce(&toosmall, &global, 1, MPI_INT, MPI_LOR, comm);
        if (global)
          DUNE_THROW(Dune::GridError,"YaspGrid does not support degrees of freedom shared by more than immediately neighboring subdomains."
                                     " Note that this also holds for DOFs on subdomain boundaries."
                                     " Increase grid elements or decrease overlap accordingly.");
      }
#endif // #if HAVE_MPI


      std::array<std::vector<ctype>,dim> newcoords;
      std::array<int, dim> offset(o_interior);

      // find the relevant part of the coords vector for this processor and copy it to newcoords
      for (int i=0; i<dim; ++i)
      {
        //define iterators on coords that specify the coordinate range to be used
        typename std::vector<ctype>::iterator begin = coords[i].begin() + o_interior[i];
        typename std::vector<ctype>::iterator end = begin + s_interior[i] + 1;

        // check whether we are not at the physical boundary. In that case overlap is a simple
        // extension of the coordinate range to be used
        if (o_interior[i] - overlap > 0)
        {
          begin = begin - overlap;
          offset[i] -= overlap;
        }
        if (o_interior[i] + s_interior[i] + overlap < _coarseSize[i])
          end = end + overlap;

        //copy the selected part in the new coord vector
        newcoords[i].resize(end-begin);
        std::copy(begin, end, newcoords[i].begin());

        // check whether we are at the physical boundary and a have a periodic grid.
        // In this case the coordinate vector has to be tweaked manually.
        if ((periodic[i]) && (o_interior[i] + s_interior[i] + overlap >= _coarseSize[i]))
        {
          // we need to add the first <overlap> cells to the end of newcoords
          typename std::vector<ctype>::iterator it = coords[i].begin();
          for (int j=0; j<overlap; ++j)
            newcoords[i].push_back(newcoords[i].back() - *it + *(++it));
        }

        if ((periodic[i]) && (o_interior[i] - overlap <= 0))
        {
          offset[i] -= overlap;

          // we need to add the last <overlap> cells to the begin of newcoords
          typename std::vector<ctype>::iterator it = coords[i].end() - 1;
          for (int j=0; j<overlap; ++j)
            newcoords[i].insert(newcoords[i].begin(), newcoords[i].front() - *it + *(--it));
        }
      }

      TensorProductCoordinates<ctype,dim> cc(newcoords, offset);

      // add level
      makelevel(cc,periodic,o_interior,overlap);
      init();
    }

  private:

    /** @brief Constructor for a tensorproduct YaspGrid with only coordinate
     *         information on this processor
     *  @param comm MPI communicator where this mesh is distributed to
     *  @param coords coordinate vectors to be used for coarse grid
     *  @param periodic tells if direction is periodic or not
     *  @param overlap size of overlap on coarsest grid (same in all directions)
     *  @param coarseSize the coarse size of the global grid
     *  @param lb pointer to an overloaded YLoadBalance instance
     *
     *  @warning The construction of overlapping coordinate ranges is
     *           an error-prone procedure. For this reason, it is kept private.
     *           You can safely use it through BackupRestoreFacility. All other
     *           use is not supported for the moment.
     */
    YaspGrid (std::array<std::vector<ctype>, std::size_t{dim}> coords,
              std::bitset<std::size_t{dim}> periodic,
              int overlap,
              Communication comm,
              std::array<int,dim> coarseSize,
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
      : ccobj(comm), _torus(comm,tag,coarseSize,lb), leafIndexSet_(*this),
        _periodic(periodic), _coarseSize(coarseSize), _overlap(overlap),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
    {
      // check whether YaspGrid has been given the correct template parameter
      static_assert(std::is_same<Coordinates,TensorProductCoordinates<ctype,dim> >::value,
                  "YaspGrid coordinate container template parameter and given constructor values do not match!");

      if (!Dune::Yasp::checkIfMonotonous(coords))
        DUNE_THROW(Dune::GridError,"Setup of a tensorproduct grid requires monotonous sequences of coordinates.");

      for (int i=0; i<dim; i++)
        _L[i] = coords[i][coords[i].size() - 1] - coords[i][0];

      _levels.resize(1);

      std::array<int,dim> o;
      std::fill(o.begin(), o.end(), 0);
      std::array<int,dim> o_interior(o);
      std::array<int,dim> s_interior(coarseSize);

      _torus.partition(_torus.rank(),o,coarseSize,o_interior,s_interior);

      // get offset by modifying o_interior according to overlap
      std::array<int,dim> offset(o_interior);
      for (int i=0; i<dim; i++)
        if ((periodic[i]) || (o_interior[i] > 0))
          offset[i] -= overlap;

      TensorProductCoordinates<ctype,dim> cc(coords, offset);

      // add level
      makelevel(cc,periodic,o_interior,overlap);

      init();
    }

    // the backup restore facility needs to be able to use above constructor
    friend struct BackupRestoreFacility<YaspGrid<dim,Coordinates> >;

    // do not copy this class
    YaspGrid(const YaspGrid&);

  public:

    /*! Return maximum level defined in this grid. Levels are numbered
          0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const
    {
      return _levels.size()-1;
    }

    //! refine the grid refCount times.
    void globalRefine (int refCount)
    {
      if (refCount < -maxLevel())
        DUNE_THROW(GridError, "Only " << maxLevel() << " levels left. " <<
                   "Coarsening " << -refCount << " levels requested!");

      // If refCount is negative then coarsen the grid
      for (int k=refCount; k<0; k++)
      {
        // create an empty grid level
        YGridLevel empty;
        _levels.back() = empty;
        // reduce maxlevel
        _levels.pop_back();

        indexsets.pop_back();
      }

      // If refCount is positive refine the grid
      for (int k=0; k<refCount; k++)
      {
        // access to coarser grid level
        YGridLevel& cg = _levels[maxLevel()];

        std::bitset<dim> ovlp_low(0ULL), ovlp_up(0ULL);
        for (int i=0; i<dim; i++)
        {
          if (cg.overlap[0].dataBegin()->origin(i) > 0 || _periodic[i])
            ovlp_low[i] = true;
          if (cg.overlap[0].dataBegin()->max(i) + 1 < globalSize(i) || _periodic[i])
            ovlp_up[i] = true;
        }

        Coordinates newcont(cg.coords.refine(ovlp_low, ovlp_up, cg.overlapSize, keep_ovlp));

        int overlap = (keep_ovlp) ? 2*cg.overlapSize : cg.overlapSize;

        //determine new origin
        iTupel o_interior;
        for (int i=0; i<dim; i++)
          o_interior[i] = 2*cg.interior[0].dataBegin()->origin(i);

        // add level
        _levels.resize(_levels.size() + 1);
        makelevel(newcont,_periodic,o_interior,overlap);

        indexsets.push_back( std::make_shared<YaspIndexSet<const YaspGrid<dim,Coordinates>, false > >(*this,maxLevel()) );
      }
    }

    /**
       \brief set options for refinement
       @param keepPhysicalOverlap [true] keep the physical size of the overlap, [false] keep the number of cells in the overlap.  Default is [true].
     */
    void refineOptions (bool keepPhysicalOverlap)
    {
      keep_ovlp = keepPhysicalOverlap;
    }

    /** \brief Marks an entity to be refined/coarsened in a subsequent adapt.

       \param[in] refCount Number of subdivisions that should be applied. Negative value means coarsening.
       \param[in] e        Entity to Entity that should be refined

       \return true if Entity was marked, false otherwise.

       \note
          -  On yaspgrid marking one element will mark all other elements of the level as well
          -  If refCount is lower than refCount of a previous mark-call, nothing is changed
     */
    bool mark( int refCount, const typename Traits::template Codim<0>::Entity & e )
    {
      assert(adaptActive == false);
      if (e.level() != maxLevel()) return false;
      adaptRefCount = std::max(adaptRefCount, refCount);
      return true;
    }

    /** \brief returns adaptation mark for given entity

       \param[in] e   Entity for which adaptation mark should be determined

       \return int adaptation mark, here the default value 0 is returned
     */
    int getMark ( const typename Traits::template Codim<0>::Entity &e ) const
    {
      return ( e.level() == maxLevel() ) ? adaptRefCount : 0;
    }

    //! map adapt to global refine
    bool adapt ()
    {
      globalRefine(adaptRefCount);
      return (adaptRefCount > 0);
    }

    //! returns true, if the grid will be coarsened
    bool preAdapt ()
    {
      adaptActive = true;
      adaptRefCount = comm().max(adaptRefCount);
      return (adaptRefCount < 0);
    }

    //! clean up some markers
    void postAdapt()
    {
      adaptActive = false;
      adaptRefCount = 0;
    }

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lbegin (int level) const
    {
      return levelbegin<cd,pitype>(level);
    }

    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lend (int level) const
    {
      return levelend<cd,pitype>(level);
    }

    //! version without second template parameter for convenience
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lbegin (int level) const
    {
      return levelbegin<cd,All_Partition>(level);
    }

    //! version without second template parameter for convenience
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lend (int level) const
    {
      return levelend<cd,All_Partition>(level);
    }

    //! return LeafIterator which points to the first entity in maxLevel
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafbegin () const
    {
      return levelbegin<cd,pitype>(maxLevel());
    }

    //! return LeafIterator which points behind the last entity in maxLevel
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafend () const
    {
      return levelend<cd,pitype>(maxLevel());
    }

    //! return LeafIterator which points to the first entity in maxLevel
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafbegin () const
    {
      return levelbegin<cd,All_Partition>(maxLevel());
    }

    //! return LeafIterator which points behind the last entity in maxLevel
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafend () const
    {
      return levelend<cd,All_Partition>(maxLevel());
    }

    // \brief obtain Entity from EntitySeed. */
    template <typename Seed>
    typename Traits::template Codim<Seed::codimension>::Entity
    entity(const Seed& seed) const
    {
      const int codim = Seed::codimension;
      YGridLevelIterator g = begin(seed.impl().level());

      typedef typename Traits::template Codim<Seed::codimension>::Entity Entity;
      typedef YaspEntity<codim,dim,const YaspGrid> EntityImp;
      typedef typename YGrid::Iterator YIterator;

      return Entity(EntityImp(g,YIterator(g->overlapfront[codim],seed.impl().coord(),seed.impl().offset())));
    }

    //! return size (= distance in graph) of overlap region
    int overlapSize (int level, [[maybe_unused]] int codim) const
    {
      YGridLevelIterator g = begin(level);
      return g->overlapSize;
    }

    //! return size (= distance in graph) of overlap region
    int overlapSize ([[maybe_unused]] int odim) const
    {
      YGridLevelIterator g = begin(maxLevel());
      return g->overlapSize;
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize ([[maybe_unused]] int level, [[maybe_unused]] int codim) const
    {
      return 0;
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize ([[maybe_unused]] int codim) const
    {
      return 0;
    }

    //! number of entities per level and codim in this process
    int size (int level, int codim) const
    {
      YGridLevelIterator g = begin(level);

      // sum over all components of the codimension
      int count = 0;
      typedef typename std::array<YGridComponent<Coordinates>, Dune::power(2,dim)>::iterator DAI;
      for (DAI it = g->overlapfront[codim].dataBegin(); it != g->overlapfront[codim].dataEnd(); ++it)
        count += it->totalsize();

      return count;
    }

    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return size(maxLevel(),codim);
    }

    //! number of entities per level and geometry type in this process
    int size (int level, GeometryType type) const
    {
      return (type.isCube()) ? size(level,dim-type.dim()) : 0;
    }

    //! number of leaf entities per geometry type in this process
    int size (GeometryType type) const
    {
      return size(maxLevel(),type);
    }

    //! \brief returns the number of boundary segments within the macro grid
    size_t numBoundarySegments () const
    {
      return nBSegments;
    }

    //! \brief returns the size of the physical domain
    const Dune::FieldVector<ctype, dim>& domainSize () const {
      return _L;
    }

    /*! The new communication interface

       communicate objects for all codims on a given level
     */
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp,DataType> & data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {
      YaspCommunicateMeta<dim,dim>::comm(*this,data,iftype,dir,level);
    }

    /*! The new communication interface

       communicate objects for all codims on the leaf grid
     */
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp,DataType> & data, InterfaceType iftype, CommunicationDirection dir) const
    {
      YaspCommunicateMeta<dim,dim>::comm(*this,data,iftype,dir,this->maxLevel());
    }

    /*! The new communication interface

       communicate objects for one codim
     */
    template<class DataHandle, int codim>
    void communicateCodim (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {
      // check input
      if (!data.contains(dim,codim)) return; // should have been checked outside

      // data types
      typedef typename DataHandle::DataType DataType;

      // access to grid level
      YGridLevelIterator g = begin(level);

      // find send/recv lists or throw error
      const YGridList<Coordinates>* sendlist = 0;
      const YGridList<Coordinates>* recvlist = 0;

      if (iftype==InteriorBorder_InteriorBorder_Interface)
      {
        sendlist = &g->send_interiorborder_interiorborder[codim];
        recvlist = &g->recv_interiorborder_interiorborder[codim];
      }
      if (iftype==InteriorBorder_All_Interface)
      {
        sendlist = &g->send_interiorborder_overlapfront[codim];
        recvlist = &g->recv_overlapfront_interiorborder[codim];
      }
      if (iftype==Overlap_OverlapFront_Interface || iftype==Overlap_All_Interface)
      {
        sendlist = &g->send_overlap_overlapfront[codim];
        recvlist = &g->recv_overlapfront_overlap[codim];
      }
      if (iftype==All_All_Interface)
      {
        sendlist = &g->send_overlapfront_overlapfront[codim];
        recvlist = &g->recv_overlapfront_overlapfront[codim];
      }

      // change communication direction?
      if (dir==BackwardCommunication)
        std::swap(sendlist,recvlist);

      int cnt;

      // Size computation (requires communication if variable size)
      std::vector<int> send_size(sendlist->size(),-1);    // map rank to total number of objects (of type DataType) to be sent
      std::vector<int> recv_size(recvlist->size(),-1);    // map rank to total number of objects (of type DataType) to be recvd
      std::vector<size_t*> send_sizes(sendlist->size(),static_cast<size_t*>(0)); // map rank to array giving number of objects per entity to be sent
      std::vector<size_t*> recv_sizes(recvlist->size(),static_cast<size_t*>(0)); // map rank to array giving number of objects per entity to be recvd

      // define type to iterate over send and recv lists
      typedef typename YGridList<Coordinates>::Iterator ListIt;

      if (data.fixedSize(dim,codim))
      {
        // fixed size: just take a dummy entity, size can be computed without communication
        cnt=0;
        for (ListIt is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(g, typename YGrid::Iterator(is->yg)));
          send_size[cnt] = is->grid.totalsize() * data.size(*it);
          cnt++;
        }
        cnt=0;
        for (ListIt is=recvlist->begin(); is!=recvlist->end(); ++is)
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(g, typename YGrid::Iterator(is->yg)));
          recv_size[cnt] = is->grid.totalsize() * data.size(*it);
          cnt++;
        }
      }
      else
      {
        // variable size case: sender side determines the size
        cnt=0;
        for (ListIt is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          // allocate send buffer for sizes per entitiy
          size_t *buf = new size_t[is->grid.totalsize()];
          send_sizes[cnt] = buf;

          // loop over entities and ask for size
          int i=0; size_t n=0;
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(g, typename YGrid::Iterator(is->yg)));
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          itend(YaspLevelIterator<codim,All_Partition,GridImp>(g, typename YGrid::Iterator(is->yg,true)));
          for ( ; it!=itend; ++it)
          {
            buf[i] = data.size(*it);
            n += buf[i];
            i++;
          }

          // now we know the size for this rank
          send_size[cnt] = n;

          // hand over send request to torus class
          torus().send(is->rank,buf,is->grid.totalsize()*sizeof(size_t));
          cnt++;
        }

        // allocate recv buffers for sizes and store receive request
        cnt=0;
        for (ListIt is=recvlist->begin(); is!=recvlist->end(); ++is)
        {
          // allocate recv buffer
          size_t *buf = new size_t[is->grid.totalsize()];
          recv_sizes[cnt] = buf;

          // hand over recv request to torus class
          torus().recv(is->rank,buf,is->grid.totalsize()*sizeof(size_t));
          cnt++;
        }

        // exchange all size buffers now
        torus().exchange();

        // release send size buffers
        cnt=0;
        for (ListIt is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          delete[] send_sizes[cnt];
          send_sizes[cnt] = 0;
          cnt++;
        }

        // process receive size buffers
        cnt=0;
        for (ListIt is=recvlist->begin(); is!=recvlist->end(); ++is)
        {
          // get recv buffer
          size_t *buf = recv_sizes[cnt];

          // compute total size
          size_t n=0;
          for (int i=0; i<is->grid.totalsize(); ++i)
            n += buf[i];

          // ... and store it
          recv_size[cnt] = n;
          ++cnt;
        }
      }


      // allocate & fill the send buffers & store send request
      std::vector<DataType*> sends(sendlist->size(), static_cast<DataType*>(0)); // store pointers to send buffers
      cnt=0;
      for (ListIt is=sendlist->begin(); is!=sendlist->end(); ++is)
      {
        // allocate send buffer
        DataType *buf = new DataType[send_size[cnt]];

        // remember send buffer
        sends[cnt] = buf;

        // make a message buffer
        MessageBuffer<DataType> mb(buf);

        // fill send buffer; iterate over cells in intersection
        typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
        it(YaspLevelIterator<codim,All_Partition,GridImp>(g, typename YGrid::Iterator(is->yg)));
        typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
        itend(YaspLevelIterator<codim,All_Partition,GridImp>(g, typename YGrid::Iterator(is->yg,true)));
        for ( ; it!=itend; ++it)
          data.gather(mb,*it);

        // hand over send request to torus class
        torus().send(is->rank,buf,send_size[cnt]*sizeof(DataType));
        cnt++;
      }

      // allocate recv buffers and store receive request
      std::vector<DataType*> recvs(recvlist->size(),static_cast<DataType*>(0)); // store pointers to send buffers
      cnt=0;
      for (ListIt is=recvlist->begin(); is!=recvlist->end(); ++is)
      {
        // allocate recv buffer
        DataType *buf = new DataType[recv_size[cnt]];

        // remember recv buffer
        recvs[cnt] = buf;

        // hand over recv request to torus class
        torus().recv(is->rank,buf,recv_size[cnt]*sizeof(DataType));
        cnt++;
      }

      // exchange all buffers now
      torus().exchange();

      // release send buffers
      cnt=0;
      for (ListIt is=sendlist->begin(); is!=sendlist->end(); ++is)
      {
        delete[] sends[cnt];
        sends[cnt] = 0;
        cnt++;
      }

      // process receive buffers and delete them
      cnt=0;
      for (ListIt is=recvlist->begin(); is!=recvlist->end(); ++is)
      {
        // get recv buffer
        DataType *buf = recvs[cnt];

        // make a message buffer
        MessageBuffer<DataType> mb(buf);

        // copy data from receive buffer; iterate over cells in intersection
        if (data.fixedSize(dim,codim))
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(g, typename YGrid::Iterator(is->yg)));
          size_t n=data.size(*it);
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          itend(YaspLevelIterator<codim,All_Partition,GridImp>(g, typename YGrid::Iterator(is->yg,true)));
          for ( ; it!=itend; ++it)
            data.scatter(mb,*it,n);
        }
        else
        {
          int i=0;
          size_t *sbuf = recv_sizes[cnt];
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(g, typename YGrid::Iterator(is->yg)));
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          itend(YaspLevelIterator<codim,All_Partition,GridImp>(g, typename YGrid::Iterator(is->yg,true)));
          for ( ; it!=itend; ++it)
            data.scatter(mb,*it,sbuf[i++]);
          delete[] sbuf;
        }

        // delete buffer
        delete[] buf; // hier krachts !
        cnt++;
      }
    }

    // The new index sets from DDM 11.07.2005
    const typename Traits::GlobalIdSet& globalIdSet() const
    {
      return theglobalidset;
    }

    const typename Traits::LocalIdSet& localIdSet() const
    {
      return theglobalidset;
    }

    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");
      return *(indexsets[level]);
    }

    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return leafIndexSet_;
    }

    /*! @brief return a communication object
     */
    const Communication& comm () const
    {
      return ccobj;
    }

  private:

    // number of boundary segments of the level 0 grid
    int nBSegments;

    // Index classes need access to the real entity
    friend class Dune::YaspIndexSet<const Dune::YaspGrid<dim, Coordinates>, true >;
    friend class Dune::YaspIndexSet<const Dune::YaspGrid<dim, Coordinates>, false >;
    friend class Dune::YaspGlobalIdSet<const Dune::YaspGrid<dim, Coordinates> >;
    friend class Dune::YaspPersistentContainerIndex<const Dune::YaspGrid<dim, Coordinates> >;

    friend class Dune::YaspIntersectionIterator<const Dune::YaspGrid<dim, Coordinates> >;
    friend class Dune::YaspIntersection<const Dune::YaspGrid<dim, Coordinates> >;
    friend class Dune::YaspEntity<0, dim, const Dune::YaspGrid<dim, Coordinates> >;

    template<int codim_, int dim_, class GridImp_, template<int,int,class> class EntityImp_>
    friend class Entity;

    template<class DT>
    class MessageBuffer {
    public:
      // Constructor
      MessageBuffer (DT *p)
      {
        a=p;
        i=0;
        j=0;
      }

      // write data to message buffer, acts like a stream !
      template<class Y>
      void write (const Y& data)
      {
        static_assert(( std::is_same<DT,Y>::value ), "DataType mismatch");
        a[i++] = data;
      }

      // read data from message buffer, acts like a stream !
      template<class Y>
      void read (Y& data) const
      {
        static_assert(( std::is_same<DT,Y>::value ), "DataType mismatch");
        data = a[j++];
      }

    private:
      DT *a;
      int i;
      mutable int j;
    };

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    YaspLevelIterator<cd,pitype,GridImp> levelbegin (int level) const
    {
      YGridLevelIterator g = begin(level);
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");

      if (pitype==Interior_Partition)
        return YaspLevelIterator<cd,pitype,GridImp>(g,g->interior[cd].begin());
      if (pitype==InteriorBorder_Partition)
        return YaspLevelIterator<cd,pitype,GridImp>(g,g->interiorborder[cd].begin());
      if (pitype==Overlap_Partition)
        return YaspLevelIterator<cd,pitype,GridImp>(g,g->overlap[cd].begin());
      if (pitype<=All_Partition)
        return YaspLevelIterator<cd,pitype,GridImp>(g,g->overlapfront[cd].begin());
      if (pitype==Ghost_Partition)
        return levelend <cd, pitype> (level);

      DUNE_THROW(GridError, "YaspLevelIterator with this codim or partition type not implemented");
    }

    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    YaspLevelIterator<cd,pitype,GridImp> levelend (int level) const
    {
      YGridLevelIterator g = begin(level);
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");

      if (pitype==Interior_Partition)
        return YaspLevelIterator<cd,pitype,GridImp>(g,g->interior[cd].end());
      if (pitype==InteriorBorder_Partition)
        return YaspLevelIterator<cd,pitype,GridImp>(g,g->interiorborder[cd].end());
      if (pitype==Overlap_Partition)
        return YaspLevelIterator<cd,pitype,GridImp>(g,g->overlap[cd].end());
      if (pitype<=All_Partition || pitype == Ghost_Partition)
        return YaspLevelIterator<cd,pitype,GridImp>(g,g->overlapfront[cd].end());

      DUNE_THROW(GridError, "YaspLevelIterator with this codim or partition type not implemented");
    }

    Communication ccobj;

    Torus<Communication,dim> _torus;

    std::vector< std::shared_ptr< YaspIndexSet<const YaspGrid<dim,Coordinates>, false > > > indexsets;
    YaspIndexSet<const YaspGrid<dim,Coordinates>, true> leafIndexSet_;
    YaspGlobalIdSet<const YaspGrid<dim,Coordinates> > theglobalidset;

    Dune::FieldVector<ctype, dim> _L;
    iTupel _s;
    std::bitset<dim> _periodic;
    iTupel _coarseSize;
    ReservedVector<YGridLevel,32> _levels;
    int _overlap;
    bool keep_ovlp;
    int adaptRefCount;
    bool adaptActive;
  };

#if __cpp_deduction_guides >= 201611
  // Class template deduction guides
  template<typename ctype, int dim>
  YaspGrid(FieldVector<ctype, dim>,
           std::array<int, std::size_t{dim}>,
           std::bitset<std::size_t{dim}> = std::bitset<std::size_t{dim}>{0ULL},
           int = 1,
           YaspCommunication = YaspCommunication(),
           const YLoadBalance<dim>* = YaspGrid< dim, EquidistantCoordinates<ctype, dim> >::defaultLoadbalancer())
    -> YaspGrid< dim, EquidistantCoordinates<ctype, dim> >;

  template<typename ctype, int dim>
  YaspGrid(FieldVector<ctype, dim>,
           FieldVector<ctype, dim>,
           std::array<int, std::size_t{dim}>,
           std::bitset<std::size_t{dim}> = std::bitset<std::size_t{dim}>{0ULL},
           int = 1,
           YaspCommunication = YaspCommunication(),
           const YLoadBalance<dim>* = YaspGrid< dim, EquidistantOffsetCoordinates<ctype, dim> >::defaultLoadbalancer())
    -> YaspGrid< dim, EquidistantOffsetCoordinates<ctype, dim> >;

  template<typename ctype, std::size_t dim>
  YaspGrid(std::array<std::vector<ctype>, dim>,
           std::bitset<dim> = std::bitset<dim>{0ULL},
           int = 1,
           YaspCommunication = YaspCommunication(),
           const YLoadBalance<int{dim}>* = YaspGrid< int{dim}, TensorProductCoordinates<ctype, int{dim}> >::defaultLoadbalancer())
    -> YaspGrid< int{dim}, TensorProductCoordinates<ctype, int{dim}> >;
#endif

  //! Output operator for multigrids
  template <int d, class CC>
  std::ostream& operator<< (std::ostream& s, const YaspGrid<d,CC>& grid)
  {
    int rank = grid.torus().rank();

    s << "[" << rank << "]:" << " YaspGrid maxlevel=" << grid.maxLevel() << std::endl;

    s << "Printing the torus: " <<std::endl;
    s << grid.torus() << std::endl;

    for (typename YaspGrid<d,CC>::YGridLevelIterator g=grid.begin(); g!=grid.end(); ++g)
    {
      s << "[" << rank << "]:   " << std::endl;
      s << "[" << rank << "]:   " << "==========================================" << std::endl;
      s << "[" << rank << "]:   " << "level=" << g->level() << std::endl;

      for (int codim = 0; codim < d + 1; ++codim)
      {
        s << "[" << rank << "]:   " << "overlapfront[" << codim << "]:    " << g->overlapfront[codim] << std::endl;
        s << "[" << rank << "]:   " << "overlap[" << codim << "]:    " << g->overlap[codim] << std::endl;
        s << "[" << rank << "]:   " << "interiorborder[" << codim << "]:    " << g->interiorborder[codim] << std::endl;
        s << "[" << rank << "]:   " << "interior[" << codim << "]:    " << g->interior[codim] << std::endl;

        typedef typename YGridList<CC>::Iterator I;
        for (I i=g->send_overlapfront_overlapfront[codim].begin();
                 i!=g->send_overlapfront_overlapfront[codim].end(); ++i)
          s << "[" << rank << "]:    " << " s_of_of[" << codim << "] to rank "
                   << i->rank << " " << i->grid << std::endl;

        for (I i=g->recv_overlapfront_overlapfront[codim].begin();
                 i!=g->recv_overlapfront_overlapfront[codim].end(); ++i)
          s << "[" << rank << "]:    " << " r_of_of[" << codim << "] to rank "
                   << i->rank << " " << i->grid << std::endl;

        for (I i=g->send_overlap_overlapfront[codim].begin();
                 i!=g->send_overlap_overlapfront[codim].end(); ++i)
          s << "[" << rank << "]:    " << " s_o_of[" << codim << "] to rank "
                   << i->rank << " " << i->grid << std::endl;

        for (I i=g->recv_overlapfront_overlap[codim].begin();
                 i!=g->recv_overlapfront_overlap[codim].end(); ++i)
          s << "[" << rank << "]:    " << " r_of_o[" << codim << "] to rank "
                   << i->rank << " " << i->grid << std::endl;

        for (I i=g->send_interiorborder_interiorborder[codim].begin();
                 i!=g->send_interiorborder_interiorborder[codim].end(); ++i)
          s << "[" << rank << "]:    " << " s_ib_ib[" << codim << "] to rank "
          << i->rank << " " << i->grid << std::endl;

        for (I i=g->recv_interiorborder_interiorborder[codim].begin();
                 i!=g->recv_interiorborder_interiorborder[codim].end(); ++i)
             s << "[" << rank << "]:    " << " r_ib_ib[" << codim << "] to rank "
             << i->rank << " " << i->grid << std::endl;

        for (I i=g->send_interiorborder_overlapfront[codim].begin();
                 i!=g->send_interiorborder_overlapfront[codim].end(); ++i)
             s << "[" << rank << "]:    " << " s_ib_of[" << codim << "] to rank "
             << i->rank << " " << i->grid << std::endl;

        for (I i=g->recv_overlapfront_interiorborder[codim].begin();
                 i!=g->recv_overlapfront_interiorborder[codim].end(); ++i)
             s << "[" << rank << "]:    " << " r_of_ib[" << codim << "] to rank "
             << i->rank << " " << i->grid << std::endl;
      }
    }

    s << std::endl;

    return s;
  }

  namespace Capabilities
  {

    /** \struct hasEntity
       \ingroup YaspGrid
     */

    /** \struct hasBackupRestoreFacilities
       \ingroup YaspGrid
     */
    template<int dim, class Coordinates>
    struct hasBackupRestoreFacilities< YaspGrid<dim, Coordinates> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid has only one geometry type for codim 0 entities
       \ingroup YaspGrid
     */
    template<int dim, class Coordinates>
    struct hasSingleGeometryType< YaspGrid<dim, Coordinates> >
    {
      static const bool v = true;
      static const unsigned int topologyId = GeometryTypes::cube(dim).id();
    };

    /** \brief YaspGrid is a Cartesian grid
        \ingroup YaspGrid
     */
    template<int dim, class Coordinates>
    struct isCartesian< YaspGrid<dim, Coordinates> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid has entities for all codimensions
       \ingroup YaspGrid
     */
    template<int dim, class Coordinates, int codim>
    struct hasEntity< YaspGrid<dim, Coordinates>, codim>
    {
      static const bool v = true;
    };

    /**
     * \brief YaspGrid can iterate over all codimensions
     * \ingroup YaspGrid
     **/
    template<int dim, class Coordinates, int codim>
    struct hasEntityIterator<YaspGrid<dim, Coordinates>, codim>
    {
      static const bool v = true;
    };

    /** \brief YaspGrid can communicate on all codimensions
     *  \ingroup YaspGrid
     */
    template<int dim, int codim, class Coordinates>
    struct canCommunicate< YaspGrid< dim, Coordinates>, codim >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is levelwise conforming
       \ingroup YaspGrid
     */
    template<int dim, class Coordinates>
    struct isLevelwiseConforming< YaspGrid<dim, Coordinates> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is leafwise conforming
       \ingroup YaspGrid
     */
    template<int dim, class Coordinates>
    struct isLeafwiseConforming< YaspGrid<dim, Coordinates> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is viewThreadSafe
       \ingroup YaspGrid
     */
    template<int dim, class Coordinates>
    struct viewThreadSafe< YaspGrid<dim, Coordinates> >
    {
      static const bool v = true;
    };

  }

} // end namespace

// Include the specialization of the StructuredGridFactory class for YaspGrid
#include <dune/grid/yaspgrid/structuredyaspgridfactory.hh>
// Include the specialization of the BackupRestoreFacility class for YaspGrid
#include <dune/grid/yaspgrid/backuprestore.hh>

#endif
