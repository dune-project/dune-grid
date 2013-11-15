// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_HH
#define DUNE_GRID_YASPGRID_HH

#include <iostream>
#include <vector>
#include <algorithm>
#include <stack>

// either include stdint.h or provide fallback for uint8_t
#if HAVE_STDINT_H
#include <stdint.h>
#else
typedef unsigned char uint8_t;
#endif

#include <dune/grid/common/grid.hh>     // the grid base classes
#include <dune/grid/yaspgrid/grids.hh>  // the yaspgrid base classes
#include <dune/grid/common/capabilities.hh> // the capabilities
#include <dune/common/shared_ptr.hh>
#include <dune/common/bigunsignedint.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/grid/common/indexidset.hh>
#include <dune/grid/common/datahandleif.hh>


#if HAVE_MPI
#include <dune/common/parallel/mpicollectivecommunication.hh>
#endif

/*! \file yaspgrid.hh
   YaspGrid stands for yet another structured parallel grid.
   It will implement the dune grid interface for structured grids with codim 0
   and dim, with arbitrary overlap, parallel features with two overlap
   models, periodic boundaries and fast a implementation allowing on-the-fly computations.
 */

namespace Dune {

  //************************************************************************
  /*! define name for floating point type used for coordinates in yaspgrid.
     You can change the type for coordinates by changing this single typedef.
   */
  typedef double yaspgrid_ctype;

  /* some sizes for building global ids
   */
  const int yaspgrid_dim_bits = 24; // bits for encoding each dimension
  const int yaspgrid_level_bits = 6; // bits for encoding level number
  const int yaspgrid_codim_bits = 4; // bits for encoding codimension


  //************************************************************************
  // forward declaration of templates

  template<int dim>                             class YaspGrid;
  template<int mydim, int cdim, class GridImp>  class YaspGeometry;
  template<int codim, int dim, class GridImp>   class YaspEntity;
  template<int codim, class GridImp>            class YaspEntityPointer;
  template<int codim, class GridImp>            class YaspEntitySeed;
  template<int codim, PartitionIteratorType pitype, class GridImp> class YaspLevelIterator;
  template<class GridImp>            class YaspIntersectionIterator;
  template<class GridImp>            class YaspIntersection;
  template<class GridImp>            class YaspHierarchicIterator;
  template<class GridImp, bool isLeafIndexSet>                     class YaspIndexSet;
  template<class GridImp>            class YaspGlobalIdSet;

  namespace FacadeOptions
  {

    template<int dim, int mydim, int cdim>
    struct StoreGeometryReference<mydim, cdim, YaspGrid<dim>, YaspGeometry>
    {
      static const bool v = false;
    };

    template<int dim, int mydim, int cdim>
    struct StoreGeometryReference<mydim, cdim, const YaspGrid<dim>, YaspGeometry>
    {
      static const bool v = false;
    };

  }

} // namespace Dune

#include <dune/grid/yaspgrid/yaspgridgeometry.hh>
#include <dune/grid/yaspgrid/yaspgridentity.hh>
#include <dune/grid/yaspgrid/yaspgridintersection.hh>
#include <dune/grid/yaspgrid/yaspgridintersectioniterator.hh>
#include <dune/grid/yaspgrid/yaspgridhierarchiciterator.hh>
#include <dune/grid/yaspgrid/yaspgridentityseed.hh>
#include <dune/grid/yaspgrid/yaspgridentitypointer.hh>
#include <dune/grid/yaspgrid/yaspgridleveliterator.hh>
#include <dune/grid/yaspgrid/yaspgridindexsets.hh>
#include <dune/grid/yaspgrid/yaspgrididset.hh>

namespace Dune {

  template<int dim>
  struct YaspGridFamily
  {
#if HAVE_MPI
    typedef CollectiveCommunication<MPI_Comm> CCType;
#else
    typedef CollectiveCommunication<Dune::YaspGrid<dim> > CCType;
#endif

    typedef GridTraits<dim,                                     // dimension of the grid
        dim,                                                    // dimension of the world space
        Dune::YaspGrid<dim>,
        YaspGeometry,YaspEntity,
        YaspEntityPointer,
        YaspLevelIterator,                                      // type used for the level iterator
        YaspIntersection,              // leaf  intersection
        YaspIntersection,              // level intersection
        YaspIntersectionIterator,              // leaf  intersection iter
        YaspIntersectionIterator,              // level intersection iter
        YaspHierarchicIterator,
        YaspLevelIterator,                                      // type used for the leaf(!) iterator
        YaspIndexSet< const YaspGrid< dim >, false >,                  // level index set
        YaspIndexSet< const YaspGrid< dim >, true >,                  // leaf index set
        YaspGlobalIdSet<const YaspGrid<dim> >,
        bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+yaspgrid_codim_bits>,
        YaspGlobalIdSet<const YaspGrid<dim> >,
        bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+yaspgrid_codim_bits>,
        CCType,
        DefaultLevelGridViewTraits, DefaultLeafGridViewTraits,
        YaspEntitySeed>
    Traits;
  };

  template<int dim, int codim>
  struct YaspCommunicateMeta {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,codim))
      {
        DUNE_THROW(GridError, "interface communication not implemented");
      }
      YaspCommunicateMeta<dim,codim-1>::comm(g,data,iftype,dir,level);
    }
  };

  template<int dim>
  struct YaspCommunicateMeta<dim,dim> {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,dim))
        g.template communicateCodim<DataHandle,dim>(data,iftype,dir,level);
      YaspCommunicateMeta<dim,dim-1>::comm(g,data,iftype,dir,level);
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

  //************************************************************************
  /*!
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief Provides a distributed structured cube mesh.
     \ingroup GridImplementations
     \ingroup YaspGrid

     YaspGrid stands for yet another structured parallel grid.
     It implements the dune grid interface for structured grids with codim 0
     and dim, with arbitrary overlap (including zero),
     periodic boundaries and fast implementation allowing on-the-fly computations.

     \tparam dim The dimension of the grid and its surrounding world

     \par History:
     \li started on July 31, 2004 by PB based on abstractions developed in summer 2003
   */
  template<int dim>
  class YaspGrid
    : public GridDefaultImplementation<dim,dim,yaspgrid_ctype,YaspGridFamily<dim> >
  {
  public:
    //! Type used for coordinates
    typedef yaspgrid_ctype ctype;

    struct Intersection {
      /** \brief The intersection as a subgrid of the local grid */
      SubYGrid<dim,ctype> grid;

      /** \brief Rank of the process where the other grid is stored */
      int rank;

      /** \brief Manhattan distance to the other grid */
      int distance;
    };

    /** \brief A single grid level within a YaspGrid
     */
    struct YGridLevel {

      /** \brief Level number of this level grid */
      int level() const
      {
        return level_;
      }

      // cell (codim 0) data
      YGrid<dim,ctype> cell_global;         // the whole cell grid on that level
      SubYGrid<dim,ctype> cell_overlap;     // we have no ghost cells, so our part is overlap completely
      SubYGrid<dim,ctype> cell_interior;    // interior cells are a subgrid of all cells

      std::deque<Intersection> send_cell_overlap_overlap;  // each intersection is a subgrid of overlap
      std::deque<Intersection> recv_cell_overlap_overlap;  // each intersection is a subgrid of overlap

      std::deque<Intersection> send_cell_interior_overlap; // each intersection is a subgrid of overlap
      std::deque<Intersection> recv_cell_overlap_interior; // each intersection is a subgrid of overlap

      // vertex (codim dim) data
      YGrid<dim,ctype> vertex_global;           // the whole vertex grid on that level
      SubYGrid<dim,ctype> vertex_overlapfront;  // all our vertices are overlap and front
      SubYGrid<dim,ctype> vertex_overlap;       // subgrid containing only overlap
      SubYGrid<dim,ctype> vertex_interiorborder; // subgrid containing only interior and border
      SubYGrid<dim,ctype> vertex_interior;      // subgrid containing only interior

      std::deque<Intersection> send_vertex_overlapfront_overlapfront; // each intersection is a subgrid of overlapfront
      std::deque<Intersection> recv_vertex_overlapfront_overlapfront; // each intersection is a subgrid of overlapfront

      std::deque<Intersection> send_vertex_overlap_overlapfront; // each intersection is a subgrid of overlapfront
      std::deque<Intersection> recv_vertex_overlapfront_overlap; // each intersection is a subgrid of overlapfront

      std::deque<Intersection> send_vertex_interiorborder_interiorborder; // each intersection is a subgrid of overlapfront
      std::deque<Intersection> recv_vertex_interiorborder_interiorborder; // each intersection is a subgrid of overlapfront

      std::deque<Intersection> send_vertex_interiorborder_overlapfront; // each intersection is a subgrid of overlapfront
      std::deque<Intersection> recv_vertex_overlapfront_interiorborder; // each intersection is a subgrid of overlapfront

      // general
      YaspGrid<dim>* mg;  // each grid level knows its multigrid
      int overlap;           // in mesh cells on this level
      /** \brief The level number within the YaspGrid level hierarchy */
      int level_;
    };

    //! define types used for arguments
    typedef FieldVector<int, dim> iTupel;
    typedef FieldVector<ctype, dim> fTupel;

    // communication tag used by multigrid
    enum { tag = 17 };

    //! return reference to torus
    const Torus<dim>& torus () const
    {
      return _torus;
    }

    //! Iterator over the grid levels
    typedef typename ReservedVector<YGridLevel,32>::const_iterator YGridLevelIterator;

    //! return iterator pointing to coarsest level
    YGridLevelIterator begin () const
    {
      return YGridLevelIterator(_levels,0);
    }

    //! return iterator pointing to given level
    YGridLevelIterator begin (int i) const
    {
      if (i<0 || i>maxLevel())
        DUNE_THROW(GridError, "level not existing");
      return YGridLevelIterator(_levels,i);
    }

    //! return iterator pointing to one past the finest level
    YGridLevelIterator end () const
    {
      return YGridLevelIterator(_levels,maxLevel()+1);
    }

    // static method to create the default load balance strategy
    static const YLoadBalance<dim>* defaultLoadbalancer()
    {
      static YLoadBalance<dim> lb;
      return & lb;
    }

  protected:
    /** \brief Make a new YGridLevel structure
     *
     * \param L           size of the whole domain in each direction
     * \param s           number of cells in each direction
     * \param periodic    indicate periodicity for each direction
     * \param o_interior  origin of interior (non-overlapping) cell decomposition
     * \param s_interior  size of interior cell decomposition
     * \param overlap     to be used on this grid level
     */
    YGridLevel makelevel (int level, fTupel L, iTupel s, std::bitset<dim> periodic, iTupel o_interior, iTupel s_interior, int overlap)
    {
      // first, lets allocate a new structure
      YGridLevel g;
      g.overlap = overlap;
      g.mg = this;
      g.level_ = level;

      // the global cell grid
      iTupel o = iTupel(0); // logical origin is always 0, that is not a restriction
      fTupel h;
      fTupel r;
      for (int i=0; i<dim; i++) h[i] = L[i]/s[i]; // the mesh size in each direction
      for (int i=0; i<dim; i++) r[i] = 0.5*h[i];  // the shift for cell centers
      g.cell_global = YGrid<dim,ctype>(o,s,h,r);     // this is the global cell grid

      // extend the cell interior grid by overlap considering periodicity
      iTupel o_overlap;
      iTupel s_overlap;
      for (int i=0; i<dim; i++)
      {
        if (periodic[i])
        {
          // easy case, extend by 2 overlaps in total
          o_overlap[i] = o_interior[i]-overlap;      // Note: origin might be negative now
          s_overlap[i] = s_interior[i]+2*overlap;    // Note: might be larger than global size
        }
        else
        {
          // nonperiodic case, intersect with global size
          int min = std::max(0,o_interior[i]-overlap);
          int max = std::min(s[i]-1,o_interior[i]+s_interior[i]-1+overlap);
          o_overlap[i] = min;
          s_overlap[i] = max-min+1;
        }
      }
      g.cell_overlap = SubYGrid<dim,ctype>(YGrid<dim,ctype>(o_overlap,s_overlap,h,r));

      // now make the interior grid a subgrid of the overlapping grid
      iTupel offset;
      for (int i=0; i<dim; i++) offset[i] = o_interior[i]-o_overlap[i];
      g.cell_interior = SubYGrid<dim,ctype>(o_interior,s_interior,offset,s_overlap,h,r);

      // compute cell intersections
      intersections(g.cell_overlap,g.cell_overlap,g.cell_global.size(),g.send_cell_overlap_overlap,g.recv_cell_overlap_overlap);
      intersections(g.cell_interior,g.cell_overlap,g.cell_global.size(),g.send_cell_interior_overlap,g.recv_cell_overlap_interior);

      // now we can do the vertex grids. They are derived completely from the cell grids
      iTupel o_vertex_global, s_vertex_global;
      for (int i=0; i<dim; i++) r[i] = 0.0;  // the shift for vertices is zero, and the mesh size is same as for cells

      // first let's make the global grid
      for (int i=0; i<dim; i++) o_vertex_global[i] = g.cell_global.origin(i);
      for (int i=0; i<dim; i++) s_vertex_global[i] = g.cell_global.size(i)+1; // one more vertices than cells ...
      g.vertex_global = YGrid<dim,ctype>(o_vertex_global,s_vertex_global,h,r);

      // now the local grid stored in this processor. All other grids are subgrids of this
      iTupel o_vertex_overlapfront;
      iTupel s_vertex_overlapfront;
      for (int i=0; i<dim; i++) o_vertex_overlapfront[i] = g.cell_overlap.origin(i);
      for (int i=0; i<dim; i++) s_vertex_overlapfront[i] = g.cell_overlap.size(i)+1; // one more vertices than cells ...
      g.vertex_overlapfront = SubYGrid<dim,ctype>(YGrid<dim,ctype>(o_vertex_overlapfront,s_vertex_overlapfront,h,r));

      // now overlap only (i.e. without front), is subgrid of overlapfront
      iTupel o_vertex_overlap;
      iTupel s_vertex_overlap;
      for (int i=0; i<dim; i++)
      {
        o_vertex_overlap[i] = g.cell_overlap.origin(i);
        s_vertex_overlap[i] = g.cell_overlap.size(i)+1;

        if (!periodic[i] && g.cell_overlap.origin(i)>g.cell_global.origin(i))
        {
          // not at the lower boundary
          o_vertex_overlap[i] += 1;
          s_vertex_overlap[i] -= 1;
        }

        if (!periodic[i] && g.cell_overlap.origin(i)+g.cell_overlap.size(i)<g.cell_global.origin(i)+g.cell_global.size(i))
        {
          // not at the upper boundary
          s_vertex_overlap[i] -= 1;
        }


        offset[i] = o_vertex_overlap[i]-o_vertex_overlapfront[i];
      }
      g.vertex_overlap = SubYGrid<dim,ctype>(o_vertex_overlap,s_vertex_overlap,offset,s_vertex_overlapfront,h,r);

      // now interior with border
      iTupel o_vertex_interiorborder;
      iTupel s_vertex_interiorborder;
      for (int i=0; i<dim; i++) o_vertex_interiorborder[i] = g.cell_interior.origin(i);
      for (int i=0; i<dim; i++) s_vertex_interiorborder[i] = g.cell_interior.size(i)+1;
      for (int i=0; i<dim; i++) offset[i] = o_vertex_interiorborder[i]-o_vertex_overlapfront[i];
      g.vertex_interiorborder = SubYGrid<dim,ctype>(o_vertex_interiorborder,s_vertex_interiorborder,offset,s_vertex_overlapfront,h,r);

      // now only interior
      iTupel o_vertex_interior;
      iTupel s_vertex_interior;
      for (int i=0; i<dim; i++)
      {
        o_vertex_interior[i] = g.cell_interior.origin(i);
        s_vertex_interior[i] = g.cell_interior.size(i)+1;

        if (!periodic[i] && g.cell_interior.origin(i)>g.cell_global.origin(i))
        {
          // not at the lower boundary
          o_vertex_interior[i] += 1;
          s_vertex_interior[i] -= 1;
        }

        if (!periodic[i] && g.cell_interior.origin(i)+g.cell_interior.size(i)<g.cell_global.origin(i)+g.cell_global.size(i))
        {
          // not at the upper boundary
          s_vertex_interior[i] -= 1;
        }

        offset[i] = o_vertex_interior[i]-o_vertex_overlapfront[i];
      }
      g.vertex_interior = SubYGrid<dim,ctype>(o_vertex_interior,s_vertex_interior,offset,s_vertex_overlapfront,h,r);

      // compute vertex intersections
      intersections(g.vertex_overlapfront,g.vertex_overlapfront,g.cell_global.size(),
                    g.send_vertex_overlapfront_overlapfront,g.recv_vertex_overlapfront_overlapfront);
      intersections(g.vertex_overlap,g.vertex_overlapfront,g.cell_global.size(),
                    g.send_vertex_overlap_overlapfront,g.recv_vertex_overlapfront_overlap);
      intersections(g.vertex_interiorborder,g.vertex_interiorborder,g.cell_global.size(),
                    g.send_vertex_interiorborder_interiorborder,g.recv_vertex_interiorborder_interiorborder);
      intersections(g.vertex_interiorborder,g.vertex_overlapfront,g.cell_global.size(),
                    g.send_vertex_interiorborder_overlapfront,g.recv_vertex_overlapfront_interiorborder);

      // return the whole thing
      return g;
    }


    struct mpifriendly_ygrid {
      mpifriendly_ygrid ()
        : origin(0), size(0), h(0.0), r(0.0)
      {}
      mpifriendly_ygrid (const YGrid<dim,ctype>& grid)
        : origin(grid.origin()), size(grid.size()), h(grid.meshsize()), r(grid.shift())
      {}
      iTupel origin;
      iTupel size;
      fTupel h;
      fTupel r;
    };

    /** \brief Construct list of intersections with neighboring processors
     *
     * \param recvgrid the grid stored in this processor
     * \param sendgrid  the subgrid to be sent to neighboring processors
     * \param size needed to shift local grid in periodic case
     * \returns two lists: Intersections to be sent and Intersections to be received
     * \note sendgrid/recvgrid may be SubYGrids. Since intersection method is virtual it should work properly
     */
    void intersections (const SubYGrid<dim,ctype>& sendgrid, const SubYGrid<dim,ctype>& recvgrid, const iTupel& size,
                        std::deque<Intersection>& sendlist, std::deque<Intersection>& recvlist)
    {
      // the exchange buffers
      std::vector<YGrid<dim,ctype> > send_recvgrid(_torus.neighbors());
      std::vector<YGrid<dim,ctype> > recv_recvgrid(_torus.neighbors());
      std::vector<YGrid<dim,ctype> > send_sendgrid(_torus.neighbors());
      std::vector<YGrid<dim,ctype> > recv_sendgrid(_torus.neighbors());

      // new exchange buffers to send simple struct without virtual functions
      std::vector<mpifriendly_ygrid> mpifriendly_send_recvgrid(_torus.neighbors());
      std::vector<mpifriendly_ygrid> mpifriendly_recv_recvgrid(_torus.neighbors());
      std::vector<mpifriendly_ygrid> mpifriendly_send_sendgrid(_torus.neighbors());
      std::vector<mpifriendly_ygrid> mpifriendly_recv_sendgrid(_torus.neighbors());

      // fill send buffers; iterate over all neighboring processes
      // non-periodic case is handled automatically because intersection will be zero
      for (typename Torus<dim>::ProcListIterator i=_torus.sendbegin(); i!=_torus.sendend(); ++i)
      {
        // determine if we communicate with this neighbor (and what)
        bool skip = false;
        iTupel coord = _torus.coord();   // my coordinates
        iTupel delta = i.delta();        // delta to neighbor
        iTupel nb = coord;               // the neighbor
        for (int k=0; k<dim; k++) nb[k] += delta[k];
        iTupel v = iTupel(0);                    // grid movement

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
          send_sendgrid[i.index()] = YGrid<dim,ctype>(iTupel(0),iTupel(0),fTupel(0.0),fTupel(0.0));
          send_recvgrid[i.index()] = YGrid<dim,ctype>(iTupel(0),iTupel(0),fTupel(0.0),fTupel(0.0));
        }
      }

      // issue send requests for sendgrid being sent to all neighbors
      for (typename Torus<dim>::ProcListIterator i=_torus.sendbegin(); i!=_torus.sendend(); ++i)
      {
        mpifriendly_send_sendgrid[i.index()] = mpifriendly_ygrid(send_sendgrid[i.index()]);
        _torus.send(i.rank(), &mpifriendly_send_sendgrid[i.index()], sizeof(mpifriendly_ygrid));
      }

      // issue recv requests for sendgrids of neighbors
      for (typename Torus<dim>::ProcListIterator i=_torus.recvbegin(); i!=_torus.recvend(); ++i)
        _torus.recv(i.rank(), &mpifriendly_recv_sendgrid[i.index()], sizeof(mpifriendly_ygrid));

      // exchange the sendgrids
      _torus.exchange();

      // issue send requests for recvgrid being sent to all neighbors
      for (typename Torus<dim>::ProcListIterator i=_torus.sendbegin(); i!=_torus.sendend(); ++i)
      {
        mpifriendly_send_recvgrid[i.index()] = mpifriendly_ygrid(send_recvgrid[i.index()]);
        _torus.send(i.rank(), &mpifriendly_send_recvgrid[i.index()], sizeof(mpifriendly_ygrid));
      }

      // issue recv requests for recvgrid of neighbors
      for (typename Torus<dim>::ProcListIterator i=_torus.recvbegin(); i!=_torus.recvend(); ++i)
        _torus.recv(i.rank(), &mpifriendly_recv_recvgrid[i.index()], sizeof(mpifriendly_ygrid));

      // exchange the recvgrid
      _torus.exchange();

      // process receive buffers and compute intersections
      for (typename Torus<dim>::ProcListIterator i=_torus.recvbegin(); i!=_torus.recvend(); ++i)
      {
        // what must be sent to this neighbor
        Intersection send_intersection;
        mpifriendly_ygrid yg = mpifriendly_recv_recvgrid[i.index()];
        recv_recvgrid[i.index()] = YGrid<dim,ctype>(yg.origin,yg.size,yg.h,yg.r);
        send_intersection.grid = sendgrid.intersection(recv_recvgrid[i.index()]);
        send_intersection.rank = i.rank();
        send_intersection.distance = i.distance();
        if (!send_intersection.grid.empty()) sendlist.push_front(send_intersection);

        Intersection recv_intersection;
        yg = mpifriendly_recv_sendgrid[i.index()];
        recv_sendgrid[i.index()] = YGrid<dim,ctype>(yg.origin,yg.size,yg.h,yg.r);
        recv_intersection.grid = recvgrid.intersection(recv_sendgrid[i.index()]);
        recv_intersection.rank = i.rank();
        recv_intersection.distance = i.distance();
        if(!recv_intersection.grid.empty()) recvlist.push_back(recv_intersection);
      }
    }

  protected:

    typedef const YaspGrid<dim> GridImp;

    void init()
    {
      setsizes();
      indexsets.push_back( make_shared< YaspIndexSet<const YaspGrid<dim>, false > >(*this,0) );
      boundarysegmentssize();
    }

    void boundarysegmentssize()
    {
      // sizes of local macro grid
      const FieldVector<int, dim> & size = begin()->cell_overlap.size();
      Dune::array<int, dim> sides;
      {
        for (int i=0; i<dim; i++)
        {
          sides[i] =
            ((begin()->cell_overlap.origin(i) == begin()->cell_global.origin(i))+
             (begin()->cell_overlap.origin(i) + begin()->cell_overlap.size(i)
                    == begin()->cell_global.origin(i) + begin()->cell_global.size(i)));
        }
      }
      nBSegments = 0;
      for (int k=0; k<dim; k++)
      {
        int offset = 1;
        for (int l=0; l<dim; l++)
        {
          if (l==k) continue;
          offset *= size[l];
        }
        nBSegments += sides[k]*offset;
      }
    }

  public:

    // define the persistent index type
    typedef bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+yaspgrid_codim_bits> PersistentIndexType;

    //! the GridFamily of this grid
    typedef YaspGridFamily<dim> GridFamily;
    // the Traits
    typedef typename YaspGridFamily<dim>::Traits Traits;

    // need for friend declarations in entity
    typedef YaspIndexSet<YaspGrid<dim>, false > LevelIndexSetType;
    typedef YaspIndexSet<YaspGrid<dim>, true > LeafIndexSetType;
    typedef YaspGlobalIdSet<YaspGrid<dim> > GlobalIdSetType;

    //! shorthand for some data types
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;
    typedef typename std::deque<Intersection>::const_iterator ISIT;

    //! The constructor of the old MultiYGrid class
    void MultiYGridSetup (
                fTupel L, iTupel s, std::bitset<dim> periodic, int overlap, const YLoadBalance<dim>* lb = defaultLoadbalancer())
    {
      _LL = L;
      _s = s;
      _periodic = periodic;
      _levels.resize(1);
      _overlap = overlap;

      // coarse cell interior  grid obtained through partitioning of global grid
#if HAVE_MPI
      iTupel o_interior;
      iTupel s_interior;
      iTupel o = iTupel(0);
      array<int,dim> sArray;
      std::copy(s.begin(), s.end(), sArray.begin());
      double imbal = _torus.partition(_torus.rank(),o,sArray,o_interior,s_interior);
      imbal = _torus.global_max(imbal);
#else
      iTupel o = iTupel(0);
      iTupel o_interior(o);
      iTupel s_interior(s);
#endif
      // add level
      _levels[0] = makelevel(0,L,s,periodic,o_interior,s_interior,overlap);
    }

    //! The constructor of the old MultiYGrid class
    void MultiYGridSetup (
      fTupel L,
      Dune::array<int,dim> s,
      std::bitset<dim> periodic,
      int overlap,
      const YLoadBalance<dim>* lb = defaultLoadbalancer())
    {
      _LL = L;
      _periodic = periodic;
      _levels.resize(1);
      _overlap = overlap;

      std::copy(s.begin(), s.end(), this->_s.begin());

      // coarse cell interior grid obtained through partitioning of global grid
      iTupel o = iTupel(0);
      iTupel o_interior(o);
      iTupel s_interior;
      std::copy(s.begin(), s.end(), s_interior.begin());
#if HAVE_MPI
      double imbal = _torus.partition(_torus.rank(),o,s,o_interior,s_interior);
      imbal = _torus.global_max(imbal);
#endif

      // add level
      _levels[0] = makelevel(0,L,_s,periodic,o_interior,s_interior,overlap);
    }

    /*! Constructor
       @param comm MPI communicator where this mesh is distributed to
       @param L extension of the domain
       @param s number of cells on coarse mesh in each direction
       @param periodic tells if direction is periodic or not
       @param overlap size of overlap on coarsest grid (same in all directions)
       @param lb pointer to an overloaded YLoadBalance instance

       \deprecated Will be removed after dune-grid 2.3.
         Use the corresponding constructor taking array<int> and std::bitset instead.
     */
    YaspGrid (Dune::MPIHelper::MPICommunicator comm,
              Dune::FieldVector<ctype, dim> L,
              Dune::FieldVector<int, dim> s,
              Dune::FieldVector<bool, dim> periodic, int overlap,
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
    DUNE_DEPRECATED_MSG("Use the corresponding constructor taking array<int> and std::bitset")
#if HAVE_MPI
      : ccobj(comm),
        _torus(comm,tag,s,lb),
#else
      :  _torus(tag,s,lb),
#endif
        leafIndexSet_(*this),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
    {
      MultiYGridSetup(L,s,std::bitset<dim>(),overlap,lb);

      // hack: copy input bitfield (in FieldVector<bool>) into std::bitset
      for (size_t i=0; i<dim; i++)
        this->_periodic[i] = periodic[i];
      init();
    }


    /*! Constructor for a sequential YaspGrid

       Sequential here means that the whole grid is living on one process even if your program is running
       in parallel.
       @see YaspGrid(Dune::MPIHelper::MPICommunicator, Dune::FieldVector<ctype, dim>, Dune::FieldVector<int, dim>,  Dune::FieldVector<bool, dim>, int)
       for constructing one parallel grid decomposed between the processors.
       @param L extension of the domain
       @param s number of cells on coarse mesh in each direction
       @param periodic tells if direction is periodic or not
       @param overlap size of overlap on coarsest grid (same in all directions)
       @param lb pointer to an overloaded YLoadBalance instance

       \deprecated Will be removed after dune-grid 2.3.
         Use the corresponding constructor taking array<int> and std::bitset instead.
     */
    YaspGrid (Dune::FieldVector<ctype, dim> L,
              Dune::FieldVector<int, dim> s,
              Dune::FieldVector<bool, dim> periodic, int overlap,
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
    DUNE_DEPRECATED_MSG("Use the corresponding constructor taking array<int> and std::bitset")
#if HAVE_MPI
      : ccobj(MPI_COMM_SELF),
        _torus(MPI_COMM_SELF,tag,s,lb),
#else
      : _torus(tag,s,lb),
#endif
        leafIndexSet_(*this),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
    {
      MultiYGridSetup(L,s,std::bitset<dim>(),overlap,lb);

      // hack: copy input bitfield (in FieldVector<bool>) into std::bitset
      for (size_t i=0; i<dim; i++)
        this->_periodic[i] = periodic[i];
      init();
    }

    /*! Constructor
       @param comm MPI communicator where this mesh is distributed to
       @param L extension of the domain
       @param s number of cells on coarse mesh in each direction
       @param periodic tells if direction is periodic or not
       @param overlap size of overlap on coarsest grid (same in all directions)
       @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (Dune::MPIHelper::MPICommunicator comm,
              Dune::FieldVector<ctype, dim> L,
              Dune::array<int, dim> s,
              std::bitset<dim> periodic,
              int overlap,
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
#if HAVE_MPI
      : ccobj(comm),
        _torus(comm,tag,s,lb),
#else
      : _torus(tag,s,lb),
#endif
        leafIndexSet_(*this),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
    {
      MultiYGridSetup(L,s,periodic,overlap,lb);

      init();
    }


    /*! Constructor for a sequential YaspGrid

       Sequential here means that the whole grid is living on one process even if your program is running
       in parallel.
       @see YaspGrid(Dune::MPIHelper::MPICommunicator, Dune::FieldVector<ctype, dim>, Dune::FieldVector<int, dim>,  Dune::FieldVector<bool, dim>, int)
       for constructing one parallel grid decomposed between the processors.
       @param L extension of the domain
       @param s number of cells on coarse mesh in each direction
       @param periodic tells if direction is periodic or not
       @param overlap size of overlap on coarsest grid (same in all directions)
       @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (Dune::FieldVector<ctype, dim> L,
              Dune::array<int, dim> s,
              std::bitset<dim> periodic,
              int overlap,
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
#if HAVE_MPI
      : ccobj(MPI_COMM_SELF),
        _torus(MPI_COMM_SELF,tag,s,lb),
#else
      : _torus(tag,s,lb),
#endif
        leafIndexSet_(*this),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
    {
      MultiYGridSetup(L,s,periodic,overlap,lb);

      init();
    }

    /*! Constructor for a sequential YaspGrid without periodicity

       Sequential here means that the whole grid is living on one process even if your program is running
       in parallel.
       @see YaspGrid(Dune::MPIHelper::MPICommunicator, Dune::FieldVector<ctype, dim>, Dune::FieldVector<int, dim>,  Dune::FieldVector<bool, dim>, int)
       for constructing one parallel grid decomposed between the processors.
       @param L extension of the domain (lower left is always (0,...,0)
       @param elements number of cells on coarse mesh in each direction
     */
    YaspGrid (Dune::FieldVector<ctype, dim> L,
              Dune::array<int, dim> elements)
#if HAVE_MPI
      : ccobj(MPI_COMM_SELF),
        _torus(MPI_COMM_SELF,tag,elements,defaultLoadbalancer()),
#else
      : _torus(tag,elements,defaultLoadbalancer()),
#endif
        leafIndexSet_(*this),
        _LL(L),
        _overlap(0),
        keep_ovlp(true),
        adaptRefCount(0), adaptActive(false)
    {
      _levels.resize(1);

      std::copy(elements.begin(), elements.end(), _s.begin());

      // coarse cell interior grid obtained through partitioning of global grid
      iTupel o = iTupel(0);
      iTupel o_interior(o);
      iTupel s_interior;
      std::copy(elements.begin(), elements.end(), s_interior.begin());
#if HAVE_MPI
      double imbal = _torus.partition(_torus.rank(),o,elements,o_interior,s_interior);
      imbal = _torus.global_max(imbal);
#endif

      // add level
      _levels[0] = makelevel(0,L,_s,_periodic,o_interior,s_interior,0);

      init();
    }

  private:
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

    //! refine the grid refCount times. What about overlap?
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

        setsizes();
        indexsets.pop_back();
      }

      // If refCount is positive refine the grid
      for (int k=0; k<refCount; k++)
      {
        // access to coarser grid level
        YGridLevel& cg = _levels[maxLevel()];

        // compute size of new global grid
        iTupel s;
        for (int i=0; i<dim; i++)
          s[i] = 2*cg.cell_global.size(i);

        // compute overlap
        int overlap = (keep_ovlp) ? 2*cg.overlap : cg.overlap;

        // the cell interior grid obtained from coarse cell interior grid
        iTupel o_interior;
        iTupel s_interior;
        for (int i=0; i<dim; i++)
          o_interior[i] = 2*cg.cell_interior.origin(i);
        for (int i=0; i<dim; i++)
          s_interior[i] = 2*cg.cell_interior.size(i);

        // add level
        _levels.push_back( makelevel(_levels.size(),_LL,s,_periodic,o_interior,s_interior,overlap) );

        setsizes();
        indexsets.push_back( make_shared<YaspIndexSet<const YaspGrid<dim>, false > >(*this,maxLevel()) );
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

    // \brief obtain EntityPointer from EntitySeed. */
    template <typename Seed>
    typename Traits::template Codim<Seed::codimension>::EntityPointer
    entityPointer(const Seed& seed) const
    {
      const int codim = Seed::codimension;
      YGridLevelIterator g = begin(this->getRealImplementation(seed).level());
      switch (codim)
      {
      case 0 :
        return YaspEntityPointer<codim,GridImp>(this,g,
                                                TSI(g->cell_overlap, this->getRealImplementation(seed).coord()));
      case dim :
        return YaspEntityPointer<codim,GridImp>(this,g,
                                                TSI(g->vertex_overlap, this->getRealImplementation(seed).coord()));
      default :
        DUNE_THROW(GridError, "YaspEntityPointer: codim not implemented");
      }
    }

    //! return size (= distance in graph) of overlap region
    int overlapSize (int level, int codim) const
    {
      YGridLevelIterator g = begin(level);
      return g->overlap;
    }

    //! return size (= distance in graph) of overlap region
    int overlapSize (int codim) const
    {
      YGridLevelIterator g = begin(maxLevel());
      return g->overlap;
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize (int level, int codim) const
    {
      return 0;
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize (int codim) const
    {
      return 0;
    }

    //! number of entities per level and codim in this process
    int size (int level, int codim) const
    {
      return sizes[level][codim];
    }

    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return sizes[maxLevel()][codim];
    }

    //! number of entities per level and geometry type in this process
    int size (int level, GeometryType type) const
    {
      return (type.isCube()) ? sizes[level][dim-type.dim()] : 0;
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
      const std::deque<Intersection>* sendlist=0;
      const std::deque<Intersection>* recvlist=0;
      if (codim==0) // the elements
      {
        if (iftype==InteriorBorder_InteriorBorder_Interface)
          return; // there is nothing to do in this case
        if (iftype==InteriorBorder_All_Interface)
        {
          sendlist = &g->send_cell_interior_overlap;
          recvlist = &g->recv_cell_overlap_interior;
        }
        if (iftype==Overlap_OverlapFront_Interface || iftype==Overlap_All_Interface || iftype==All_All_Interface)
        {
          sendlist = &g->send_cell_overlap_overlap;
          recvlist = &g->recv_cell_overlap_overlap;
        }
      }
      if (codim==dim) // the vertices
      {
        if (iftype==InteriorBorder_InteriorBorder_Interface)
        {
          sendlist = &g->send_vertex_interiorborder_interiorborder;
          recvlist = &g->recv_vertex_interiorborder_interiorborder;
        }

        if (iftype==InteriorBorder_All_Interface)
        {
          sendlist = &g->send_vertex_interiorborder_overlapfront;
          recvlist = &g->recv_vertex_overlapfront_interiorborder;
        }
        if (iftype==Overlap_OverlapFront_Interface || iftype==Overlap_All_Interface)
        {
          sendlist = &g->send_vertex_overlap_overlapfront;
          recvlist = &g->recv_vertex_overlapfront_overlap;
        }
        if (iftype==All_All_Interface)
        {
          sendlist = &g->send_vertex_overlapfront_overlapfront;
          recvlist = &g->recv_vertex_overlapfront_overlapfront;
        }
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
      if (data.fixedsize(dim,codim))
      {
        // fixed size: just take a dummy entity, size can be computed without communication
        cnt=0;
        for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          send_size[cnt] = is->grid.totalsize() * data.size(*it);
          cnt++;
        }
        cnt=0;
        for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          recv_size[cnt] = is->grid.totalsize() * data.size(*it);
          cnt++;
        }
      }
      else
      {
        // variable size case: sender side determines the size
        cnt=0;
        for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          // allocate send buffer for sizes per entitiy
          size_t *buf = new size_t[is->grid.totalsize()];
          send_sizes[cnt] = buf;

          // loop over entities and ask for size
          int i=0; size_t n=0;
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
          for ( ; it!=tsubend; ++it)
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
        for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
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
        for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          delete[] send_sizes[cnt];
          send_sizes[cnt] = 0;
          cnt++;
        }

        // process receive size buffers
        cnt=0;
        for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
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
      for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
      {
        // allocate send buffer
        DataType *buf = new DataType[send_size[cnt]];

        // remember send buffer
        sends[cnt] = buf;

        // make a message buffer
        MessageBuffer<DataType> mb(buf);

        // fill send buffer; iterate over cells in intersection
        typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
        it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
        typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
        tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
        for ( ; it!=tsubend; ++it)
          data.gather(mb,*it);

        // hand over send request to torus class
        torus().send(is->rank,buf,send_size[cnt]*sizeof(DataType));
        cnt++;
      }

      // allocate recv buffers and store receive request
      std::vector<DataType*> recvs(recvlist->size(),static_cast<DataType*>(0)); // store pointers to send buffers
      cnt=0;
      for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
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
      for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
      {
        delete[] sends[cnt];
        sends[cnt] = 0;
        cnt++;
      }

      // process receive buffers and delete them
      cnt=0;
      for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
      {
        // get recv buffer
        DataType *buf = recvs[cnt];

        // make a message buffer
        MessageBuffer<DataType> mb(buf);

        // copy data from receive buffer; iterate over cells in intersection
        if (data.fixedsize(dim,codim))
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          size_t n=data.size(*it);
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
          for ( ; it!=tsubend; ++it)
            data.scatter(mb,*it,n);
        }
        else
        {
          int i=0;
          size_t *sbuf = recv_sizes[cnt];
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
          for ( ; it!=tsubend; ++it)
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

#if HAVE_MPI
    /*! @brief return a collective communication object
     */
    const CollectiveCommunication<MPI_Comm>& comm () const
    {
      return ccobj;
    }
#else
    /*! @brief return a collective communication object
     */
    const CollectiveCommunication<YaspGrid>& comm () const
    {
      return ccobj;
    }
#endif

  private:

    // number of boundary segments of the level 0 grid
    int nBSegments;

    // Index classes need access to the real entity
    friend class Dune::YaspIndexSet<const Dune::YaspGrid<dim>, true >;
    friend class Dune::YaspIndexSet<const Dune::YaspGrid<dim>, false >;
    friend class Dune::YaspGlobalIdSet<const Dune::YaspGrid<dim> >;

    friend class Dune::YaspIntersectionIterator<const Dune::YaspGrid<dim> >;
    friend class Dune::YaspIntersection<const Dune::YaspGrid<dim> >;
    friend class Dune::YaspEntity<0, dim, const Dune::YaspGrid<dim> >;

    template <int codim_, class GridImp_>
    friend class Dune::YaspEntityPointer;

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
        dune_static_assert(( is_same<DT,Y>::value ), "DataType mismatch");
        a[i++] = data;
      }

      // read data from message buffer, acts like a stream !
      template<class Y>
      void read (Y& data) const
      {
        dune_static_assert(( is_same<DT,Y>::value ), "DataType mismatch");
        data = a[j++];
      }

    private:
      DT *a;
      int i;
      mutable int j;
    };

    void setsizes ()
    {
      for (YGridLevelIterator g=begin(); g!=end(); ++g)
      {
        // codim 0 (elements)
        sizes[g->level()][0] = 1;
        for (int i=0; i<dim; ++i)
          sizes[g->level()][0] *= g->cell_overlap.size(i);

        // codim 1 (faces)
        if (dim>1)
        {
          sizes[g->level()][1] = 0;
          for (int i=0; i<dim; ++i)
          {
            int s=g->cell_overlap.size(i)+1;
            for (int j=0; j<dim; ++j)
              if (j!=i)
                s *= g->cell_overlap.size(j);
            sizes[g->level()][1] += s;
          }
        }

        // codim dim-1 (edges)
        if (dim>2)
        {
          sizes[g->level()][dim-1] = 0;
          for (int i=0; i<dim; ++i)
          {
            int s=g->cell_overlap.size(i);
            for (int j=0; j<dim; ++j)
              if (j!=i)
                s *= g->cell_overlap.size(j)+1;
            sizes[g->level()][dim-1] += s;
          }
        }

        // codim dim (vertices)
        sizes[g->level()][dim] = 1;
        for (int i=0; i<dim; ++i)
          sizes[g->level()][dim] *= g->vertex_overlapfront.size(i);
      }
    }

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    YaspLevelIterator<cd,pitype,GridImp> levelbegin (int level) const
    {
      dune_static_assert( cd == dim || cd == 0 ,
                          "YaspGrid only supports Entities with codim=dim and codim=0");
      YGridLevelIterator g = begin(level);
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");
      if (pitype==Ghost_Partition)
        return levelend <cd, pitype> (level);
      if (cd==0)   // the elements
      {
        if (pitype<=InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->cell_interior.tsubbegin());
        if (pitype<=All_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->cell_overlap.tsubbegin());
      }
      if (cd==dim)   // the vertices
      {
        if (pitype==Interior_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->vertex_interior.tsubbegin());
        if (pitype==InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->vertex_interiorborder.tsubbegin());
        if (pitype==Overlap_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->vertex_overlap.tsubbegin());
        if (pitype<=All_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->vertex_overlapfront.tsubbegin());
      }
      DUNE_THROW(GridError, "YaspLevelIterator with this codim or partition type not implemented");
    }

    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    YaspLevelIterator<cd,pitype,GridImp> levelend (int level) const
    {
      dune_static_assert( cd == dim || cd == 0 ,
                          "YaspGrid only supports Entities with codim=dim and codim=0");
      YGridLevelIterator g = begin(level);
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");
      if (cd==0)   // the elements
      {
        if (pitype<=InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->cell_interior.tsubend());
        if (pitype<=All_Partition || pitype == Ghost_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->cell_overlap.tsubend());
      }
      if (cd==dim)   // the vertices
      {
        if (pitype==Interior_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->vertex_interior.tsubend());
        if (pitype==InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->vertex_interiorborder.tsubend());
        if (pitype==Overlap_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->vertex_overlap.tsubend());
        if (pitype<=All_Partition || pitype == Ghost_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g->vertex_overlapfront.tsubend());
      }
      DUNE_THROW(GridError, "YaspLevelIterator with this codim or partition type not implemented");
    }

#if HAVE_MPI
    CollectiveCommunication<MPI_Comm> ccobj;
#else
    CollectiveCommunication<YaspGrid> ccobj;
#endif

    Torus<dim> _torus;

    std::vector< shared_ptr< YaspIndexSet<const YaspGrid<dim>, false > > > indexsets;
    YaspIndexSet<const YaspGrid<dim>, true> leafIndexSet_;
    YaspGlobalIdSet<const YaspGrid<dim> > theglobalidset;

    fTupel _LL;
    iTupel _s;
    std::bitset<dim> _periodic;
    ReservedVector<YGridLevel,32> _levels;
    int _overlap;
    int sizes[32][dim+1]; // total number of entities per level and codim
    bool keep_ovlp;
    int adaptRefCount;
    bool adaptActive;
  };

  //! Output operator for multigrids

  template <int d>
  inline std::ostream& operator<< (std::ostream& s, YaspGrid<d>& grid)
  {
    int rank = grid.torus().rank();

    s << "[" << rank << "]:" << " YaspGrid maxlevel=" << grid.maxLevel() << std::endl;

    for (typename YaspGrid<d>::YGridLevelIterator g=grid.begin(); g!=grid.end(); ++g)
    {
      s << "[" << rank << "]:   " << std::endl;
      s << "[" << rank << "]:   " << "==========================================" << std::endl;
      s << "[" << rank << "]:   " << "level=" << g->level() << std::endl;
      s << "[" << rank << "]:   " << "cell_global=" << g->cell_global << std::endl;
      s << "[" << rank << "]:   " << "cell_overlap=" << g->cell_overlap << std::endl;
      s << "[" << rank << "]:   " << "cell_interior=" << g->cell_interior << std::endl;
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->send_cell_overlap_overlap.begin();
           i!=g->send_cell_overlap_overlap.end(); ++i)
      {
        s << "[" << rank << "]:    " << " s_c_o_o "
          << i->rank << " " << i->grid << std::endl;
      }
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->recv_cell_overlap_overlap.begin();
           i!=g->recv_cell_overlap_overlap.end(); ++i)
      {
        s << "[" << rank << "]:    " << " r_c_o_o "
          << i->rank << " " << i->grid << std::endl;
      }
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->send_cell_interior_overlap.begin();
           i!=g->send_cell_interior_overlap.end(); ++i)
      {
        s << "[" << rank << "]:    " << " s_c_i_o "
          << i->rank << " " << i->grid << std::endl;
      }
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->recv_cell_overlap_interior.begin();
           i!=g->recv_cell_overlap_interior.end(); ++i)
      {
        s << "[" << rank << "]:    " << " r_c_o_i "
          << i->rank << " " << i->grid << std::endl;
      }

      s << "[" << rank << "]:   " << "-----------------------------------------------"  << std::endl;
      s << "[" << rank << "]:   " << "vertex_global="         << g->vertex_global << std::endl;
      s << "[" << rank << "]:   " << "vertex_overlapfront="   << g->vertex_overlapfront << std::endl;
      s << "[" << rank << "]:   " << "vertex_overlap="        << g->vertex_overlap << std::endl;
      s << "[" << rank << "]:   " << "vertex_interiorborder=" << g->vertex_interiorborder << std::endl;
      s << "[" << rank << "]:   " << "vertex_interior="       << g->vertex_interior << std::endl;
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->send_vertex_overlapfront_overlapfront.begin();
           i!=g->send_vertex_overlapfront_overlapfront.end(); ++i)
      {
        s << "[" << rank << "]:    " << " s_v_of_of "
          << i->rank << " " << i->grid << std::endl;
      }
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->recv_vertex_overlapfront_overlapfront.begin();
           i!=g->recv_vertex_overlapfront_overlapfront.end(); ++i)
      {
        s << "[" << rank << "]:    " << " r_v_of_of "
          << i->rank << " " << i->grid << std::endl;
      }
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->send_vertex_overlap_overlapfront.begin();
           i!=g->send_vertex_overlap_overlapfront.end(); ++i)
      {
        s << "[" << rank << "]:    " << " s_v_o_of "
          << i->rank << " " << i->grid << std::endl;
      }
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->recv_vertex_overlapfront_overlap.begin();
           i!=g->recv_vertex_overlapfront_overlap.end(); ++i)
      {
        s << "[" << rank << "]:    " << " r_v_of_o "
          << i->rank << " " << i->grid << std::endl;
      }
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->send_vertex_interiorborder_interiorborder.begin();
           i!=g->send_vertex_interiorborder_interiorborder.end(); ++i)
      {
        s << "[" << rank << "]:    " << " s_v_ib_ib "
          << i->rank << " " << i->grid << std::endl;
      }
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->recv_vertex_interiorborder_interiorborder.begin();
           i!=g->recv_vertex_interiorborder_interiorborder.end(); ++i)
      {
        s << "[" << rank << "]:    " << " r_v_ib_ib "
          << i->rank << " " << i->grid << std::endl;
      }
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->send_vertex_interiorborder_overlapfront.begin();
           i!=g->send_vertex_interiorborder_overlapfront.end(); ++i)
      {
        s << "[" << rank << "]:    " << " s_v_ib_of "
          << i->rank << " " << i->grid << std::endl;
      }
      for (typename std::deque<typename YaspGrid<d>::Intersection>::const_iterator i=g->recv_vertex_overlapfront_interiorborder.begin();
           i!=g->recv_vertex_overlapfront_interiorborder.end(); ++i)
      {
        s << "[" << rank << "]:    " << " s_v_of_ib "
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

    /** \brief YaspGrid has only one geometry type for codim 0 entities
       \ingroup YaspGrid
     */
    template<int dim>
    struct hasSingleGeometryType< YaspGrid<dim> >
    {
      static const bool v = true;
      static const unsigned int topologyId = GenericGeometry :: CubeTopology< dim > :: type :: id ;
    };

    /** \brief YaspGrid is a Cartesian grid
        \ingroup YaspGrid
     */
    template<int dim>
    struct isCartesian< YaspGrid<dim> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid has codim=0 entities (elements)
       \ingroup YaspGrid
     */
    template<int dim>
    struct hasEntity< YaspGrid<dim>, 0 >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid has codim=dim entities (vertices)
       \ingroup YaspGrid
     */
    template<int dim>
    struct hasEntity< YaspGrid<dim>, dim >
    {
      static const bool v = true;
    };

    template< int dim >
    struct canCommunicate< YaspGrid< dim >, 0 >
    {
      static const bool v = true;
    };

    template< int dim >
    struct canCommunicate< YaspGrid< dim >, dim >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is parallel
       \ingroup YaspGrid
     */
    template<int dim>
    struct isParallel< YaspGrid<dim> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is levelwise conforming
       \ingroup YaspGrid
     */
    template<int dim>
    struct isLevelwiseConforming< YaspGrid<dim> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is leafwise conforming
       \ingroup YaspGrid
     */
    template<int dim>
    struct isLeafwiseConforming< YaspGrid<dim> >
    {
      static const bool v = true;
    };

  }

} // end namespace


#endif
