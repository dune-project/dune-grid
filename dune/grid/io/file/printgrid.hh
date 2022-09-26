// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PRINTGRID_HH
#define DUNE_PRINTGRID_HH

#include <fstream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {

  namespace {

    template<int dim>
    struct ElementDataLayout
    {
      bool contains (Dune::GeometryType gt)
      {
        return gt.dim()==dim;
      }
    };

    template<int dim>
    struct NodeDataLayout
    {
      bool contains (Dune::GeometryType gt)
      {
        return gt.dim()==0;
      }
    };

    // Move a point closer to basegeo's center by factor scale (used for drawing relative to the element)
    template <typename B, typename C>
    C centrify (const B& basegeo, const C& coords, const double scale) {
      C ret = coords;
      ret -= basegeo.center();
      ret *= scale;
      ret += basegeo.center();
      return ret;
    }

    // Add a line to the plotfile from p1 to p2
    template <typename Coord>
    void draw_line (std::ofstream &plotfile, const Coord &p1, const Coord &p2, std::string options) {
      plotfile << "set object poly from ";
      plotfile << p1[0] << "," << p1[1] << " to ";
      plotfile << p2[0] << "," << p2[1] << " to ";
      plotfile << p1[0] << "," << p1[1];
      plotfile << " " << options << std::endl;
    }

  }

  /** \brief Print a grid as a gnuplot for testing and development
   *  \tparam GridType the type of grid to work with
   *  \param grid the grid to print
   *  \param helper an MPIHelper to create unique output file names in parallel case
   *  \param output_file the base of the output filename
   *  \param size size of the plot in pixels; increase if plot is too cramped
   *  \param execute_plot whether to execute gnuplot automatically
   *  \param png whether to use PNG or SVG as the output format
   *  \param local_corner_indices whether to show local corner indices
   *  \param local_intersection_indices whether to show local intersection indices
   *  \param outer_normals whether to show outer normals of intersections
   *  Creates a gnuplot (one per process if parallel) showing the grid structure with indices, intersection types etc.
   */
  template <typename GridType>
  void printGrid (const GridType& grid, const Dune::MPIHelper& helper, std::string output_file = "printgrid",
                  int size = 2000, bool execute_plot = true, bool png = true, bool local_corner_indices = true,
                  bool local_intersection_indices = true, bool outer_normals = true)
  {

    // Create output file
    output_file = output_file + "_" + std::to_string(helper.rank());
    std::string plot_file_name = output_file + ".gnuplot";
    std::ofstream plotfile (plot_file_name, std::ios::out | std::ios::trunc);
    if (!plotfile.is_open()) {
      DUNE_THROW(Dune::IOError, "Could not create plot file " << output_file << "!");
      return;
    }

    // Basic plot settings
    plotfile << "set size ratio -1" << std::endl;
    if (png) {
      plotfile << "set terminal png size " << size << "," << size << std::endl;
      plotfile << "set output '" << output_file << ".png'" << std::endl;
    } else {
      plotfile << "set terminal svg size " << size << "," << size << " enhanced background rgb 'white'" << std::endl;
      plotfile << "set output '" << output_file << ".svg'" << std::endl;
    }

    // Get GridView
    typedef typename GridType::LeafGridView GV;
    const GV gv = grid.leafGridView();

    // Create mappers used to retrieve indices
    typedef typename Dune::MultipleCodimMultipleGeomTypeMapper<GV> Mapper;
    const Mapper elementmapper(gv, mcmgElementLayout());
    const Mapper nodemapper(gv, mcmgVertexLayout());

    // Create iterators
    typedef typename GV::template Codim<0 >::Iterator LeafIterator;
    typedef typename GV::IntersectionIterator IntersectionIterator;

    LeafIterator it = gv.template begin<0>();

    // Will contain min/max coordinates. Needed for scaling of the plot
    Dune::FieldVector<typename GridType::ctype,2> max_coord (it->geometry().center()), min_coord (max_coord);

    // Iterate over elements
    for (; it != gv.template end<0>(); ++it) {

      const auto& entity = *it;
      auto geo = entity.geometry();

      // Plot element index
      int element_id = elementmapper.index(entity);
      plotfile << "set label at " << geo.center()[0] << "," << geo.center()[1] << " '"
            << element_id << "' center" << std::endl;

      for (int i = 0; i < geo.corners(); ++i) {
        // Plot corner indices
        const int globalNodeNumber1 = nodemapper.subIndex(entity, i, 2);
        auto labelpos = centrify (geo, geo.corner(i), 0.7);
        plotfile << "set label at " << labelpos[0] << "," << labelpos[1] << " '" << globalNodeNumber1;
        if (local_corner_indices)
          plotfile << "(" << i << ")";
        plotfile << "' center" << std::endl;

        // Adapt min / max coordinates
        for (int dim = 0; dim < 2; ++dim) {
          if (geo.corner(i)[dim] < min_coord[dim])
            min_coord[dim] = geo.corner(i)[dim];
          else if (geo.corner(i)[dim] > max_coord[dim])
            max_coord[dim] = geo.corner(i)[dim];
        }
      }

      // Iterate over intersections
      for (IntersectionIterator is = gv.ibegin(entity); is != gv.iend(entity); ++is) {

        const auto& intersection = *is;
        auto igeo = intersection.geometry();

        // Draw intersection line
        draw_line (plotfile, igeo.corner(0), igeo.corner(1), "fs empty border 1");

        // Plot local intersection index
        if (local_intersection_indices) {
          auto label_pos = centrify (geo, igeo.center(), 0.8);
          plotfile << "set label at " << label_pos[0] << "," << label_pos[1]
                << " '" << intersection.indexInInside() << "' center" << std::endl;
        }

        // Plot outer normal
        if (outer_normals) {
          auto intersection_pos = igeo.center();
          auto normal = intersection.centerUnitOuterNormal();
          normal *= 0.15 * igeo.volume();
          auto normal_end = intersection_pos + normal;
          plotfile << "set arrow from " << intersection_pos[0] << "," << intersection_pos[1]
                << " to " << normal_end[0] << "," << normal_end[1] << " lt rgb \"gray\"" << std::endl;
        }

        // Get corners for inner intersection representation
        auto inner_corner1 = centrify (geo, igeo.corner(0), 0.5);
        auto inner_corner2 = centrify (geo, igeo.corner(1), 0.5);

        // Thick line in case of boundary()
        if (intersection.boundary())
          draw_line (plotfile, inner_corner1, inner_corner2, "fs empty border 3 lw 4");

        // Thin line with color according to neighbor()
        if (intersection.neighbor())
          draw_line (plotfile, inner_corner1, inner_corner2, "fs empty border 2");
        else
          draw_line (plotfile, inner_corner1, inner_corner2, "fs empty border 1");
      }

    }

    // Finish plot, pass extend of the grid
    Dune::FieldVector<typename GridType::ctype,2> extend (max_coord - min_coord);

    extend *= 0.2;
    min_coord -= extend;
    max_coord += extend;
    plotfile << "plot [" << min_coord[0] << ":" << max_coord[0] << "] [" << min_coord[1]
          << ":" << max_coord[1] << "] NaN notitle" << std::endl;
    plotfile.close();

    if (execute_plot) {
      std::string cmd = "gnuplot -p '" + plot_file_name + "'";
      if (std::system (cmd.c_str()) != 0)
        DUNE_THROW(Dune::Exception,"Error running GNUPlot: " << cmd);
    }
  }

}

#endif // #ifndef DUNE_PRINTGRID_HH
