// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \cond */
namespace Dune {
  template <int dim>
  inline UGGrid<dim>*
  MacroGrid::Impl<UGGrid<dim> >::generate(MacroGrid& mg,
                                          const char* filename,
                                          MPICommunicatorType )
  {
    mg.element=General;
    std::ifstream gridin(filename);
    if (! gridin)
    {
      DUNE_THROW(DGFException,
                 "Macrofile " << filename << " not found");
    }

    std::string str(filename);

    if( mg.readDuneGrid( gridin, dim, dim ) )
    {
      mg.setOrientation(0,1);

      // get grid parameter block
      dgf :: UGGridParameterBlock gridParam(gridin);

      // create grid here to set heap size
      // create grid factory (passed grid is returned by createGrid method)
      GridFactory<UGGrid<dim> > factory( new UGGrid<dim> ( gridParam.heapSize() ) );

      for (int n=0; n<mg.nofvtx; n++)
      {
        FieldVector<double,dim> v;
        for (int j=0; j<dim; j++)
          v[j] = mg.vtx[n][j];

        factory.insertVertex(v);
      }

      // eval 2^dim
      size_t two_power_dim = 1;
      for(int i=0; i<mg.dimw; ++i) two_power_dim *= 2;

      for (int n=0; n<mg.nofelements; n++)
      {
        std::vector<unsigned int> el;
        for (size_t j=0; j<mg.elements[n].size(); j++)
        {
          el.push_back((mg.elements[n][j]));
        }

        // simplices
        if ((int) el.size()== mg.dimw+1)
          factory.insertElement(GeometryType(GeometryType::simplex,dim),el);
        // cubes
        else if (el.size() == two_power_dim)
          factory.insertElement(GeometryType(GeometryType::cube,dim),el);
#ifdef EXPERIMENTAL_GRID_EXTENSIONS
        // pyramid
        else if (el.size() == 5 )
          factory.insertElement(GeometryType(GeometryType::pyramid,dim),el);
        // prisms
        else if (el.size() == 6 )
          factory.insertElement(GeometryType(GeometryType::prism,dim),el);
#endif
        else
          DUNE_THROW(DGFException, "Wrong number of vertices for element");
      }

      UGGrid<dim> *resultGrid = factory.createGrid();

      // set closure type to none if parameter say so
      if( gridParam.noClosure() )
      {
        resultGrid->setClosureType(UGGrid<dim> :: NONE);
      }
      if ( !gridParam.noCopy() )
      {
        resultGrid->setRefinementType(UGGrid<dim>::COPY);
      }
      return resultGrid;
    }
    DUNE_THROW(DGFException,"Macrofile " << filename << " is not a valid DGF file\n"
                                         << "No alternative File-Format implemented!");
    return 0;
  }
}
/** \endcond */