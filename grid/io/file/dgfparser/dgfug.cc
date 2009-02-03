// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  template <int dim>
  inline UGGrid<dim>*
  MacroGrid::Impl<UGGrid<dim> >::generate(MacroGrid& mg,
                                          const char* filename,
                                          MPICommunicatorType )
  {
    mg.element=General;
    std::ifstream gridin(filename);

    std::string str(filename);

    if(mg.readDuneGrid(gridin))
    {
      if (mg.dimw != dim) {
        DUNE_THROW(DGFException,
                   "Macrofile " << filename << " is for dimension " << mg.dimw
                                << " and cannot be used to initialize an "
                                << "UGGrid of dimension "
                                << dim);
      }

      mg.setOrientation(0,1);

      GridFactory<UGGrid<dim> > factory;

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
        else
          DUNE_THROW(DGFException, "Wrong number of vertices for element");
      }

      UGGrid<dim> *resultGrid = factory.createGrid();

      // get grid parameter block
      dgf :: GridParameterBlock gridParam(gridin, false);

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
