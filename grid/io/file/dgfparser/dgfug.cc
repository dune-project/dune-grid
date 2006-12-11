// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  template <int dim>
  inline UGGrid<dim>*
  MacroGrid::Impl<UGGrid<dim> >::generate(MacroGrid& mg,
                                          const char* filename, MPICommunicatorType ) {
    mg.element=General;
    std::ifstream gridin(filename);

    std::string str(filename);

    if(mg.readDuneGrid(gridin))
    {
      if (mg.dimw != dim) {
        DUNE_THROW(DGFException,
                   "Macrofile " << filename << " is for dimension " << mg.dimw
                                << " and connot be used to initialize an "
                                << "UGGrid of dimension "
                                << dim);
      }
      mg.setOrientation(0,1);
      UGGrid<dim> *grid = new UGGrid<dim>();
      grid->createBegin();
      for (int n=0; n<mg.nofvtx; n++) {
        FieldVector<double,dim> v;
        for (int j=0; j<dim; j++)
          v[j] = mg.vtx[n][j];
        grid->insertVertex(v);
      }
      for (int n=0; n<mg.nofelements; n++) {
        std::vector<unsigned int> el;
        for (size_t j=0; j<mg.elements[n].size(); j++)
          el.push_back((mg.elements[n][j]));
        if ((int) el.size()==mg.dimw+1)
          grid->insertElement(GeometryType(GeometryType::simplex,dim),el);
        else if (el.size()==pow(2,mg.dimw))
          grid->insertElement(GeometryType(GeometryType::cube,dim),el);
        else
          DUNE_THROW(DGFException, "Wrong number of verticies for element");
      }
      grid->createEnd();
      return grid;
    }
    DUNE_THROW(DGFException,"Macrofile " << filename << " is not a valid DGF file\n"
                                         << "No alternative File-Format implemented!");
    return 0;
  }
}
