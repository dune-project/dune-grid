// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  template <int dim,int dimworld>
  inline UGGrid<dim,dimworld>*
  MacroGrid::Impl<UGGrid<dim,dimworld> >::generate(MacroGrid& mg,
                                                   const char* filename, MPICommunicatorType ) {
    mg.element=General;
    std::ifstream gridin(filename);

    std::string str(filename);

    if(mg.readDuneGrid(gridin))
    {
      if (mg.dimw != dimworld) {
        DUNE_THROW(DGFException,
                   "Macrofile " << filename << " is for dimension " << mg.dimw
                                << " and connot be used to initialize an "
                                << "UGGrid of dimension "
                                << dimworld);
      }
      mg.setOrientation(0);
      UGGrid<dim,dimworld> *grid = new UGGrid<dim,dimworld>();
      grid->createBegin();
      for (int n=0; n<mg.nofvtx; n++) {
        FieldVector<double,dimworld> v;
        for (int j=0; j<dimworld; j++)
          v[j] = mg.vtx[n][j];
        grid->insertVertex(v);
      }
      for (int n=0; n<mg.nofelements; n++) {
        std::vector<unsigned int> el;
        for (size_t j=0; j<mg.elements[n].size(); j++)
          el.push_back((mg.elements[n][j]));
        if ((int) el.size()==mg.dimw+1)
          grid->insertElement(GeometryType(GeometryType::simplex,dimworld),el);
        else if (el.size()==pow(2,mg.dimw))
          grid->insertElement(GeometryType(GeometryType::cube,dimworld),el);
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
