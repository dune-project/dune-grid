  if(INKSCAPE_FOUND)
    inkscape_generate_eps_from_svg(DPI 180 OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR} ${EPS_IMGS})
  endif(INKSCAPE_FOUND)

if(LATEX_USABLE)
  set(EPS_IMGS simplexvertexindex.eps referencetokuhn0.eps kuhntriangulation.eps)
  ADD_LATEX_DOCUMENT(refinement.tex FATHER_TARGET doc
    DEFAULT_SAFEPDF INPUTS poosc08_test.cc
    IMAGES ${EPS_IMGS})
  set(PDFDOC refinement.pdf)
endif(LATEX_USABLE)

create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/refinement.pdf ${CMAKE_INSTALL_DOCDIR}/refinement)
