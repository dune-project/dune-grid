
if(LATEX_USABLE)
  set(EPS_IMGS simplexvertexindex.eps referencetokuhn0.eps kuhntriangulation.eps)
  ADD_LATEX_DOCUMENT(refinement.tex FATHER_TARGET doc
    DEFAULT_SAFEPDF INPUTS poosc08_test.cc
    IMAGES ${EPS_IMGS})
  set(PDFDOC refinement.pdf)
endif(LATEX_USABLE)

set(DOCS  referencetokuhn0.png simplexvertexindex.png)

if(INKSCAPE_FOUND)
  set(DOCS ${DOCS} referencetokuhn0.png simplexvertexindex.png)
  inkscape_generate_eps_from_png(DPI 180 OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR} ${EPS_IMGS})
endif(INKSCAPE_FOUND)

install(FILES ${PDFDOCS} DESTINATION ${CMAKE_INSTALL_DOCDIR}/refinement)
