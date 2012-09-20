
if(LATEX_USABLE AND INKSCAPE)
  set(DOCS ${DOCS} referencetokuhn0.png simplexvertexindex.png)
  set(EPS_IMGS simplexvertexindex.eps referencetokuhn0.eps kuhntriangulation.eps)
  inkscape_generate_eps_from_svg(DPI 180
    OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
    INPUT_DIR ${CMAKE_CURRENT_BINARY_DIR} ${EPS_IMGS})
  ADD_LATEX_DOCUMENT(refinement.tex FATHER_TARGET doc
    DEFAULT_SAFEPDF
    DEPENDS ${EPS_IMGS})
  set(PDFDOCS ${CMAKE_CURRENT_BINARY_DIR}/refinement.pdf)
endif(LATEX_USABLE AND INKSCAPE)

set(DOCS  referencetokuhn0.png simplexvertexindex.png)

if(INKSCAPE)

endif(INKSCAPE)
create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/refinement.pdf ${CMAKE_INSTALL_DOCDIR}/refinement)
