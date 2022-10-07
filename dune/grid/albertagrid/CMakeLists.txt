# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(test)

if(Alberta_FOUND)
  foreach(_dim ${ALBERTA_WORLD_DIMS})
    target_sources(dunealbertagrid${_dim}d PRIVATE
      dgfparser.cc
      elementinfo.cc
      indexsets.cc
      macrodata.cc
      meshpointer.cc)
  endforeach(_dim)
endif()

#install header
set(HEADERS
  agrid.hh
  albertagrid.cc
  albertaheader.hh
  indexsets.hh
  indexstack.hh
  datahandle.hh
  misc.hh
  macroelement.hh
  elementinfo.hh
  geometrycache.hh
  meshpointer.hh
  macrodata.hh
  dofadmin.hh
  dofvector.hh
  refinement.hh
  coordcache.hh
  level.hh
  undefine-2.0.hh
  undefine-3.0.hh
  entity.hh
  entity.cc
  entityseed.hh
  hierarchiciterator.hh
  algebra.hh
  geometry.hh
  geometry.cc
  projection.hh
  transformation.hh
  leveliterator.hh
  leafiterator.hh
  treeiterator.hh
  intersection.hh
  intersection.cc
  intersectioniterator.hh
  capabilities.hh
  gridfactory.hh
  dgfparser.hh
  albertareader.hh
  gridfamily.hh
  gridview.hh
  persistentcontainer.hh
  backuprestore.hh
  structuredgridfactory.hh
  geometryreference.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/albertagrid/)
