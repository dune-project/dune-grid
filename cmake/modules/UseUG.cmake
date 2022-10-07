# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

# If dune-uggrid was found this module adds the dgf magic to config.h
# and makes add_dune_ug_flags available.
#

# Add dgf magic to config.h and register flags
if(dune-uggrid_FOUND)
  dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE UGGRID
    ASSERTION "GRIDDIM == WORLDDIM"
    DUNETYPE "Dune::UGGrid< dimgrid >"
    HEADERS "dune/grid/uggrid.hh" "dune/grid/io/file/dgfparser/dgfug.hh")
endif()

function(add_dune_ug_flags)
  message(DEPRECATION "The function add_dune_ug_flags() is deprecated and will be removed "
                      "after 2.9. It is not necessary to set any additional UG flags and "
                      "to link against the duneuggrid library, since this is done automatically "
                      "when linking the dunegrid library.")
endfunction(add_dune_ug_flags)
