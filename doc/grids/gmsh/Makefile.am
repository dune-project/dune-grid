GRIDS = circle1storder.msh  curved2d.msh         pyramid4.msh \
telescope2ndorder.msh circle2ndorder.msh        pyramid.geo  \
telescope.geo circle.geo          pyramid1storder.msh  pyramid.msh \
telescope.msh curved2d.geo        pyramid2ndorder.msh  telescope1storder.msh \
sphere.msh  oned-testgrid.msh

gmshdir=$(docdir)/grids/gmsh
gmsh_DATA = $(GRIDS)
EXTRA_DIST = CMakeLists.txt $(GRIDS)

include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/global-rules
