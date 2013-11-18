#
# Dune's experimental grid extensions are enabled by adding
# -DDUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS=TRUE to the
# CMAKE_FLAGS
# This file adds the summary line whether it is set.
#

include(FeatureSummary)
add_feature_info("Experimental grid extensions" DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS "Enables additional grid features.")
