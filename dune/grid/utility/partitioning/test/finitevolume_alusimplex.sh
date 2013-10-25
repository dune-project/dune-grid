#!/bin/sh
set -ex
rm -f alusimplex_seq_plain_concentration-00001.vtu      \
      alusimplex_seq_opt_concentration-00001.vtu        \
      alusimplex_tbb_plain_concentration-00001.vtu      \
      alusimplex_tbb_opt_concentration-00001.vtu
timeout 20 ./finitevolume_alusimplex "${REFINES:-6}"
cmp alusimplex_seq_plain_concentration-00001.vtu        \
    alusimplex_seq_opt_concentration-00001.vtu
cmp alusimplex_seq_plain_concentration-00001.vtu        \
    alusimplex_tbb_plain_concentration-00001.vtu
cmp alusimplex_seq_opt_concentration-00001.vtu  \
    alusimplex_tbb_opt_concentration-00001.vtu
