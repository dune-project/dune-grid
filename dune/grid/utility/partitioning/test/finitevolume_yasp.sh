#!/bin/sh
set -ex
rm -f yasp_seq_plain_concentration-00001.vtu   \
      yasp_seq_opt_concentration-00001.vtu     \
      yasp_tbb_plain_concentration-00001.vtu   \
      yasp_tbb_opt_concentration-00001.vtu
timeout 20 ./finitevolume_yasp "${REFINES:-7}"
cmp yasp_seq_plain_concentration-00001.vtu     \
    yasp_seq_opt_concentration-00001.vtu
cmp yasp_seq_plain_concentration-00001.vtu      \
    yasp_tbb_plain_concentration-00001.vtu
cmp yasp_seq_opt_concentration-00001.vtu        \
    yasp_tbb_opt_concentration-00001.vtu
