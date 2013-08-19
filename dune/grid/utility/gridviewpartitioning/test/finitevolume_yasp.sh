#!/bin/sh
set -ex
rm -f yasp_seq_plain_concentration-00001.vtu   \
      yasp_seq_opt_concentration-00001.vtu     \
      yasp_ost_plain_concentration-00001.vtu   \
      yasp_ost_opt_concentration-00001.vtu
./finitevolume_yasp "${REFINES:-7}"
cmp yasp_seq_plain_concentration-00001.vtu     \
    yasp_seq_opt_concentration-00001.vtu
cmp yasp_seq_plain_concentration-00001.vtu     \
    yasp_ost_plain_concentration-00001.vtu
cmp yasp_seq_opt_concentration-00001.vtu  \
    yasp_ost_opt_concentration-00001.vtu
