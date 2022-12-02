# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

from dune.grid import gridFunction
from picklefunc import localF, TimeDependent
func = TimeDependent()

def error(gv,t,df,dfs):
    if t is None:
        return [ gridFunction(gv, name="error", order=3)(
                   lambda e,x: abs(df(e,x)-localF(e,x))
               )]
    else:
        func.t = t
        return [ gridFunction(gv, name="error", order=3)(
                   lambda e,x: abs(df(e,x)-func(e,x))
               )]
register = [error]
