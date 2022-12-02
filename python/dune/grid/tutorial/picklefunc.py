# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import numpy
def globalF(x):
    x -= [0.1,0.1,0.1]
    return numpy.sin( (x[0]*x[1]*x[2])*numpy.pi )
def localF(e,y):
    x = e.geometry.toGlobal(y)
    x -= [0.1,0.1,0.1]
    return numpy.sin( (x[0]*x[1]*x[2])*numpy.pi )

class TimeDependent:
    def __init__(self):
        self.t = 0
    def __call__(self,e,y):
        x = e.geometry.toGlobal(y)
        x -= [0.1,0.1,0.1]
        return numpy.sin( (x[0]*x[1]+self.t*x[1])*numpy.pi )
