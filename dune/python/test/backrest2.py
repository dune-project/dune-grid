# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import pickle
import dune.generator # needed to set path to 'dune.generated' - could be moved to 'import dune'
[a,string,otherHGrid,value] = pickle.load(open("dumpA","rb"))

otherGrid = otherHGrid.leafView

otherGrid.plot()

print("leaf after refine", otherGrid.size(0))
print("level 1 after refine",
      otherGrid.hierarchicalGrid.levelView(1).size(0))
otherGrid.hierarchicalGrid.globalRefine(-2)
print("coarsen other", otherGrid.size(0))
print("numpy vector",a)

test = pickle.load(open("dumpB","rb"))
print("test:",test.run())
