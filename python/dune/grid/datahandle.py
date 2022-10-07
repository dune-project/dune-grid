# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

from dune.generator.generator import SimpleGenerator

from dune.common.hashit import hashIt

generator = SimpleGenerator("DataHandle", "Dune::Python")

def load(includes, typeName, *args):
    includes = includes + ["dune/python/common/numpycommdatahandle.hh"]
    moduleName = "numpycommdatahandle_" + hashIt(typeName)
    return generator.load(includes, typeName, moduleName, *args)

def dataHandle(mapper,array,function):
    typeName = "Dune::Python::NumPyCommDataHandle< " + mapper.cppTypeName + ", double, std::function<double(double,double)> >"
    includes = mapper.cppIncludes
    return load(includes, typeName).DataHandle(mapper,array,function)
