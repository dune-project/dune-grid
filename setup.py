from skbuild import setup

from dune.dunepackaging import metaData

print(metaData()[1])

setup(**metaData()[1])
