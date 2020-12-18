import os

path = os.path.join( os.path.dirname(__file__), "tutorial" )
execute  = "cp -r " + path + " "
execute += "grid_tutorial"
status = os.system(execute)
if status != 0: raise RuntimeError(status)

print("##################################################################")
print("## An example script is now located in the 'grid_tutorial' folder.")
try:
    import matplotlib
except ImportError:
    print("## Note: the examples requires the installation of 'matplotlib'.")
print("##################################################################")
