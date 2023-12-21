import eel
import os
import glob

eel.init("interface")

@eel.expose
def find_paths(path):
    return glob.glob('interface/models/*.ply') 
    # path to models (create a symbolic link at interface/models pointing to the directory where the ply's are stored)



#def close_callback(route, websockets):
#    if not websockets:
#        exit()

eel.start("index.html")

exit()