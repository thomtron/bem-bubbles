import eel
import os
import glob

eel.init("interface")

@eel.expose
def find_paths(path):
    return glob.glob('interface/models/*.ply')



#def close_callback(route, websockets):
#    if not websockets:
#        exit()

eel.start("index.html")

exit()