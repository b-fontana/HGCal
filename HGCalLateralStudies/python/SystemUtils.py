import errno
import os

def createDir(name):
    if not os.path.isdir(name):
        os.mkdir(name)
        
