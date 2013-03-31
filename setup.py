from distutils.core import setup
import py2exe

#setup(windows=['main.py'])
includes = ["encodings", "encodings.*"]    

options = {"py2exe":

    {"compressed": 1,
     "optimize": 2,
     "ascii": 1,
     "includes":includes,
     "bundle_files": 1
     }
    }
setup(     
    options = options,      
    zipfile=None,   
    windows=[{"script": "main.py", "icon_resources": [(1, "xvx.ico")] }]
    )
