import importlib.util

if importlib.util.find_spec('tkinter') is not None:
    from .exporter import *

from .fitfuncs import *

if importlib.util.find_spec('IPython') is not None:
    from .jupyter_helper import *

from .process import *