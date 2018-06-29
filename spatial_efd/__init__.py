import sys
if(sys.version[0] == 2):
	from spatial_efd import *
else:
	from .spatial_efd import *
