# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.9-EF1 replay file
# Internal Version: 2009_10_21-15.49.55 96721
# Run by lcharleux on Mon Nov 17 16:02:28 2014
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=102.02082824707, 
    height=109.766662597656)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
