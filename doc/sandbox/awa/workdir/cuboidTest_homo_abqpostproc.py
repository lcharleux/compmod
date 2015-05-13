# ABQPOSTPROC.PY
# Warning: executable only in abaqus abaqus viewer -noGUI,... not regular python.
import sys
from abapy.postproc import GetFieldOutput_byRpt as gfo
from abapy.postproc import GetVectorFieldOutput_byRpt as gvfo
from abapy.postproc import GetTensorFieldOutput_byRpt as gtfo
from abapy.postproc import GetHistoryOutputByKey as gho
from abapy.postproc import GetMesh
from abapy.indentation import Get_ContactData
from abapy.misc import dump
from odbAccess import openOdb
from abaqusConstants import JOB_STATUS_COMPLETED_SUCCESSFULLY



# Odb opening  
file_name = 'cuboidTest_homo'
odb = openOdb(file_name + '.odb')
data = {}

# Check job status:
job_status = odb.diagnosticData.jobStatus

if job_status == JOB_STATUS_COMPLETED_SUCCESSFULLY:
  data['completed'] = True# History Outputs
  data['history'] = {} 
  ho = data['history']
  ho['disp'] =   gho(odb,'U2')
  ho['force'] =   gho(odb,'RF2')
  ho['allse'] =   gho(odb,'ALLSE').values()[0]
  ho['allpd'] =   gho(odb,'ALLPD').values()[0]
  ho['allwk'] =   gho(odb,'ALLWK').values()[0]
  ho['volume'] =  gho(odb,'EVOL')
  
  # Mesh 
  data['mesh'] = GetMesh(odb, "ISAMPLE")
  
else:
  data['completed'] = False
# Closing and dumping
odb.close()
dump(data, file_name+'.pckl')