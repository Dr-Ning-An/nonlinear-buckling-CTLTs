###################################################
### Created By Ning An, anning003@stu.xjtu.edu.cn
### 2021/1
### Created in Abaqus Version 2020
###################################################
# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import os


#############################################################################################################################
# Created by Ning An
# 2019/12
# http://www.anning003.com/create-virtual-nodes
# Created in Abaqus Version 2017

# Function for creating virtual nodes
# mdb: model database
# NameModel: A string with the name of your model
# NameRef1 & NameRef2: A string with the name of your two virtual nodes.
#############################################################################################################################
def VirtualNodes(mdb, NameModel, NameRef1, x, y, z):
    from part import THREE_D, DEFORMABLE_BODY
    #Create reference parts and assemble
    mdb.models[NameModel].Part(dimensionality=THREE_D, name=NameRef1, type=
        DEFORMABLE_BODY)
    mdb.models[NameModel].parts[NameRef1].ReferencePoint(point=(x,y,z))
    mdb.models[NameModel].rootAssembly.Instance(dependent=ON, name=NameRef1, 
        part=mdb.models[NameModel].parts[NameRef1])


    #Create set of reference points
    mdb.models[NameModel].rootAssembly.Set(name=NameRef1, referencePoints=(
        mdb.models[NameModel].rootAssembly.instances[NameRef1].referencePoints[1],))


Mdb()

# execfile("Parameters_for_NonlinearBucklingCTLTs_X.py")
h_total = 8.0 + 10.6*105.0/180.0*pi
meshSize = h_total/50.0

l = 504.0 # Length
w = 5.5 # Web
L = 3.5 # lumbus
r = (h_total - 2.0*w - L)*3.0/4.0/pi  # Radius



session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

LBModel = 'LinearBucklingModel'
mdb.models.changeKey(fromName='Model-1', toName=LBModel)

mdb.models[LBModel].ConstrainedSketch(name='__profile__', sheetSize=400.0)
mdb.models[LBModel].sketches['__profile__'].Line(
    point1=(0.0, 0.0), 
    point2=(w, 0.0))
mdb.models[LBModel].sketches['__profile__'].ArcByCenterEnds(
    center=(w, r), 
    direction=COUNTERCLOCKWISE, 
    point1=(w, 0.0), 
    point2=(w + r*sin(60.0*pi/180.0), r*(1-cos(60.0*pi/180.0))))
mdb.models[LBModel].sketches['__profile__'].ArcByCenterEnds(
    center=(w + r*sin(60.0*pi/180.0)*2.0, 0.0), 
    direction=CLOCKWISE, 
    point1=(w + r*sin(60.0*pi/180.0), r*(1-cos(60.0*pi/180.0))), 
    point2=(w + r*sin(60.0*pi/180.0)*2.0, r))
mdb.models[LBModel].sketches['__profile__'].Line(
    point1=(w + r*sin(60.0*pi/180.0)*2.0, r), 
    point2=(w + r*sin(60.0*pi/180.0)*2.0 + L/2.0, r))

mdb.models[LBModel].sketches['__profile__'].ConstructionLine(
	point1=(w + r*sin(60.0*pi/180.0)*2.0 + L/2.0, r), point2=(w + r*sin(60.0*pi/180.0)*2.0 + L/2.0, r/2.0))
mdb.models[LBModel].sketches['__profile__'].copyMirror(
    mirrorLine=mdb.models[LBModel].sketches['__profile__'].geometry.findAt((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0, 15.0), ), 
    objectList=(
    mdb.models[LBModel].sketches['__profile__'].geometry.findAt((w/2.0, 0.0), ), 
    mdb.models[LBModel].sketches['__profile__'].geometry.findAt((w + r*sin(30.0*pi/180.0), r*(1-cos(30.0*pi/180.0))), ), 
    mdb.models[LBModel].sketches['__profile__'].geometry.findAt((w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0), r*sin(60.0*pi/180.0)), ), 
    mdb.models[LBModel].sketches['__profile__'].geometry.findAt((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0, r), )))

mdb.models[LBModel].sketches['__profile__'].ConstructionLine(
    point1=(0.0, 0.0), point2=(10.0, 0.0))
mdb.models[LBModel].sketches['__profile__'].copyMirror(
    mirrorLine=mdb.models[LBModel].sketches['__profile__'].geometry.findAt((10.0, 0.0), ), 
    objectList=(
    mdb.models[LBModel].sketches['__profile__'].geometry.findAt((w + r*sin(30.0*pi/180.0), r*(1-cos(30.0*pi/180.0))), ), 
    mdb.models[LBModel].sketches['__profile__'].geometry.findAt((w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0), r*sin(60.0*pi/180.0)), ), 
    mdb.models[LBModel].sketches['__profile__'].geometry.findAt((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0, r), ), 
    mdb.models[LBModel].sketches['__profile__'].geometry.findAt(((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(30.0*pi/180.0)), r*(1-cos(30.0*pi/180.0))), ), 
    mdb.models[LBModel].sketches['__profile__'].geometry.findAt(((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0)), r*sin(60.0*pi/180.0)), ), 
    mdb.models[LBModel].sketches['__profile__'].geometry.findAt(((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0), r), )))


mdb.models[LBModel].ConstrainedSketch(name='Sketch-1', objectToCopy=
    mdb.models[LBModel].sketches['__profile__'])

mdb.models[LBModel].Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models[LBModel].parts['Part-1'].BaseShellExtrude(depth=l, sketch=
    mdb.models[LBModel].sketches['__profile__'])


mdb.models[LBModel].parts['Part-1'].Set(faces=
    mdb.models[LBModel].parts['Part-1'].faces.findAt(((w/2.0, 0.0, 0.0), ), 
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - w/2.0, 0.0, 0.0), ), 
                                                       ), name='Set-Web')
mdb.models[LBModel].parts['Part-1'].Set(faces=
    mdb.models[LBModel].parts['Part-1'].faces.findAt(((w + r*sin(30.0*pi/180.0), r*(1-cos(30.0*pi/180.0)), 0.0), ), 
                                                       ((w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0), r*sin(60.0*pi/180.0), 0.0), ), 
                                                       ((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0, r, 0.0), ), 
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(30.0*pi/180.0)), r*(1-cos(30.0*pi/180.0)), 0.0), ),
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0)), r*sin(60.0*pi/180.0), 0.0), ),
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0), r, 0.0), ), 
                                                       ((w + r*sin(30.0*pi/180.0), -r*(1-cos(30.0*pi/180.0)), 0.0), ), 
                                                       ((w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0), -r*sin(60.0*pi/180.0), 0.0), ), 
                                                       ((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0, -r, 0.0), ), 
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(30.0*pi/180.0)), -r*(1-cos(30.0*pi/180.0)), 0.0), ),
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0)), -r*sin(60.0*pi/180.0), 0.0), ),
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0), -r, 0.0), ), 
                                                       ), name='Set-Flanges')

mdb.models[LBModel].Material(name='Material-1')
mdb.models[LBModel].materials['Material-1'].Elastic(
    table=((128.0e3, 6.5e3, 0.35, 7.5e3, 7.5e3, 7.5e3), ), type=LAMINA)

mdb.models[LBModel].parts['Part-1'].MaterialOrientation(
    additionalRotationType=ROTATION_NONE, axis=AXIS_2, fieldName='', localCsys=
    None, orientationType=GLOBAL, region=
    mdb.models[LBModel].parts['Part-1'].sets['Set-Flanges'])
mdb.models[LBModel].parts['Part-1'].MaterialOrientation(
    additionalRotationType=ROTATION_NONE, axis=AXIS_2, fieldName='', localCsys=
    None, orientationType=GLOBAL, region=
    mdb.models[LBModel].parts['Part-1'].sets['Set-Web'])

t = 0.071
mdb.models[LBModel].CompositeShellSection(
    idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, 
    layup=(SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1'), 
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'),  
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'), 
           SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1')), 
    name='Section-Flanges', 
    poissonDefinition=DEFAULT, preIntegrate=OFF, symmetric=False, temperature=GRADIENT, thicknessModulus=
    None, thicknessType=UNIFORM, useDensity=OFF)

mdb.models[LBModel].CompositeShellSection(
    idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, 
    layup=(SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1'), 
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'),  
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'), 
           SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1'),
           SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1'),
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'),
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'),
           SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1')), 
    name='Section-Web', 
    poissonDefinition=DEFAULT, preIntegrate=OFF, symmetric=False, temperature=GRADIENT, thicknessModulus=
    None, thicknessType=UNIFORM, useDensity=OFF)


mdb.models[LBModel].parts['Part-1'].SectionAssignment(
    offset=0.0, 
    offsetField='', 
    offsetType=MIDDLE_SURFACE, 
    region=mdb.models[LBModel].parts['Part-1'].sets['Set-Web'], 
    sectionName='Section-Web', 
    thicknessAssignment=FROM_SECTION)


mdb.models[LBModel].parts['Part-1'].SectionAssignment(
    offset=0.0, 
    offsetField='', 
    offsetType=BOTTOM_SURFACE, 
    region=mdb.models[LBModel].parts['Part-1'].sets['Set-Flanges'], 
    sectionName='Section-Flanges', 
    thicknessAssignment=FROM_SECTION)

mdb.models[LBModel].rootAssembly.Instance(dependent=OFF, name='Part-1-1', 
    part=mdb.models[LBModel].parts['Part-1'])

mdb.models[LBModel].rootAssembly.seedPartInstance(deviationFactor=0.1, 
    minSizeFactor=0.1, regions=(
    mdb.models[LBModel].rootAssembly.instances['Part-1-1'], ), size=meshSize)
mdb.models[LBModel].rootAssembly.generateMesh(regions=(
    mdb.models[LBModel].rootAssembly.instances['Part-1-1'], ))

Bottom_Nodes = mdb.models[LBModel].rootAssembly.instances['Part-1-1'].nodes.getByBoundingBox(-200.0, -200.0, -0.0001, 200.0, 200.0, 0.0001) 
mdb.models[LBModel].rootAssembly.Set(name='Bottom_Nodes', nodes=Bottom_Nodes)
Top_Nodes = mdb.models[LBModel].rootAssembly.instances['Part-1-1'].nodes.getByBoundingBox(-200.0, -200.0, l-0.0001, 200.0, 200.0, l+0.0001) 
mdb.models[LBModel].rootAssembly.Set(name='Top_Nodes', nodes=Top_Nodes)

# # # #
VirtualNodes(mdb, LBModel, 'Ref-1', (w + r*sin(60.0*pi/180.0)*2.0 + L/2.0), 0.0, 0.0)
VirtualNodes(mdb, LBModel, 'Ref-2', (w + r*sin(60.0*pi/180.0)*2.0 + L/2.0), 0.0, l)
mdb.models[LBModel].Coupling(controlPoint=
    mdb.models[LBModel].rootAssembly.sets['Ref-1'], couplingType=KINEMATIC, 
    influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-1', 
    surface=mdb.models[LBModel].rootAssembly.sets['Bottom_Nodes'], u1=ON, u2=ON, 
    u3=ON, ur1=ON, ur2=ON, ur3=ON)
mdb.models[LBModel].Coupling(controlPoint=
    mdb.models[LBModel].rootAssembly.sets['Ref-2'], couplingType=KINEMATIC, 
    influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-2', 
    surface=mdb.models[LBModel].rootAssembly.sets['Top_Nodes'], u1=ON, u2=ON, 
    u3=ON, ur1=ON, ur2=ON, ur3=ON)

mdb.models[LBModel].BuckleStep(blockSize=DEFAULT, eigensolver=LANCZOS, 
    maxBlocks=DEFAULT, minEigen=0.0, name='LinearBucklingStep', numEigen=20, previous=
    'Initial')
# # # #
mdb.models[LBModel].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
    region=mdb.models[LBModel].rootAssembly.sets['Ref-1'], u1=SET, u2=SET, 
    u3=SET, ur1=SET, ur2=SET, ur3=SET)
mdb.models[LBModel].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', 
    region=mdb.models[LBModel].rootAssembly.sets['Ref-2'], u1=SET, u2=UNSET, 
    u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET)

mdb.models[LBModel].Moment(cm1=-1000.0, createStepName='LinearBucklingStep', 
    distributionType=UNIFORM, field='', localCsys=None, name='Load-1', region=
    mdb.models[LBModel].rootAssembly.sets['Ref-2'])

# if not os.path.exists('FEModelFiles'):
#     os.mkdir('FEModelFiles')
# os.chdir('FEModelFiles')
ParName = str(floor(w)).split('.')[0] + '_' + \
              str(float("%.2f" % floor(w*10-10*floor(w)))).split('.')[0] + \
              str(float("%.2f" % abs(w*100-100*floor(w) - 10*floor(w*10-10*floor(w))))).split('.')[0] + '+' +\
              str(floor(L)).split('.')[0] + '_' + \
              str(float("%.2f" % floor(L * 10 - 10 * floor(L)))).split('.')[0] + \
              str(float("%.2f" % abs(L * 100 - 100 * floor(L) - 10 * floor(L * 10 - 10 * floor(L))))).split('.')[0]
jobName_Linear = 'CTLTs_LinearBuckling_X-' + ParName
mdb.Job(model=LBModel, name=jobName_Linear, numCpus=12, numDomains=12, numGPUs=1)
mdb.saveAs(pathName=jobName_Linear)
mdb.jobs[jobName_Linear].submit()
mdb.jobs[jobName_Linear].waitForCompletion()


##############################################################
## Extract EigenValues from Linear Buckling Analysis
## http://www.anning003.com/extract-eigenvalues/
##############################################################

stepName_Linear = 'LinearBucklingStep'

from odbAccess import*
from abaqusConstants import*
import string
import numpy as np
import os
import re

odb_Linear = openOdb(path = jobName_Linear + '.odb')
outfile = open(jobName_Linear + '.csv', 'w')
outfile.write('Mode Order'+ ',' + 'Linear Critical Load [N m]'+ '\n')
for fm in range(1, len(odb_Linear.steps[stepName_Linear].frames)):
    text = odb_Linear.steps[stepName_Linear].frames[fm].description
    numeric_const_pattern = '[-+]? (?: (?: \\d* \\. \\d+ ) | (?: \\d+ \\.? ) )(?: [Ee] [+-]? \\d+ ) ?'
    rx = re.compile(numeric_const_pattern, re.VERBOSE)
    EigenValue = [float(i) for i in rx.findall(text)]
    outfile.write(str(int(EigenValue[0])) + ',' + str(EigenValue[1]) + '\n')
outfile.close()
odb_Linear.close()
print("=================================================")
print("Linear Buckling Analysis Completed Successfully!!")
print("=================================================")

#=========================================================================
# Delete unnecessary steps
#=========================================================================
#---------------------------------------------------------------
NLBModel = 'NonLinearBucklingModel'
mdb.Model(name=NLBModel, objectToCopy=mdb.models[LBModel])
del mdb.models[NLBModel].steps['LinearBucklingStep']

mdb.models[NLBModel].StaticStep(initialInc=0.01, maxNumInc=1000, name=
    'AppliedPreLoad', nlgeom=ON, previous='Initial', maxInc=0.01)
mdb.models['NonLinearBucklingModel'].steps['AppliedPreLoad'].Restart(frequency=
    1, numberIntervals=0, overlay=OFF, timeMarks=OFF)
mdb.models[NLBModel].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'U', 'RM'))
mdb.models[NLBModel].BuckleStep(blockSize=DEFAULT, eigensolver=LANCZOS, 
    maxBlocks=DEFAULT, minEigen=0.0, name='NonLinearBucklingStep', numEigen=20, previous=
    'AppliedPreLoad')

# # # #
mdb.models[NLBModel].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
    region=mdb.models[NLBModel].rootAssembly.sets['Ref-1'], u1=SET, u2=SET, 
    u3=SET, ur1=SET, ur2=SET, ur3=SET)
mdb.models[NLBModel].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', 
    region=mdb.models[NLBModel].rootAssembly.sets['Ref-2'], u1=SET, u2=UNSET, 
    u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET)

import numpy as np
CriticalEigenValue = np.loadtxt(open(jobName_Linear + ".csv", "rb"), delimiter=",",
                      skiprows=1)

mdb.models[NLBModel].Moment(cm1=-1000.0*CriticalEigenValue[0,1], createStepName='AppliedPreLoad', 
    distributionType=UNIFORM, field='', localCsys=None, name='Load-1', region=
    mdb.models[NLBModel].rootAssembly.sets['Ref-2'])

mdb.models[NLBModel].Moment(cm1=-1000.0, createStepName='NonLinearBucklingStep', 
    distributionType=UNIFORM, field='', localCsys=None, name='Load-2', region=
    mdb.models[NLBModel].rootAssembly.sets['Ref-2'])

jobName_NonLinear = 'CTLTs_NonLinearBuckling_X-' + ParName

mdb.Job(model=NLBModel, name=jobName_NonLinear, numCpus=12, numDomains=12, numGPUs=1)
mdb.saveAs(pathName=jobName_NonLinear)
mdb.jobs[jobName_NonLinear].submit()
mdb.jobs[jobName_NonLinear].waitForCompletion()
print('======================================================')
print("NonLinear Buckling Analysis Completed First Running!!!")
print('======================================================')
#


# ##############################################################
# ## Perform NonLinear Buckling Analysis and Extract EigenValues
# ## http://www.anning003.com/extract-nonlinear-eigenvalues/
# ##############################################################

stepName_Applied = 'AppliedPreLoad'

from odbAccess import*
from abaqusConstants import*
import string
import numpy as np
import os
import re
from os import path

odb_NonLinear = openOdb(path = jobName_NonLinear+'.odb')
num_Increment = len(odb_NonLinear.steps[stepName_Applied].frames)
i_num = 0
jobName_Restart = 'CTLTs_NonLinearBuckling_X-Restart-' + ParName
while True:
    if path.exists(jobName_NonLinear + '.lck') == True:
        os.remove(jobName_NonLinear + '.lck')
    if path.exists(jobName_Restart + '.lck') == True:
        os.remove(jobName_Restart + '.lck')
    num_Increment -= 1
    i_num += 1
    mdb.models[NLBModel].setValues(restartIncrement=num_Increment, restartJob=
        jobName_NonLinear, restartStep=stepName_Applied)
    mdb.Job(model=NLBModel, name=jobName_Restart, numCpus=12, numDomains=12, numGPUs=1, type=RESTART)
    mdb.saveAs(pathName=jobName_Restart)
    mdb.jobs[jobName_Restart].submit()
    mdb.jobs[jobName_Restart].waitForCompletion()
    print("NonLinear Buckling Analysis Restarted " + str(i_num) + " Times!!")
    odb_Restart = openOdb(path = jobName_Restart + '.odb')
#
    try:
        if len(odb_Restart.steps['NonLinearBucklingStep'].frames) > 1:
            break
    except:
        print("NonLinear Buckling Analysis Restarted " + str(i_num) + " Failedly.....")
    odb_Restart.close() 

odb_Restart = openOdb(path = jobName_Restart + '.odb')
stepName_NonLinear = 'NonLinearBucklingStep'
outfile = open(jobName_Restart + '.csv', 'w')
outfile.write('Mode Order' + ',' + 'NonLinear Critical Load [N m]'+ '\n')
for fm in range(1, len(odb_Restart.steps[stepName_NonLinear].frames)):
    text = odb_Restart.steps[stepName_NonLinear].frames[fm].description
    numeric_const_pattern = '[-+]? (?: (?: \\d* \\. \\d+ ) | (?: \\d+ \\.? ) )(?: [Ee] [+-]? \\d+ ) ?'
    rx = re.compile(numeric_const_pattern, re.VERBOSE)
    EigenValue = [float(i) for i in rx.findall(text)]
    outfile.write(str(int(EigenValue[0])) + ',' + str(EigenValue[1]) + '\n')
outfile.close()
odb_Restart.close()


outfile = open( jobName_NonLinear + '.csv', 'w')
outfile.write('Diaplacement UR1 [rad]' + ',' + 'Reaction Moment RM1 [N m]' + '\n')

for fm in range(0, len(odb_NonLinear.steps[stepName_Applied].frames)-i_num+1):
    timeFrame = odb_NonLinear.steps[stepName_Applied].frames[fm]
    readNode1 = odb_NonLinear.rootAssembly.nodeSets['REF-1']
    ReForce = timeFrame.fieldOutputs['RM']
    readNodeRF = ReForce.getSubset(region=readNode1)
    readNodeRFValues = readNodeRF.values
    ReactionForce = np.zeros(len(odb_NonLinear.steps[stepName_Applied].frames))
    ReactionForce[fm] = readNodeRFValues[0].data[0]/1000.0
    readNode2 = odb_NonLinear.rootAssembly.nodeSets['REF-2']
    Disp = timeFrame.fieldOutputs['UR']
    Displacement = np.zeros(len(odb_NonLinear.steps[stepName_Applied].frames))
    readNodeDisp = Disp.getSubset(region=readNode2)
    readNodeDispValues = readNodeDisp.values
    Displacement[fm] = -readNodeDispValues[0].data[0] # 0-X Direction; 1-Y Direction; 2-Z Direction
    outfile.write(str(Displacement[fm]) + ',' + str(ReactionForce[fm]) + '\n')

outfile.close()

outfile = open( jobName_NonLinear + '_ALLIncrements.csv', 'w')
outfile.write('Diaplacement UR1 [rad]' + ',' + 'Reaction Moment RM1 [N m]' + '\n')

for fm in range(0, len(odb_NonLinear.steps[stepName_Applied].frames)):
    timeFrame = odb_NonLinear.steps[stepName_Applied].frames[fm]
    readNode1 = odb_NonLinear.rootAssembly.nodeSets['REF-1']
    ReForce = timeFrame.fieldOutputs['RM']
    readNodeRF = ReForce.getSubset(region=readNode1)
    readNodeRFValues = readNodeRF.values
    ReactionForce = np.zeros(len(odb_NonLinear.steps[stepName_Applied].frames))
    ReactionForce[fm] = readNodeRFValues[0].data[0]/1000.0
    readNode2 = odb_NonLinear.rootAssembly.nodeSets['REF-2']
    Disp = timeFrame.fieldOutputs['UR']
    Displacement = np.zeros(len(odb_NonLinear.steps[stepName_Applied].frames))
    readNodeDisp = Disp.getSubset(region=readNode2)
    readNodeDispValues = readNodeDisp.values
    Displacement[fm] = -readNodeDispValues[0].data[0] # 0-X Direction; 1-Y Direction; 2-Z Direction
    outfile.write(str(Displacement[fm]) + ',' + str(ReactionForce[fm]) + '\n')

outfile.close()

odb_NonLinear.close()
print('======================================================')
print("   NonLinear Buckling Analysis Completed Finally!!    ")
print('======================================================')


import numpy as np
FEData = np.loadtxt(open(jobName_NonLinear + ".csv", "rb"), delimiter=",",
                      skiprows=1)
K = np.polyfit(FEData[:,0],FEData[:,1], 1) # K[0] - Stiffness

NLData = np.loadtxt(open(jobName_Restart + ".csv", "rb"), delimiter=",",
                      skiprows=1)
outfile = open( jobName_Restart + '.csv', 'w')
outfile.write('Mode Order' + ',' + 'NonLinear Critical Load [N m]'+ '\n')
for i in range(0, len(NLData)):
    outfile.write(str(int(NLData[i,0])) + ',' + str(str(NLData[i,1] + FEData[FEData.shape[0]-1, 1]) + '\n'))

outfile.close()

Data = np.loadtxt(open(jobName_Restart + ".csv", "rb"), delimiter=",",
                      skiprows=1)
outfile = open('CTLTs_NonLinearBuckling_X-' + ParName + '-FinalResults.csv', 'w')
outfile.write('Stiffness [N m/rad]' + ',' + 'NonlinearCriticalLoad [N m]' + '\n')
outfile.write(str(K[0]) + ',' + str(Data[0,1]) + '\n')
outfile.close()

######################################################################
# First PostBuckling Analysis
######################################################################
print('======================================================')
print("   First NonLinear PostBuckling Analysis Started!!    ")
print('======================================================')

jobName_PostBuckling = 'CTLTs_NonLinearPostBuckling_X-' + ParName

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

NPBModel = 'NonLinearPostBucklingModel'
mdb.Model(name=NPBModel, objectToCopy=mdb.models[NLBModel])


#/////////////////////////////////////////////////////////////////////////
#
# IMPERFECTION  need to be started here!!!!!!!!!!
#
#/////////////////////////////////////////////////////////////////////////  
#=========================================================================
# IMPERFECTION
#=========================================================================
# Load NonLinearBuckling Analysis Results
#---------------------------------------------------------------
ODB = openOdb(path = jobName_Restart + '.odb')
pbpPartNodes=ODB.rootAssembly.nodeSets[' ALL NODES']

#CREATE A MATRIX TO SAVE THE NEW COORDINATES OF ALL NODES
NewCoord=np.zeros((len(mdb.models[NPBModel].rootAssembly.instances['Part-1-1'].nodes), 3))

#SELECT THE WANTED BUCKLING MODE AND IMPERFECTION WEIGHT
ImpFrames  = [1, 2]
ImpWeights = [0.05*t, 0.05*t]
for CImp in range(len(ImpFrames)):
    cframe = ImpFrames[CImp]
    firstFrame = ODB.steps['NonLinearBucklingStep'].frames[cframe]
    displacement = firstFrame.fieldOutputs['U']
    pbpDispField = displacement.getSubset(region=pbpPartNodes)
    pbpDisp = pbpDispField.values
# Imperfection Using NonLinearBuckling Analysis Results
#---------------------------------------------------------------
    ind=0;
    IMP = ImpWeights[CImp]
#CALCULATE THE MODIFIED COORDINATES
    for i in mdb.models[NPBModel].rootAssembly.instances['Part-1-1'].nodes:
        NewCoord[ind][0]=i.coordinates[0]+IMP*pbpDisp[ind].data[0]
        NewCoord[ind][1]=i.coordinates[1]+IMP*pbpDisp[ind].data[1]
        NewCoord[ind][2]=i.coordinates[2]+IMP*pbpDisp[ind].data[2]
        ind=ind+1

#SET THE NEW COORDINATES
    mdb.models[NPBModel].rootAssembly.editNode(
    nodes=mdb.models[NPBModel].rootAssembly.instances['Part-1-1'].nodes,
    coordinates=NewCoord)

mdb.models[NPBModel].rootAssembly.regenerate()
print('======================================================')
print("Imperfection has been introduced successfully!!!")
print('======================================================')

mdb.models[NPBModel].loads['Load-1'].setValues(cm1=-1000.0*CriticalEigenValue[0,1]*1.5, 
    distributionType=UNIFORM, field='')

mdb.saveAs(pathName=jobName_PostBuckling)
mdb.Job(model=NPBModel, name=jobName_PostBuckling, numCpus=12, numDomains=12, numGPUs=1)
mdb.saveAs(pathName=jobName_PostBuckling)
mdb.jobs[jobName_PostBuckling].submit()
mdb.jobs[jobName_PostBuckling].waitForCompletion()


# ##############################################################
# ## Perform NonLinear Buckling Analysis and Extract EigenValues
# ## http://www.anning003.com/extract-nonlinear-eigenvalues/
# ##############################################################

stepName_Applied = 'AppliedPreLoad'

from odbAccess import*
from abaqusConstants import*
import string
import numpy as np
import os
import re
from os import path

odb_PostBuckling = openOdb(path = jobName_PostBuckling + '.odb')
num_Increment = len(odb_PostBuckling.steps[stepName_Applied].frames)
i_num = 0
jobName_Restart_P = 'CTLTs_NonLinearPostBuckling_X-Restart-' + ParName
while True:
    if path.exists(jobName_PostBuckling + '.lck') == True:
        os.remove(jobName_PostBuckling + '.lck')
    if path.exists(jobName_Restart_P + '.lck') == True:
        os.remove(jobName_Restart_P + '.lck')
    num_Increment -= 1
    i_num += 1
    mdb.models[NPBModel].setValues(restartIncrement=num_Increment, restartJob=
        jobName_PostBuckling, restartStep=stepName_Applied)
    mdb.Job(model=NPBModel, name=jobName_Restart_P, numCpus=12, numDomains=12, numGPUs=1, type=RESTART)
    mdb.saveAs(pathName=jobName_Restart_P)
    mdb.jobs[jobName_Restart_P].submit()
    mdb.jobs[jobName_Restart_P].waitForCompletion()
    print("NonLinear Buckling Analysis Restarted " + str(i_num) + " Times!!")
    odb_Restart_P = openOdb(path = jobName_Restart_P + '.odb')
    try:
        if len(odb_Restart_P.steps['NonLinearBucklingStep'].frames) > 1:
            break
    except:
        print("NonLinear Buckling Analysis Restarted " + str(i_num) + " Failedly.....")
    odb_Restart_P.close() 

odb_Restart_P = openOdb(path = jobName_Restart_P + '.odb')
stepName_NonLinear = 'NonLinearBucklingStep'
outfile = open(jobName_Restart_P + '.csv', 'w')
outfile.write('Mode Order' + ',' + 'NonLinear Critical Load [N m]'+ '\n')
for fm in range(1, len(odb_Restart_P.steps[stepName_NonLinear].frames)):
    text = odb_Restart_P.steps[stepName_NonLinear].frames[fm].description
    numeric_const_pattern = '[-+]? (?: (?: \\d* \\. \\d+ ) | (?: \\d+ \\.? ) )(?: [Ee] [+-]? \\d+ ) ?'
    rx = re.compile(numeric_const_pattern, re.VERBOSE)
    EigenValue = [float(i) for i in rx.findall(text)]
    outfile.write(str(int(EigenValue[0])) + ',' + str(EigenValue[1]) + '\n')
outfile.close()
odb_Restart_P.close()

outfile = open( jobName_PostBuckling + '.csv', 'w')
outfile.write('Diaplacement UR1 [rad]' + ',' + 'Reaction Moment RM1 [N m]' + '\n')

for fm in range(0, len(odb_PostBuckling.steps[stepName_Applied].frames)-i_num+1):
    timeFrame = odb_PostBuckling.steps[stepName_Applied].frames[fm]
    readNode1 = odb_PostBuckling.rootAssembly.nodeSets['REF-1']
    ReForce = timeFrame.fieldOutputs['RM']
    readNodeRF = ReForce.getSubset(region=readNode1)
    readNodeRFValues = readNodeRF.values
    ReactionForce = np.zeros(len(odb_PostBuckling.steps[stepName_Applied].frames))
    ReactionForce[fm] = readNodeRFValues[0].data[0]/1000.0
    readNode2 = odb_PostBuckling.rootAssembly.nodeSets['REF-2']
    Disp = timeFrame.fieldOutputs['UR']
    Displacement = np.zeros(len(odb_PostBuckling.steps[stepName_Applied].frames))
    readNodeDisp = Disp.getSubset(region=readNode2)
    readNodeDispValues = readNodeDisp.values
    Displacement[fm] = -readNodeDispValues[0].data[0] # 0-X Direction; 1-Y Direction; 2-Z Direction
    outfile.write(str(Displacement[fm]) + ',' + str(ReactionForce[fm]) + '\n')

outfile.close()

outfile = open( jobName_PostBuckling + '_ALLIncrements.csv', 'w')
outfile.write('Diaplacement UR1 [rad]' + ',' + 'Reaction Moment RM1 [N m]' + '\n')

for fm in range(0, len(odb_PostBuckling.steps[stepName_Applied].frames)):
    timeFrame = odb_PostBuckling.steps[stepName_Applied].frames[fm]
    readNode1 = odb_PostBuckling.rootAssembly.nodeSets['REF-1']
    ReForce = timeFrame.fieldOutputs['RM']
    readNodeRF = ReForce.getSubset(region=readNode1)
    readNodeRFValues = readNodeRF.values
    ReactionForce = np.zeros(len(odb_PostBuckling.steps[stepName_Applied].frames))
    ReactionForce[fm] = readNodeRFValues[0].data[0]/1000.0
    readNode2 = odb_PostBuckling.rootAssembly.nodeSets['REF-2']
    Disp = timeFrame.fieldOutputs['UR']
    Displacement = np.zeros(len(odb_PostBuckling.steps[stepName_Applied].frames))
    readNodeDisp = Disp.getSubset(region=readNode2)
    readNodeDispValues = readNodeDisp.values
    Displacement[fm] = -readNodeDispValues[0].data[0] # 0-X Direction; 1-Y Direction; 2-Z Direction
    outfile.write(str(Displacement[fm]) + ',' + str(ReactionForce[fm]) + '\n')

outfile.close()

odb_PostBuckling.close()
print('===============================================================================================')
print("   NonLinear Buckling Analysis Restarted from Firt PostBuckling Results Completed Finally!!    ")
print('===============================================================================================')

import numpy as np
FEData = np.loadtxt(open(jobName_PostBuckling + ".csv", "rb"), delimiter=",",
                      skiprows=1)
K = np.polyfit(FEData[:,0],FEData[:,1], 1) # K[0] - Stiffness

NLData = np.loadtxt(open(jobName_Restart_P + ".csv", "rb"), delimiter=",",
                      skiprows=1)
outfile = open( jobName_Restart_P + '.csv', 'w')
outfile.write('Mode Order' + ',' + 'NonLinear Critical Load [N m]'+ '\n')
for i in range(0, len(NLData)):
    outfile.write(str(int(NLData[i,0])) + ',' + str(str(NLData[i,1] + FEData[FEData.shape[0]-1, 1]) + '\n'))

outfile.close()

Data = np.loadtxt(open(jobName_Restart_P + ".csv", "rb"), delimiter=",",
                      skiprows=1)
outfile = open('CTLTs_NonLinearPostBuckling_X-' + ParName + '-FinalResults.csv', 'w')
outfile.write('Stiffness [N m/rad]' + ',' + 'NonlinearCriticalLoad [N m]' + '\n')
outfile.write(str(K[0]) + ',' + str(Data[0,1]) + '\n')
outfile.close()

