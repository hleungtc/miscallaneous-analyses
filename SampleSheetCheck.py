import sys

import numpy

filein = open(sys.argv[1], 'r')

with filein as samplesheet:
    cycles = samplesheet.readlines()[12:14]
    forward = cycles[0].rstrip(',\n')
    reverse = cycles[1].rstrip(',\n')
    print("The sequencing run is", forward, "by", reverse, "cycles. Is this correct?")
    
from numpy import genfromtxt
SampleSheet = genfromtxt(sys.argv[1], delimiter=',', skip_header = 21, dtype='str')

unique_projects = numpy.unique(SampleSheet[:,8])
unique_identifiers = numpy.unique(SampleSheet[:,0])



print("Checking for duplicate sample identifiers...")
if len(unique_identifiers) == len(SampleSheet[:,0]):
    print("Good! Sample identifiers are all unique")
else:
    print("Duplicate sample identifiers found. Check Sample_ID column!")
    
print("Displaying unique projects in the sequencing run...")
for project in unique_projects:
    if project.startswith('P-001Z') or project.startswith('_'):
        continue
    else:
        print(project)
        
print("Displaying positive controls in the sequencing run...")
for sample in range(len(SampleSheet)):
    if SampleSheet[sample,1].startswith('S001Z'):
        print(SampleSheet[sample,1], "is on plate", SampleSheet[sample,2], ", well", SampleSheet[sample,3], ", and amplicon", SampleSheet[sample, 8].split("_")[1])



print("Displaying negative controls in the sequencing run...")        
for sample in range(len(SampleSheet)):
    if SampleSheet[sample,1].startswith('DNA'):
        print(SampleSheet[sample,1], "is on plate", SampleSheet[sample,2], ", well", SampleSheet[sample,3], ", and amplicon", SampleSheet[sample, 8].split("_")[1])


print("Please verify that this information is correct")    