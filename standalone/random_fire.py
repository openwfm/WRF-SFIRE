#!/usr/bin/env python

from math import sqrt
from random import normalvariate as randn
from random import randint,seed
import subprocess as sp

sigma=1.
nx=256
ny=256
dx=6.
dy=6.
igntime=2.
runtime=30*60

slpsigma=sigma*.1/sqrt(2)
windsigma=sigma*5./sqrt(2)

cenx=nx*dx/2.
ceny=ny*dy/2.
dz=sqrt(dx**2+dy**2)
ignr=dz*2
historys=runtime
timestep=dz*1.0/6.

ignxsigma=sigma*dx*nx/10.
ignysigma=sigma*dy*ny/10.

seed()
slpx=randn(0,slpsigma)
slpy=randn(0,slpsigma)
windx=randn(0,windsigma)
windy=randn(0,windsigma)
fuel=randint(1,13)
ignx=randn(cenx,ignxsigma)
igny=randn(ceny,ignysigma)

args=['--nx',nx,'--ny',ny,'--dx',dx,'--dy',dy,
      '--windx',windx,'--windy',windy,'--slopex',slpx,'--slopey',slpy,
      '--fuelcat',fuel,'--timestep',timestep,'--runtime',runtime,'--history',historys,
      '--ignx1',ignx,'--ignx2',ignx,'--igny1',igny,'--igny2',igny,
      '--ignr',ignr,'--ignt1',igntime,'--ignt2',igntime]

args=[ str(a) for a in args]

print ' '.join(['fire.py'] +args)
p=sp.Popen(['python','fire.py']+args)
p.communicate()

f=open('params.txt','w')
for i in xrange(0,len(args),2):
    f.write('%s=%s\n' % (args[i][2:],args[i+1]))
f.close()
