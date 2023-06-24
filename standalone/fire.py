#!/usr/bin/env python

from optparse import OptionParser
import subprocess as sp

nmls=('namelist.input','namelist.fire')

usage='usage: %prog [options]'
o=OptionParser(usage)

o.add_option('--nx',action='store',type='int',dest='nx')
o.add_option('--ny',action='store',type='int',dest='ny')
o.add_option('--dx',action='store',type='float',dest='dx')
o.add_option('--dy',action='store',type='float',dest='dy')

o.add_option('--windx',action='store',type='float',dest='windx')
o.add_option('--windy',action='store',type='float',dest='windy')
o.add_option('--slopex',action='store',type='float',dest='slopex')
o.add_option('--slopey',action='store',type='float',dest='slopey')

o.add_option('--runtime',action='store',type='int',dest='runseconds')
o.add_option('--history',action='store',type='int',dest='historys')
o.add_option('--timestep',action='store',type='float',dest='timestep')

for s in ('ignros','ignx1','igny1','ignx2','igny2','ignr','ignt1','ignt2'):
    o.add_option('--'+s,action='store',type='float',dest=s)

fcats=[ str(i) for i in xrange(1,14) ]
o.add_option('--fuelcat',action='store',choices=fcats,dest='fuelcat')

o.add_option('--debug',action='store_true',dest='debug')

o.set_defaults(nx=200,ny=200,dx=6,dy=6,windx=2,windy=1,
               slopex=.5,slopey=-.25,runseconds=60*60,historys=60,
               timestep=.25,ignros=.1,ignx1=500,ignx2=500,
               igny1=500,igny2=500,ignr=10,ignt1=2,ignt2=2,
               fuelcat=3,debug=False)

(opts,args)=o.parse_args()
opts=opts.__dict__

time=opts.pop('timestep')
opts['times']=int(time)
opts['timen']=int( 1000*(time-int(time)) )
opts['timed']=1000

debug=opts.pop('debug')
if debug:
    opts['debuglvl']=2
else:
    opts['debuglvl']=0

opts['fuelcat']=int(opts['fuelcat'])

for n in nmls:
    s=open(n+'.template','r').read()
    open(n,'w').write(s % opts)

p=sp.Popen('./init.exe')
p.communicate()

p=sp.Popen('./fire.exe')
p.communicate()
