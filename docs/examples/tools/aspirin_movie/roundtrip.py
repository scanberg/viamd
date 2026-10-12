import os
import sys, math, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from movie_lib import *
import story
txt=open(os.path.join(ROOT, 'docs/examples/aspirin_binding_movie.via')).read()
lines=[l for l in txt.splitlines() if l.startswith('KeyframeV3=')]
print(len(lines),'keys; field counts', {len(l.split('=')[1].split(',')) for l in lines})
tool=[str(len(lines))]
for l in lines:
    v=[float(x) for x in l.split('=')[1].split(',')]
    assert all(math.isfinite(x) for x in v)
    tool.append(' '.join('%.7g'%x if i not in (10,12,13,14,15,16,20) else '%d'%round(x) for i,x in enumerate(v[:21]))+' 0')
tool.append('0 600 60 0 1 0 0 1'); tool.append('1 %d'%NF); tool.append('0 0 0'+' 0 0 0'*(NF-1)); tool.append('0 60 0.25')
out=subprocess.run([TOOL],input='\n'.join(tool),capture_output=True,text=True).stdout
rows=np.array([[float(x) for x in l.split()] for l in out.strip().split('\n')])
ref=run_tool(story.keys, story.ANCH,(0,60,0.25))
print(rows.shape, 'finite', np.isfinite(rows).all())
print('max diff eye', np.abs(rows[:,5:8]-ref[:,5:8]).max(), 'fov', np.abs(rows[:,8]-ref[:,8]).max(), 'frame', np.abs(rows[:,1]-ref[:,1]).max())
sp=np.linalg.norm(np.diff(rows[:,5:8],axis=0),axis=1)/0.25
print('eye speed max %.1f A/s at t=%.2f'%(sp.max(), rows[sp.argmax(),0]))
for t0,t1 in ((0,4),(4,16),(16,25),(25,30.5),(30.5,44),(44,57),(57,60)):
    m=(rows[:-1,0]>=t0)&(rows[:-1,0]<t1); print(t0,t1,'mean eye speed %.1f  max %.1f'%(sp[m].mean(), sp[m].max()))
