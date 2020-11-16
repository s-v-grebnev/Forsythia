#
# (c) 2020 Sergey Grebnev
#
from gfp2 import p, GFp2element
from montgomery import MontyCurve
from montgomery import MontyPoint
from montgomery import isogen2,isogen3,isoex2,isoex3
from parameters import toy_example

e2 = 0
e3 = 0
f = 0
xp2 = GFp2element(0)
xq2 = GFp2element(0)
xr2 = GFp2element(0)
xp3 = GFp2element(0)
xq3 = GFp2element(0)
xr3 = GFp2element(0)

def ParseParameters(params):
    global p, e2, e3, f, xp2, xq2, xr2, xp3, xq3, xr3
    e2 = params['eA']
    e3 = params['eB']
    f = params['f']
    p = (2 ** e2) * (3 ** e3) * f - 1
    xp2 = GFp2element(params['xp2'][0], params['xp2'][1])
    xq2 = GFp2element(params['xq2'][0], params['xq2'][1])
    xr2 = GFp2element(params['xr2'][0], params['xr2'][1])
    xp3 = GFp2element(params['xp3'][0], params['xp3'][1])
    xq3 = GFp2element(params['xq3'][0], params['xq3'][1])
    xr3 = GFp2element(params['xr3'][0], params['xr3'][1])

global p
print(p)

ParseParameters(toy_example)
print(p)

sk2 = 1234567890
sk3 = 9876543210

pkAlice = isogen2(sk2, e2, xp3, xq3, xr3)
pkBob = isogen3(sk3, e3, xp2, xq2, xq3)

print('pkAlice = ', pkAlice)
print('pkBob = ', pkBob)

j1 = isoex2(sk2, e2, pkBob)
j2 = isoex3(sk3, e3, pkAlice)

print('jAlice = ', j1)
print('jBob = ', j2)