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
A = GFp2element(0)
C = GFp2element(1)
E0 = MontyCurve(A, C)

def ParseParameters(params):
    global p, e2, e3, f, xp2, xq2, xr2, xp3, xq3, xr3, A, C
    e2 = params['eA']
    e3 = params['eB']
    f = params['f']
    p = (2 ** e2) * (3 ** e3) * f - 1
    A = GFp2element(params['A'][0], params['A'][1])
#    C = GFp2element(params['C'][0], params['C'][1])
    E0 = MontyCurve(A, GFp2element(1))
    xp2 = GFp2element(params['xp2'][0], params['xp2'][1]) // GFp2element(params['zp2'][0], params['zp2'][1])
    xq2 = GFp2element(params['xq2'][0], params['xq2'][1]) // GFp2element(params['zq2'][0], params['zq2'][1])
    yp2 = GFp2element(params['yp2'][0], params['yp2'][1]) // GFp2element(params['zp2'][0], params['zp2'][1])
    yq2 = GFp2element(params['yq2'][0], params['yq2'][1]) // GFp2element(params['zq2'][0], params['zq2'][1])
#    xr2 = GFp2element(params['xr2'][0], params['xr2'][1])
    l = (yp2 + yq2) // (xp2 - xq2)
    xr2 = l * l - (xp2 + xq2) - A

    xp3 = GFp2element(params['xp3'][0], params['xp3'][1]) // GFp2element(params['zp3'][0], params['zp3'][1])
    xq3 = GFp2element(params['xq3'][0], params['xq3'][1]) // GFp2element(params['zq3'][0], params['zq3'][1])
    yp3 = GFp2element(params['yp3'][0], params['yp3'][1]) // GFp2element(params['zp3'][0], params['zp3'][1])
    yq3 = GFp2element(params['yq3'][0], params['yq3'][1]) // GFp2element(params['zq3'][0], params['zq3'][1])

#    xr3 = GFp2element(params['xr3'][0], params['xr3'][1])
    l = (yp3 + yq3) // (xp3 - xq3)
    xr3 = l * l - (xp3 + xq3) - A

global p
print(p)

ParseParameters(toy_example)
print(p)

sk2 = 123
sk3 = 321

pkAlice = isogen2(sk2, e2, xp3, xq3, xr3)
print('pkAlice = ', pkAlice)

pkBob = isogen3(sk3, e3, xp2, xq2, xr2)
print('pkBob = ', pkBob)

j1 = isoex2(sk2, e2, pkBob)
j2 = isoex3(sk3, e3, pkAlice)

print('jAlice = ', j1)
print('jBob = ', j2)