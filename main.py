#
# implementation of the Forsythia key exchange protocol
# (c) 2020 Sergey Grebnev, s.v.grebnev@yandex.ru
#

from ecver.p import p
from gfp2 import GFp2element
from montgomery import MontgomeryCurve
from montgomery import isogen2, isogen3, isoex2, isoex3
from parameters import toy_example, costello, gost_128
import time
from random import randint

e2 = 0
e3 = 0
f = 0
xp2 = None
xq2 = None
xr2 = None
xp3 = None
xq3 = None
xr3 = None
A = None
C = None
E0 = None

def ParseParameters(params):
    global p, e2, e3, f, xp2, xq2, xr2, xp3, xq3, xr3, A, C, E0
    e2 = params['eA']
    e3 = params['eB']
    f = params['f']
    p = (2 ** e2) * (3 ** e3) * f - 1
    A = GFp2element(params['A'][0], params['A'][1])
#    C = GFp2element(params['C'][0], params['C'][1])
    E0 = MontgomeryCurve(A, GFp2element(1))
    xp2 = GFp2element(params['xp2'][0], params['xp2'][1]) // GFp2element(params['zp2'][0], params['zp2'][1])
    xq2 = GFp2element(params['xq2'][0], params['xq2'][1]) // GFp2element(params['zq2'][0], params['zq2'][1])
    yp2 = GFp2element(params['yp2'][0], params['yp2'][1]) // GFp2element(params['zp2'][0], params['zp2'][1])
    yq2 = GFp2element(params['yq2'][0], params['yq2'][1]) // GFp2element(params['zq2'][0], params['zq2'][1])
# Encode basis points as x-coordinates of P2, Q2, R2=P2-Q2
    l = (yp2 + yq2) // (xp2 - xq2)
    xr2 = l * l - (xp2 + xq2) - A

    xp3 = GFp2element(params['xp3'][0], params['xp3'][1]) // GFp2element(params['zp3'][0], params['zp3'][1])
    xq3 = GFp2element(params['xq3'][0], params['xq3'][1]) // GFp2element(params['zq3'][0], params['zq3'][1])
    yp3 = GFp2element(params['yp3'][0], params['yp3'][1]) // GFp2element(params['zp3'][0], params['zp3'][1])
    yq3 = GFp2element(params['yq3'][0], params['yq3'][1]) // GFp2element(params['zq3'][0], params['zq3'][1])
# Encode basis points as x-coordinates of P3, Q3, R3=P3-Q3
    l = (yp3 + yq3) // (xp3 - xq3)
    xr3 = l * l - (xp3 + xq3) - A

ParseParameters(gost_128)
print("p =", p)
print('E0:', E0, '; j(E0) =', E0.jinv())

start = time.time()
sk2 = randint(0, 2 ** e2)
print('skAlice =', sk2)
pkAlice = isogen2(E0, sk2, e2, xp2, xq2, xr2, xp3, xq3, xr3)
print('pkAlice = ', pkAlice)

#c = MontgomeryCurve(1, 1)
#c.seta(pkAlice[0], pkAlice[1], pkAlice[2])
#print('Alice image curve', c)

sk3 = randint(0, 3 ** e3)
print('skBob =', sk3)
pkBob = isogen3(E0, sk3, e3, xp2, xq2, xr2, xp3, xq3, xr3)
print('pkBob = ', pkBob)

#c.seta(pkBob[0], pkBob[1], pkBob[2])
#print('Bob image curve', c)

j1 = isoex2(sk2, e2, pkBob)
j2 = isoex3(sk3, e3, pkAlice)

print('jAlice = ', j1)
print('jBob = ', j2)

end = time.time()
print('Time elapsed:', end - start, 's')
