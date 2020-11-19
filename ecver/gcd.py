"""
GCD, Extended GCD and Modular inversion implementation
(c) 2018 Sergey Grebnev, s.v.grebnev@yandex.ru
"""

def gcd(a,b):
    while not b == 0:
        r = a % b
        a = b
        b = r
    return a

def egcd(a, b):
    u = 1
    d = a
    if b == 0:
        v = 0
        return (d, u, v)
    else:
        v1 = 0
        v3 = b
    while v3 != 0:
        q = d // v3
        t3 = d % v3
        t1 = u - q * v1
        u = v1
        d = v3
        v1 = t1
        v3 = t3
    v = (d - a * u) // b
    return (d, u, v)

def modinv(a, m):
    d, u, v = egcd(a %m, m)
    if d != 1:
        raise ValueError('Modular inverse does not exist for ' + hex(a).lstrip('0x')[:10] + '... mod ' +
                         hex(m).lstrip('0x')[:10] + '...')
    else:
        return u % m