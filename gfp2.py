"""
Class GFp2element: represents GF(p^2)
"""

from ecver.gcd import modinv

p = 0

class GFp2element:
    a = 0
    b = 0

    def __init__(self, a=0, b=0,  bs=16):
        if isinstance(a, int):
            self.a = a
        else:
            self.a = int(a, bs)
        if isinstance(b, int):
            self.b = b
        else:
            self.b = int(b, bs)

    def __str__(self):
        return str(self.a) + ' + ' + str(self.b) + " * i"

    def __add__(self, other):
#        assert(p == other.p)
        return GFp2element((self.a + other.a) % p, (self.b + other.b) % p,  16)

    def __sub__(self, other):
#        assert(p == other.p)
        return GFp2element((self.a - other.a) % p, (self.b - other.b) % p,  16)

    def __mul__(self, other):
        if isinstance(other, int):
            return GFp2element((self.a * other) % p, (self.b * other) % p,  16)
        assert(p == other.p)
        return GFp2element((self.a * other.a - self.b * other.b) % p,
                           (self.a * other.b + self.b * other.a) % p,  16)

    def modinv(self):
        j = modinv((self.a * self.a + self.b * self.b) % p, p)
        return GFp2element((self.a * j) % p, (-self.b * j) % p,  16)

    def __truediv__(self, other):
#        assert(p == other.p)
        if isinstance(other, int):
            j = modinv(other, p)
        else:
            j = other.modinv()
        return self * j

    def __eq__(self, other):
        if isinstance(other, int):
            return (self.b == 0) and (self.a == other)
#        assert(p == other.p)
        return self.a == other.a and self.b == other.b

    def __ne__(self, other):
        return not (self == other)

def __init__():
    global p
    p = 241