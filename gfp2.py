"""
Class GFp2element: represents GF(p^2)
"""

from ecver.gcd import modinv


class GFp2element:
    a = 0
    b = 0
    p = 0

    def __init__(self, a=0, b=0, p=0, bs=16):
        if isinstance(a, int):
            self.a = a
        else:
            self.a = int(a, bs)
        if isinstance(b, int):
            self.b = b
        else:
            self.b = int(b, bs)
        if isinstance(p, int):
            self.p = p
        else:
            self.p = int(p, bs)

    def __str__(self):
        return str(self.a) + ' + ' + str(self.b) + " * i"

    def __add__(self, other):
        assert(self.p == other.p)
        return GFp2element((self.a + other.a) % self.p, (self.b + other.b) % self.p, self.p, 16)

    def __sub__(self, other):
        assert(self.p == other.p)
        return GFp2element((self.a - other.a) % self.p, (self.b - other.b) % self.p, self.p, 16)

    def __mul__(self, other):
        if isinstance(other, int):
            return GFp2element((self.a * other) % self.p, (self.b * other) % self.p, self.p, 16)
        assert(self.p == other.p)
        return GFp2element((self.a * other.a - self.b * other.b) % self.p,
                           (self.a * other.b + self.b * other.a) % self.p, self.p, 16)

    def modinv(self):
        j = modinv((self.a * self.a + self.b * self.b) % self.p, self.p)
        return GFp2element((self.a * j) % self.p, (-self.b * j) % self.p, self.p, 16)

    def __truediv__(self, other):
        assert(self.p == other.p)
        j = modinv(other)
        return self * j

    def __eq__(self, other):
        if isinstance(other, int):
            return (self.b == 0) and (self.a == other)
        assert(self.p == other.p)
        return self.a == other.a and self.b == other.b

    def __ne__(self, other):
        return not (self == other)
