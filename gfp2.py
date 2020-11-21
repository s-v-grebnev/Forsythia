"""
Class GFp2element: represents GF(p^2)
(c) 2020 Sergey Grebnev, s.v.grebnev@yandex.ru
"""

from ecver.gcd import modinv
from ecver.globals import getp


class GFp2element:
    a = 0
    b = 0
    global p

    def __init__(self, a=0, b=0):
        assert (isinstance(a, int))
        assert (isinstance(b, int))
        self.a = a
        self.b = b

    def __str__(self):
        return str(self.a) + ' + ' + str(self.b) + " * i"

    def __repr__(self):
        return str(self.a) + ' + ' + str(self.b) + " * i"

    def __add__(self, other):
        #        assert(p == other.p)
        if isinstance(other, int):
            return GFp2element(self.a + other, self.b)
        return GFp2element((self.a + other.a) % getp(), (self.b + other.b) % getp())

    def __sub__(self, other):
        #        assert(p == other.p)
        if isinstance(other, int):
            return GFp2element(self.a - other, self.b)
        return GFp2element((self.a - other.a) % getp(), (self.b - other.b) % getp())

    def __mul__(self, other):
        if isinstance(other, int):
            return GFp2element((self.a * other) % getp(), (self.b * other) % getp())
        return GFp2element((self.a * other.a - self.b * other.b) % getp(),
                           (self.a * other.b + self.b * other.a) % getp())

    def modinv(self):
        j = modinv((self.a * self.a + self.b * self.b) % getp(), getp())
        return GFp2element((self.a * j) % getp(), (-self.b * j) % getp())

    def __truediv__(self, other):
        #        assert(p == other.p)
        if isinstance(other, int):
            j = modinv(other, getp())
        else:
            j = other.modinv()
        return self * j

    def __floordiv__(self, other):
        #        assert(p == other.p)
        if isinstance(other, int):
            j = modinv(other, getp())
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

    def iszero(self):
        return self.a == self.b == 0
