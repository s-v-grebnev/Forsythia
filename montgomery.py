"""
Class Montgomery curves
"""
from gfp2 import GFp2element


class MontyCurve:
    A = None  # GFp2element(0, 0, 0, 16)
    B = None  # GFp2element(0, 0, 0, 16)
    C = None  # GFp2element(0, 0, 0, 16)
    p = 0
    Ap24 = None
    Am24 = None
    C24 = None
    ap24 = None

    def __init__(self, A, B, C, p):
        self.p = p
        if isinstance(A, GFp2element):
            self.A = A
        else:
            self.A = GFp2element(A, 0, self.p, 16)
        if isinstance(B, GFp2element):
            self.B = B
        else:
            self.B = GFp2element(B, 0, self.p, 16)
        if isinstance(C, GFp2element):
            self.C = C
        else:
            self.C = GFp2element(C, 0, self.p, 16)
        self.Ap24 = self.A + self.C * 2
        self.Am24 = self.A - self.C * 2
        self.C24 = self.C * 4
        self.ap24 = self.Ap24 / self.C

    def __str__(self):
        return str(self.B / self.C) + ' * y^2 = ' + "x^3 + (" + str(self.A / self.C) + ") * x^2 + x"


class MontyPoint:
    X = GFp2element(0, 0, 0)
#    Y = GFp2element(0, 0, 0)
    Z = GFp2element(0, 0, 0)
    parent = None

    def __init__(self, X, Z, parent):
        self.parent = parent
        self.X = X
#        self.Y = Y
        if not Z is None:
            self.Z = Z
        else:
            self.Z = GFp2element(1, 0, parent.p, 16)

    def __str__(self):
        return ('(' + str(self.X) + ':' + str(self.Z) + ')')

    def __add__(self, other):
        pass


    def mul2(self):
        t0 = self.X - self.Z
        t1 = self.X + self.Z
        t0 = t0 * t0
        t1 = t1 * t1
        z = t0 * self.parent.C24
        x = z * t1
        t1 = t1 - t0
        t0 = t1 * self.parent.Ap24
        z = z + t0
        z = z * t1
        return MontyPoint(x, z, self.parent)

    def mul3(self):
        pass

    def mul2e(self, e):
        pass

    def mul3e(self, e):
        pass

    def diffadd(self, other, diff):
        pass

    def scalar(self, k):
        pass

    def __mul__(self, other):
        return self.scalar(other)
