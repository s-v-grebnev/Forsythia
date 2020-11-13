"""
Class Montgomery curves
"""
from gfp2 import GFp2element
from gfp2 import p


class MontyCurve:
    A = None  # GFp2element(0, 0, 0, 16)
#    B = None  # GFp2element(0, 0, 0, 16)
    C = None  # GFp2element(0, 0, 0, 16)
    Ap24 = None
    Am24 = None
    C24 = None
    ap24 = None

    def __init__(self, A, C): #, Ap24 = None, Am24 = None, C24 = None, ap24 = None
        if isinstance(A, GFp2element):
            self.A = A
        else:
            self.A = GFp2element(A, 0, 16)
#        if isinstance(B, GFp2element):
#           self.B = B
#        else:
#            self.B = GFp2element(B, 0, 16)
        if isinstance(C, GFp2element):
            self.C = C
        else:
            self.C = GFp2element(C, 0, 16)
        self.Ap24 = self.A + self.C * 2
        self.Am24 = self.A - self.C * 2
        self.C24 = self.C * 4
        self.ap24 = self.Ap24 / self.C

    def __str__(self):
        return 'y^2 = ' + "x^3 + (" + str(self.A / self.C) + ") * x^2 + x"

    def jinv(self):
        j = self.A * self.A
        t1 = self.C * self.C
        t0 = t1 + t1
        t0 = j - t0
        t0 = t0 - t1
        j = t0 - t1
        t1 = t1 * t1
        j = j * t1
        t0 = t0 + t0
        t0 = t0 + t0
        t1 = t0 * t0
        t0 = t1 * t0
        t0 = t0 + t0
        t0 = t0 + t0
        j = j.modinv()
        j = t0 * j
        return j

    def xdbladd(self, P, Q, D):
        assert(P.parent == Q.parent == D.parent == self)
        t0 = P.X + P.Z
        t1 = P.X - P.Z
        x2p = t0 * t0
        t2 = Q.X - Q.Z
        xpq = Q.X + Q.Z
        t0 = t0 * t2
        z2p = t1 * t1
        t1 = t1 * xpq
        t2 = x2p - z2p
        x2p = x2p * z2p
        xpq = self.ap24 * t2
        zpq = t0 - t1
        z2p = xpq + z2p
        xpq = t0 + t1
        z2p = z2p * t2
        zpq = zpq * zpq
        xpq = xpq * xpq
        zpq = D.X * zpq
        xpq = D.Z * xpq
        return [MontyPoint(x2p, z2p, self), MontyPoint(xpq, zpq, self)]

    def ladder3pt(self, m, xP, xQ, xD):
        p0 = MontyPoint(xP, GFp2element(1, 0, 16), self)
        p1 = MontyPoint(xQ, GFp2element(1, 0, 16), self)
        p2 = MontyPoint(xD, GFp2element(1, 0, 16), self)
        self.ap24 = (self.A + 2) / 4
        while m > 0  :
            if m % 2 == 1:
                [p0, p1] = self.xdbladd(p0, p1, p2)
            else:
                [p0, p2] = self.xdbladd(p0, p2, p1)
            m = m // 2
        return p1

    def geta(self, p, q, d):
        t1 = p.X + q.X
        t0 = p.X * q.X
        A = d.X * t1
        A = A + t0
        t0 = t0 * d.X
        A = A - 1
        t0 = t0 + t0
        t1 = t1 + d.X
        t0 = t0 + t0
        A = A * A
        t0 = t0.modinv()
        A = A * t0
        A = A - t1
        self.A = A
        self.C = GFp2element(1, 0, 16)
        self.Ap24 = self.A + self.C * 2
        self.Am24 = self.A - self.C * 2
        self.C24 = self.C * 4
        self.ap24 = self.Ap24 / self.C

    def iso2_curve(self, P2):
        Ap24 = P2.X * P2.X
        C24 = P2.Z * P2.Z
        Ap24 = C24 - Ap24
        A = Ap24 * 4 - C24 * 2
        return MontyCurve(A, C24)

    def iso2_eval(self, P2, Q, image):
        t0 = P2.X + P2.Z
        t1 = P2.X - P2.Z
        t2 = Q.X + Q.Z
        t3 = Q.X - Q.Z
        t0 = t0 * t3
        t1 = t1 * t2
        t2 = t0 + t1
        t3 = t0 - t1
        XQP = Q.X * t2
        ZQP = Q.Z * t3
        return MontyPoint(XQP, ZQP, image)

    def iso4_curve(self, P4):
        K2 = P4.X - P4.Z
        K3 = P4.X + P4.Z
        K1 = P4.Z * P4.Z
        K1 = K1 + K1
        C24 = K1 * K1
        K1 = K1 + K1
        Ap24 = P4.X * P4.X
        Ap24 = Ap24 + Ap24
        Ap24 = Ap24 * Ap24
        A = Ap24 * 4 - C24 * 2
        curve = MontyCurve(A, C24)
        return [curve, K1, K2, K3]

    def iso4_eval(self, K1, K2, K3, Q, image):
        t0 = Q.X + Q.Z
        t1 = Q.X - Q.Z
        XPQ = t0 * K2
        ZPQ = t1 * K3
        t0 = t0 * t1
        t0 = t0 * K1
        t1 = XPQ + ZPQ
        ZPQ = XPQ-ZPQ
        t1 = t1 * t1
        ZPQ = ZPQ * ZPQ
        XPQ = t0 + t1
        t0 = ZPQ - t1
        XPQ = XPQ * t1
        ZPQ = ZPQ * t0
        return MontyPoint(XPQ, ZPQ, image)

    def iso3_curve(self, P3):
        K1 = P3.X - P3.Z
        t0 = K1 * K1
        K2 = P3.X + P3.Z
        t1 = K2 * K2
        t2 = t0 + t1
        t3 = K1 + K2
        t3 = t3 * t3
        t3 = t3 - t2
        t2 = t1 + t3
        t3 = t3 + t0
        t4 = t3 + t0
        t4 = t4 + t4
        t4 = t1 + t4
        Am24 = t2 * t4
        t4 = t1 + t2
        t4 = t4 + t4
        t4 = t0 + t4
        Ap24 = t3 * t4
        A = Ap24 * 2 + Am24 * 2
        C = Ap24 - Am24
        curve = MontyCurve(A, C)
        return [curve, K1, K2]

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
            self.Z = GFp2element(1, 0, 16)

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
        t0 = self.X - self.Z
        t2 = t0 * t0
        t1 = self.X + self.Z
        t3 = t1 * t1
        t4 = t1 + t0
        t0 = t1 - t0
        t1 = t4 * t4
        t1 = t1 - t3
        t1 = t1 - t2
        t5 = t3 * self.parent.Ap24
        t3 = t5 * t3
        t6 = t2 * self.parent.Am24
        t2 = t2 * t6
        t3 = t2 - t3
        t2 = t5 - t6
        t1 = t2 * t1
        t2 = t3 + t1
        t2 = t2 * t2
        x = t2 * t4
        t1 = t3 - t1
        t1 = t1 * t1
        z = t1 * t0
        return MontyPoint(x, z, self.parent)

    def mul2e(self, e):
        res = MontyPoint(self.X, self.Z, self.parent)
        for i in range(1, e):
            res = res.mul2()
        return res

    def mul3e(self, e):
        res = MontyPoint(self.X, self.Z, self.parent)
        for i in range(1, e):
            res = res.mul3()
        return res

    def ladder3pt(self, m, P, Q, D):
        p0 = MontyPoint()



    def get_A(self, P, Q, D):
        assert (P.parent == Q.parent == D.parent)
        return A

    def diffadd(self, other, diff):
        pass

    def scalar(self, k):
        pass

    def __mul__(self, other):
        return self.scalar(other)
