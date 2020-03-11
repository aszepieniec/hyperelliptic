class HyperellipticCurveElement:
    def __init__( self, curve, U, V ):
        self.curve = curve
        self.U = U
        self.V = V

    @staticmethod
    def Cantor( curve, U1, V1, U2, V2 ):
        # 1.
        g, a, b = xgcd(U1, U2)   # a*U1 + b*U2 == g
        d, c, h3 = xgcd(g, V1+V2) # c*g + h3*(V1+V2) = d
        h2 = c*b
        h1 = c*a
        # h1 * U1 + h2 * U2 + h3 * (V1+V2) = d = gcd(U1, U2, V1-V2)

        # 2.
        V0 = (U1 * V2 * h1 + U2 * V1 * h2 + (V1*V2 + curve.f) * h3) / d
        R = U1.parent()
        V0 = R(V0)

        # 3.
        U = U1 * U2 / d^2
        U = R(U)
        V = V0 % U

        while U.degree() > curve.genus:
            # 4.
            U_ = (curve.f - V^2) / U
            U_ = R(U_)
            V_ = -V % U_

            # 5.
            U, V = U_.monic(), V_
        # (6.)

        # 7.
        return U, V

    def parent( self ):
        return self.curve

    def __mul__( self, other ):
        U, V = HyperellipticCurveElement.Cantor(self.curve, self.U, self.V, other.U, other.V)
        return HyperellipticCurveElement(self.curve, U, V)

    def inverse( self ):
        return HyperellipticCurveElement(self.curve, self.U, -self.V)

    def __pow__( self, exp ):
        R = self.U.parent()

        if exp < 0:
            return (self.__pow__(-exp)).inverse()
        if exp == 0:
            return HyperellipticCurveElement(self.curve, R(1), R(0))
        if exp == 1:
            return self

        acc = HyperellipticCurveElement(self.curve, R(1), R(0))

        for b in range(0, ceil(log(exp,2)) + 1):
            B = ceil(log(exp, 2)) - b
            acc = acc * acc
            if exp & ( 1 << B ) != 0:
                acc = acc * self

        return acc
    
    def __eq__( self, other ):
        if self.curve == other.curve and self.V == other.V and self.U == other.U:
            return True
        else:
            return False

def legendre_symbol( a ):
    F = a.parent()
    p = F.order()
    return a^((p - 1) / 2)

def is_quadratic_residue( n ):
    F = n.parent()
    p = F.order()
    if legendre_symbol(n) == F(1) or n == F(0):
        return True
    else:
        return False
 
def tonelli_shanks_sqrt( n ):
    F = n.parent()
    p = F.order()

    assert is_quadratic_residue(n), "not a square (mod p)"

    q = p - 1
    s = 0
    while q % 2 == 0:
        q /= 2
        s += 1

    # p = q * 2^s + 1

    if s == 1:
        return n^((p + 1) // 4)

    for Z in range(2, p):
        z = F(Z)
        if not is_quadratic_residue(z):
            break

    c = z^q
    r = n^((q + 1) / 2)
    t = n^q
    m = s
    t2 = 0
    while (t - 1) != F(0):
        t2 = t * t
        for i in range(1, m):
            if t2 == F(1):
                break
            t2 = t2 * t2
        b = c^( 1 << (m - i - 1) )
        r = r * b
        c = b * b
        t = t * c
        m = i
    return r

class HyperellipticCurve:
    def __init__( self, f ):
        self.R = f.parent()
        self.F = self.R.base_ring()
        self.x = self.R.gen()
        self.f = f
        self.genus = floor((f.degree()-1) / 2)

    @staticmethod
    def Example( ):
        F = FiniteField(2003)
        R.<x> = PolynomialRing(F)
        f = x^5 + 1184*x^3 + 1846*x^2 + 956*x + 560
        return HyperellipticCurve(f)

    @staticmethod
    def Random( field, genus ):
        R.<x> = PolynomialRing(field)
        degree = 2 * genus + 1
        coeffs = [field.random_element() for d in range(degree)]
        f = x^degree + sum(x^i * coeffs[i] for i in range(len(coeffs)))
        while not f.is_irreducible():
            coeffs = [field.random_element() for d in range(degree)]
            f = x^degree + sum(x^i * coeffs[i] for i in range(len(coeffs)))
        return HyperellipticCurve(f)

    def identity( self ):
        return HyperellipticCurveElement(self, self.R(1), self.R(0))

    def random_element( self ):
        roots = []
        while len(roots) != self.genus:
            xi = self.F.random_element()
            yi2 = self.f(xi)
            if not is_quadratic_residue(yi2):
                continue
            roots.append(xi)
            roots = list(set(roots))
        signs = [ZZ(Integers(2).random_element()) for r in roots]

        U = self.R(1)
        for r in roots:
            U = U * (self.x - r)

        V = self.R(0)
        for i in range(len(roots)):
            y = (-1)^(ZZ(Integers(2).random_element())) * tonelli_shanks_sqrt(self.f(roots[i]))
            lagrange = self.R(1)
            for j in range(len(roots)):
                if j == i:
                    continue
                lagrange = lagrange * (self.x - roots[j])/(roots[i] - roots[j])
            V += y * lagrange

        return HyperellipticCurveElement(self, U, V)
        
