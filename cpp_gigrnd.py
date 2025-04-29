"""
Random variate generator for the generalized inverse Gaussian distribution.
Reference: L Devroye. Random variate generation for the generalized inverse Gaussian distribution.
           Statistics and Computing, 24(2):239-246, 2014.
"""

import math
from numpy import random
from time import perf_counter
import gigrnd_imp as gi

def psi(x, alpha, lam):
    f = -alpha*(math.cosh(x)-1.0)-lam*(math.exp(x)-x-1.0)
    return f

def dpsi(x, alpha, lam):
    f = -alpha*math.sinh(x)-lam*(math.exp(x)-1.0)
    return f

def g(x, sd, td, f1, f2):
    if (x >= -sd) and (x <= td):
        f = 1.0
    elif x > td:
        f = f1
    elif x < -sd:
        f = f2

    return f

def gigrnd(p, a, b):
    # setup -- sample from the two-parameter version gig(lam,omega)
    p = float(p)
    a = float(a)
    b = float(b)
    lam = p
    omega = math.sqrt(a*b)

    if lam < 0:
        lam = -lam
        swap = True
    else:
        swap = False

    alpha = math.sqrt(math.pow(omega,2)+math.pow(lam,2))-lam

    # find t
    x = -psi(1.0, alpha, lam)
    if (x >= 0.5) and (x <= 2.0):
        t = 1.0
    elif x > 2.0:
        if (alpha == 0) and (lam == 0):
            t = 1.0
        else:
            t = math.sqrt(2.0/(alpha+lam))
    elif x < 0.5:
        if (alpha == 0) and (lam == 0):
            t = 1.0
        else:
            t = math.log(4.0/(alpha+2.0*lam))

    # find s
    x = -psi(-1.0, alpha, lam)
    if (x >= 0.5) and (x <= 2.0):
        s = 1.0
    elif x > 2.0:
        if (alpha == 0) and (lam == 0):
            s = 1.0
        else:
            s = math.sqrt(4.0/(alpha*math.cosh(1)+lam))
    elif x < 0.5:
        if (alpha == 0) and (lam == 0):
            s = 1.0
        elif alpha == 0:
            s = 1.0/lam
        elif lam == 0:
            s = math.log(1.0+1.0/alpha+math.sqrt(1.0/math.pow(alpha,2)+2.0/alpha))
        else:
            s = min(1.0/lam, math.log(1.0+1.0/alpha+math.sqrt(1.0/math.pow(alpha,2)+2.0/alpha)))

    # find auxiliary parameters
    eta = -psi(t, alpha, lam)
    zeta = -dpsi(t, alpha, lam)
    theta = -psi(-s, alpha, lam)
    xi = dpsi(-s, alpha, lam)

    p = 1.0/xi
    r = 1.0/zeta

    td = t-r*eta
    sd = s-p*theta
    q = td+sd

    # random variate generation
    while True:
        U = random.random()
        V = random.random()
        W = random.random()
        if U < q/(p+q+r):
            rnd = -sd+q*V
        elif U < (q+r)/(p+q+r):
            rnd = td-r*math.log(V)
        else:
            rnd = -sd+p*math.log(V)

        f1 = math.exp(-eta-zeta*(rnd-t))
        f2 = math.exp(-theta+xi*(rnd+s))
        if W*g(rnd, sd, td, f1, f2) <= math.exp(psi(rnd, alpha, lam)):
            break

    # transform back to the three-parameter version gig(p,a,b)
    rnd = math.exp(rnd)*(lam/omega+math.sqrt(1.0+math.pow(lam,2)/math.pow(omega,2)))
    if swap:
        rnd = 1.0/rnd

    rnd = rnd/math.sqrt(a/b)
    return rnd

def test(fn, name):
    start = perf_counter()
    result = fn()
    duration = perf_counter() - start
    print('{} took {:.0f} microseconds\n\n'.format(name, duration * 1e6))
    print(result[:50])

if __name__ == "__main__":
    # Python Implementation
    random.seed(42)
    
    print(f"testing psi(0, 1, 1): {psi(0, 1, 1)}")
    print(f"testing psi(1, 1, 1): {psi(1, 1, 1)}")
    print(f"testing dpsi(0, 1, 1): {dpsi(0, 1, 1)}")
    print(f"testing dpsi(1, 1, 1): {dpsi(1, 1, 1)}")
    ## gigrnd
    print(f"testing g(1.0, 2.0, 0.5, 0.3, 0.4): {g(1.0, 2.0, 0.5, 0.3, 0.4)}")
    print(f"testing gigrnd(0.5, 2.0, 0.5): {gigrnd(0.5, 2.0, 0.5)}")
    test(lambda: [gigrnd(0.5, 2.0, 0.5) for _ in range(5000)], '(Python implementation)')

    random.seed(42)

    # C++ Implementation
    gi.set_seed(42)

    print(f"testing gi.psi(0, 1, 1): {gi.psi(0, 1, 1)}")
    print(f"testing gi.psi(1, 1, 1): {gi.psi(1, 1, 1)}")
    print(f"testing gi.dpsi(0, 1, 1): {gi.dpsi(0, 1, 1)}")
    print(f"testing gi.dpsi(1, 1, 1): {gi.dpsi(1, 1, 1)}")
    print(f"testing gi.g(1.0, 2.0, 0.5, 0.3, 0.4): {gi.g(1.0, 2.0, 0.5, 0.3, 0.4)}")
    ## gigrnd
    print(f"testing gi.gigrnd(0.5, 2.0, 0.5): {gi.gigrnd(0.5, 2.0, 0.5)}")
    test(lambda: [gi.gigrnd(0.5, 2.0, 0.5) for _ in range(5000)], '(C++ implementation)')
