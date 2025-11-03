import numpy as np

def primeDecomposition(n):
    """
    renvoie la liste  des couples (p,v_p(n)), p premier divisant n
    """
    l = []
    p = 2
    while p * p <= n and n > 1:
        if n % p == 0:
            alpha = 1
            n = n // p
            while n % p == 0:
                alpha += 1
                n = n // p
            l.append((p, alpha))
        p = 3 if p == 2 else p + 2
    if n > 1:
        l.append((n, 1))
    return l

def sumDivisors(n):
    """
    renvoie la somme des diviseurs de n (https://cp-algorithms.com/algebra/divisors.html)
    """
    s = 1
    for p, a in primeDecomposition(n):
        s *= (p ** (a + 1) - 1) // (p - 1)
    return s

def Eratosthene(n):
    """
    renvoie un tableau numpy t tel que, si p <= n, t[p] = (p est premier)
    """
    t = np.ones(n + 1, dtype = np.bool)
    t[0] = False
    t[1] = False
    m = int(n ** .5) + 1
    p = 2
    while p <= m:
        for k in range(2, n // p + 1):
            t[k * p] = False
        p += 1
        while p <= n and not t[p]:
            p += 1
    return t

def eratosthene(n):
    """
    renvoie la liste ordonnee des nombres premiers <= n
    """
    return np.nonzero(Eratosthene(n))[0].tolist()

def primeOfRank(n):
    """
    renvoie le nieme nombre premier
    """
    m = 2
    primes = [2] 
    while len(primes) < n:
        m *= 2
        primes = eratosthene(m)
    return primes[n - 1]

def gcd(a, b, bezout = False):
    """
    renvoie le pgcd d de a et b si bezout = False
    sinon renvoie (d, u, v) tq d = ua + vb
    """
    if b:
        if bezout:
            q, r = divmod(a, b)
            d, u, v = gcd(b, r, bezout = True)
            return d, v, u - v * q
        else:
            return gcd(b, a % b)
    else:
        return (a, 1, 0) if bezout else a

def phi(n):
    x = n
    for p, _ in primeDecomposition(n):
        x -= x // p
    return x

def isPalindrome(n):
    digits = str(n)
    return digits == digits[::-1]

def pythagoreanTriples(M):
    """
    Un triplet pythagoricien ou TP est (a,b,c) entiers avec 0 < a < b < c et a^2 + b^2 = c^2.  

    pythagoreanTriples(M) engendre les TPs tq a < b <= M

    Méthode
    Un TP (a, b, c) est primitif (TPP) si (a,b) = (1)
    Pour tout TPP (a', b' ,c), le triplet (a, b, c), obtenu en échangeant a' et b' si a' est pair, s'écrit de manière unique sous la forme
    a = m^2 - n^2, b = 2mn, c = m^2 + n^2
    où m > n > 0, (m, n) = (1) et m et n de parités différentes.
    Les TPs sont les (ka, kb, kc) où (a, b, c) est un TPP et k > 0.
    Preuve : https://en.wikipedia.org/wiki/Pythagorean_triple
    """
    for m in range(2, M // 2 + 1):
        nmax = min(M // (2 * m), m - 1)
        for n in range(nmax - (nmax - m + 1) % 2, 0, -2):
            x = m * m - n * n
            if x > M: break
            if gcd(m, n) == 1:
                y = 2 * m * n
                if x > y: x, y = y, x
                z = m * m + n * n
                for k in range(1, M // y + 1):
                    yield k * x, k * y, k * z
