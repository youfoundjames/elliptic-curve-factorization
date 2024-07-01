import numpy as np

# BORING FUNCTIONS USED IN LENSTRA:

def EEA(a, b): # EXTENDED EUCLIDEAN ALGORITHM
    if b == 0:
        return 1, 0, a

    r1, r2, gcd = EEA(b, a % b)
    s1 = r2
    s2 = r1 - (a // b) * r2

    return s1, s2, gcd

def fastExp(base, power, modulus): # FAST (MODULAR) EXPONENTIATION
    result = 1
    base %= modulus

    while power > 0:
        # If the power is odd, multiply the base to the result
        if power % 2 == 1:
            result = (result * base) % modulus

        # Square the base and halve the power
        base = (base * base) % modulus
        power //= 2

    return result

def getInverse(a,m): # FINDING MULTIPLICATIVE INVERSES MOD m
    (r,s,gcd) = EEA(a,m)
    if gcd != 1:
        raise Exception('a is not a unit!')
    return r

# ELLIPTIC CURVE FUNCTIONS:

def EC_addition(E,P,Q):
# E = [A,B,p] corresponds to y^2 = x^3+Ax+B modulo p
# P,Q are points which are inputted as [x,y]
# Output is P+Q, or evidence that E is singular

    [A, B, p] = [E[0] % E[2], E[1] % E[2], E[2]];
    
    if None in P or Q == 0: # if we already have a witness or infinity
        return P
    if None in Q or P == 0:
        return Q
    if (fastExp(4*A,3,p) + fastExp(27*B,2,p)) % p == 0:
        print("The elliptic curve entered is singular.")
        return None
    if P!=0:
        [xP, yP]  = [P[0]%p, P[1]%p];
        if fastExp(yP,2,p) != (fastExp(xP,3,p) + A*xP + B) % p:
            print(str(P)+' is not on the elliptic curve E')
            return None
    if Q!=0:
        [xQ, yQ]  = [Q[0]%p, Q[1]%p];
        if fastExp(yQ,2,p) != (fastExp(xQ,3,p) + A*xQ + B) % p:
            print(str(Q)+' is not on the elliptic curve E')
            return None
    if P != Q:
        if xP == xQ:
            return 0
    if P == Q:
        if yP == 0: 
            return 0

    return EC_addition_nonsymmetric(E,P,Q)

def EC_addition_nonsymmetric(E,P,Q):
    
    [A, B, p] = [E[0] % E[2], E[1] % E[2], E[2]];
    [Px, Py]  = [P[0] % p, P[1] % p];
    [Qx, Qy]  = [Q[0] % p, Q[1] % p];

    if P != Q: 
        num = (Qy - Py) % p
        den = (Qx - Px) % p

    if P == Q:
        num = (3*fastExp(Px,2,p) + A) % p
        den = (2 * Py) % p
    
    witness = EEA(den,p)[2]
    if witness != 1:
        return [None, witness]
    
    m = (num*getInverse(den,p)) % p;
    b = (Py - m * Px) % p;
    
    Rx = m**2 - Px - Qx % p # equating quadratic coefficients
    Ry = (-1) * (m * Rx + b) % p
    return [Rx, Ry]


def EC_fast_multiplication(E,P,n):
    
    expansion = n.digits(2) #returns a list containing the base 2 expansion of n
    k = len(expansion)
    
    def point_doubler(P,k): #takes in P and returns [P, 2P, 4P, ..., (2^k)P]
        accumulator = P
        terms = [P]
        while k>0:
            accumulator = EC_addition(E,accumulator,accumulator)
            terms.append(accumulator)
            k = k-1
        return terms
    
    terms = point_doubler(P,k)
    actualTermsToSum = [a * b for a, b in zip(expansion, terms)] #multiply the ith element of expansion by the ith element of terms
    actualTermsToSum = [term for term in actualTermsToSum if term] #remove empty arrays which appear when that point is "zeroed" out

    accumulator = actualTermsToSum[0]
    for term in actualTermsToSum[1:]:
        accumulator = EC_addition(E,accumulator,term)
        
    return accumulator

def EC_random(p):
# Generates a random EC, and a random point P on this EC.
    [x, y] = [randint(0, p - 1), randint(0, p - 1)];
    A = randint(0, p - 1);
    B = (fastExp(y,2,p) - fastExp(x,3,p) - A * x) % p;
    return([[A,B,p],[x,y]])

def lenstraFactorial(E,P,b):
# Implementation of Lenstra's factoring algorithm;  b specifies how hard to try
# In other words, attempt to compute P, (2!)P, (3!)P, ... (b!)P in an attempt to factor p
# Already assumes that P lies on E
    i=1;
    Q=P;
    while Q != 0 and Q[0] != None and i<b:
        i = i + 1 
        Q = EC_fast_multiplication(E,Q,i)
        #print(Q)
    if i == b or Q == 0:    
        #print('No luck.  Try again, maybe increase b?')
        return None
    return [E,P,i,Q[1]] 

def lenstraDoubling(E,P,b):
    # Already assumes that P is on E
    i=1;
    Q=P;
    while Q and Q != 0 and Q[0] != None and i<b:
        i = i + 1 
        Q = EC_addition(E,Q,Q)
        #print(Q)
    if i == b or Q == 0 or not Q:    
        #print('No luck.  Try again, maybe increase b?')
        return None
    return [E,P,i,Q[1]] 

def random_lenstra(n,b):
    (E,P) = EC_random(n)
    #  CHANGE THE FUNCTION CALL HERE FOR
    #     DIFFERENT IMPLEMENTATIONS:
    return lenstraDoubling(E,P,b)

def lenstraExperimental(n,b,r):
    Bcount = 0
    Cvals = []
    factorProducedWasNItself = 0
    for i in range(0,r):
        lens = random_lenstra(n,b)
        if lens:
            if lens[3] == n:
                factorProducedWasNItself += 1
            Cvals.append(lens[2])
        else:
            Bcount += 1
    return [np.mean(Cvals),Bcount,factorProducedWasNItself]
    
Ns = [2006099159,2319263231,5298998713,7158148691]
n = Ns[1]
b = 1000
r = 100
lenstraExperimental(35,b,r)
