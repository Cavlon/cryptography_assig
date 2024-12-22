import numpy as np
import math

def repeated_squaring(generator, target_exponent, modulo):
    current_exponentiation = generator
    result = 1

    # The binary of the exponent shows which powers of 2 should be included in the result
    binary = bin(target_exponent)[2:][::-1]

    for i in range(len(binary)):
        # Only include the required powers
        if binary[i] == '1':
            result = (result * current_exponentiation) % modulo
        
        # Repeatedly square
        current_exponentiation = (current_exponentiation * current_exponentiation) % modulo
    
    return result

# Identical to the standard recursive euclidean algorithm but includes the extra equations
def extended_euclidean(a, b, s1 = 1, s2 = 0, t1 = 0, t2 = 1):
    if b == 0:
        return a, s1, t1
    else:
        q = a // b
        return extended_euclidean(b, a % b, s2, s1-(q*s2), t2, t1-(q*t2))

# Miller-Rabin Primality test from Wikipedia
def primality_test(n, k):
    # print(n)

    if n < 4:
        return n == 2 or n == 3
    elif n % 2 == 0:
        return False
    
    d = n-1
    s = 0
    while (d % 2 == 0):
        d = d >> 1
        s += 1
    
    # print(d, s)
    
    for i in range(k):
        a = np.random.randint(2, n-1)
        x = repeated_squaring(a, d, n)
        y = 1
        for j in range(s):
            y = (x * x) % n
            if y == 1 and x != 1 and x != n-1:
                return False
            x = y
        if y != 1:
            return False
    return True

def rho_function(val):
    return val * val + 1

# Possibly use Brent's algorithm instead
def pollard_rho_factorise(n, x0):

    if n % 2 == 0:
        return 2

    x = x0
    y = rho_function(x0) % n
    p = extended_euclidean(x-y, n)[0]
    # print(x, y, p)
    while p == 1:
        x = rho_function(x) % n
        y = rho_function(rho_function(y)) % n
        p = extended_euclidean(x-y, n)[0]
        # print(x, y, p)
    if p < n:
        return p
    else:
        return False

def prime_factorise(val):
    factors = []
    n = val
    while n > 1:

        prime_factor = n

        # Repeatedly find factors of factors until a prime factor is found
        while primality_test(prime_factor, 4) == False:
            factor = prime_factor
            prime_factor = pollard_rho_factorise(factor, np.random.randint(1, factor))

            # Repeat with a different x0 if Pollard Rho fails
            while prime_factor == False:
                prime_factor = pollard_rho_factorise(factor, np.random.randint(1, factor))
        
        # print(prime_factor)

        # Calculate the power factor
        prev_power_factor = 1
        power_factor = prime_factor
        while n % power_factor == 0:
            prev_power_factor = power_factor
            power_factor *= prime_factor
        
        # Amend the list of power factors
        factors.append(prev_power_factor)
        n = n // prev_power_factor

        # print(factors)
        # print(n)
    
    return factors

def shank_discrete_log(p, g, A, cardinality):

    # for i in range(cardinality):
    #     print(repeated_squaring(g,i,p))
    
    # If the size is 2 then the sqrt rounds to 2 which messes up the modulo calculations
    if cardinality == 2:
        if A == 1:
            return 0
        if A == g:
            return 1

    m = math.ceil(math.sqrt(cardinality))
    target = A
    # print(m)
    big_steps = {1:0}

    val = 1
    step = repeated_squaring(g, m, p)

    # print(step)

    for i in range(1, m):
        val = (val * step) % p
        big_steps[val] = i*m
    
    # print(big_steps)
    
    for i in range(m):
        # print(A)
        if target in big_steps:
            return (big_steps[target] - i) % (p-1)
        target = (target * g) % p



def DiffieHellman(p, g, B):
    a = np.random.randint(2, p-2)
    A = repeated_squaring(g, a, p)
    K = repeated_squaring(B, a, p)
    return A, K

def ElGamalEncrypt(p, g, A, x):
    k = np.random.randint(1, p-1)
    y1 = repeated_squaring(g, k, p)
    y2 = (x * repeated_squaring(A, k, p)) % p
    return y1, y2

def ElGamalDecrypt(p, a, y1, y2):
    power = repeated_squaring(y1, a, p)
    inverse = extended_euclidean(power, p)[1] % p
    return (y2 * inverse) % p

def DiscreteLog(p, g, A):

    # Use Pholig Hellman to split it into several smaller discrete log instances
    factors = prime_factorise(p-1)

    print(factors)

    res = 0

    for factor in factors:
        power = (p-1) // factor

        inverse = extended_euclidean(power, factor)[1] % factor

        # Use Shank's algorithm to solve the smaller instances
        instance_result = shank_discrete_log(p, repeated_squaring(g, power, p), repeated_squaring(A, power, p), factor)

        # Combine the results by the Chinese remainder theorem
        res = (res + (power * inverse * instance_result)) % (p-1)
    

    return res

def DiffieHellmanCrack(p, g, A, B):
    pass

def ElGamalCrack(p, g, A, y1, y2):
    pass


# print(DiffieHellman(9, 13, 17))
# print(extended_euclidean(56, 15))

# a = 17
# A = repeated_squaring(2, 17, 17)

# e = ElGamalEncrypt(17, 2, A, 7)
# print(e)
# # print(repeated_squaring(e[0], a, 17))
# # print(extended_euclidean(16, 17)[1] % 17)
# print(ElGamalDecrypt(17, 17, e[0], e[1]))

# print(pollard_rho(52, np.random.randint(1, 52)))
# print(primality_test(39, 4))
# print(primality_test(35, 4))
# print(prime_factorise(756329432))

g = 7
a = 8
p = 23
A = repeated_squaring(g,a,p)

# print(shank_discrete_log(p,g,A))
print(DiscreteLog(p, g, A))