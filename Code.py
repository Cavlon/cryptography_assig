import numpy as np
import math
import random
import time

# def modular_pow(a, b, p):
#     power = b
#     current_exponentiation = a
#     result = 1

#     # The binary of the exponent shows which powers of 2 should be included in the result
#     binary = bin(a)

#     for i in range(len(binary)-2):
#         # Only include the required powers
#         if power % 2 == 1:
#             result = (result * current_exponentiation) % p
        
#         # Repeatedly square
#         current_exponentiation = (current_exponentiation * current_exponentiation) % p
#         power = power >> 1
    
#     return result

# Identical to the standard recursive euclidean algorithm but includes the extra equations
def extended_euclidean(a, b, s1 = 1, s2 = 0, t1 = 0, t2 = 1):
    if b == 0:
        return a, s1, t1
    else:
        q = a // b
        return extended_euclidean(b, a % b, s2, s1-(q*s2), t2, t1-(q*t2))

# Miller-Rabin Primality test following Wikipedia pseudocode
def primality_test(n, k):

    # These cases are hord-coded or else it breaks
    if n < 4:
        return n == 2 or n == 3
    elif n % 2 == 0:
        return False
    
    d = n-1
    s = 0
    while (d % 2 == 0):
        d = d >> 1
        s += 1
    
    for i in range(k):
        # a = random.randrange(2, n-1)
        a = int(np.random.random() * (n-3)) + 2
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue

        for i in range(s - 1):

            x = (x * x) % n

            if x == 1:
                return False 
            
            if x == n - 1:
                break
        
        if x != n - 1:
            return False
    return True

def rho_function(val, mod):
    return (pow(val, 2, mod) + 1) % mod

# Algorithm follows the paper's pseudocode
def brent_pollard_rho_factorise(n, x0, m):
    x = x0
    y = x0
    ys = 1

    r = 1
    q = 1

    G = 1

    while G == 1:
        x = y

        for i in range(r):
            y = rho_function(y, n)
        
        k = 0
        while k < r and G == 1:
            ys = y
            count = min(m, r-k)
            for i in range(count):
                y = rho_function(y, n)
                q = (q * abs(x-y)) % n
            G = extended_euclidean(q, n)[0]
            k += m

        r *= 2
    
    if G == n:
        ys = rho_function(ys, n)
        G = extended_euclidean(abs(x - ys), n)[0]
        while G == 1:
            ys = rho_function(ys, n)
            G = extended_euclidean(abs(x - ys), n)[0]
    
    if G == n:
        return False
    else:
        return G

def prime_factorise(val, m):
    small_primes = [2, 3, 5, 7, 11, 13, 17, 19]
    factors = []
    n = val

    # Try some small primes first
    for prime in small_primes:
        if n % prime == 0:
            power = 0

            while n % prime == 0:
                n = n // prime
                power += 1
            
            factors.append((prime, power))

    # Repeatedly find prime factors until they are all found
    while n > 1:

        prime_factor = n

        # Repeatedly find factors of factors until a prime factor is found
        while primality_test(prime_factor, 4) == False:
            factor = prime_factor
            prime_factor = brent_pollard_rho_factorise(factor, int(np.random.random() * (factor-1)) + 1, m)

            # Repeat with a different x0 if Pollard Rho fails
            while prime_factor == False:
                prime_factor = brent_pollard_rho_factorise(factor, int(np.random.random() * (factor-1)) + 1, m)

        # Calculate the power factor
        power = 0
        while n % prime_factor == 0:
            n = n // prime_factor
            power += 1

        # Amend the list of power factors
        factors.append((prime_factor, power))
    
    return factors

def shank_discrete_log(p, g, A, cardinality):
    print(A, cardinality)

    start = time.time()

    target = A

    # Cardinality is used since the sub-group has lower cardinality but the modulus is the same
    m = math.sqrt(cardinality)
    m_large = math.ceil(m * 2)
    m_small = math.ceil(m / 2)

    # Set for fast existence check, array for index search
    big_steps_set = set()
    big_steps_set.add(1)
    big_steps = [1]

    val = 1

    # A big step in the sub-group
    step = pow(g, m_large, p)

    # Compute all the big steps
    for i in range(1, m_small):
        val = (val * step) % p
        big_steps.append(val)
        big_steps_set.add(val)

    print(time.time() - start)
    
    # Iterate through the small steps until a collision
    for i in range(m_large):
        if target in big_steps_set:
            print(time.time() - start)

            return ((big_steps.index(target) * m_large) - i) % (p-1)
        target = (target * g) % p

def DiffieHellman(p, g, B):
    a = np.random.randint(3, p-2)
    A = pow(g, a, p)
    K = pow(B, a, p)
    return A, K

def ElGamalEncrypt(p, g, A, x):
    k = np.random.randint(0, p-1)
    y1 = pow(g, k, p)
    y2 = (x * pow(A, k, p)) % p
    return y1, y2

def ElGamalDecrypt(p, a, y1, y2):
    power = pow(y1, a, p)
    inverse = extended_euclidean(power, p)[1] % p
    return (y2 * inverse) % p

def DiscreteLog(p, g, A):

    # Use Silver-Pholig-Hellman to split it into several smaller discrete log instances
    factors = prime_factorise(p-1, 800)

    print(factors)

    # The final result of the discrete log
    res = 0

    # generator_inverse = extended_euclidean(g, p)[1] % p
    generator_inverse = pow(g, p-2, p)
    p_min_1 = p-1

    for factor, order in factors:

        # These change according to the order of the iteration
        power = p_min_1
        factor_power = 1

        # The generators are constant so they are pre-computed
        generator = pow(g, power // factor, p)

        # The value being logged at each iteration
        beta = A

        # The instance result being built as powers are lifted
        factor_res = 0

        # Iteratively raise the order of the factor until the full power factor is considered
        for i in range(order):
            # Update the current power
            power = power // factor

            # Initialise the beta
            instance_beta = pow(beta, power, p)

            # If the size is 2 then the sqrt rounds to 2 which messes up the modulo calculations
            if factor == 2:
                if instance_beta == 1:
                    coefficient = 0
                else:
                    coefficient = 1
            else:
                # Solve the smaller discrete log instance
                coefficient = shank_discrete_log(p, generator, instance_beta, factor)

            coefficient_result = coefficient * factor_power

            # Update variables
            factor_res += coefficient_result
            factor_power *= factor
            beta *= pow(generator_inverse, coefficient_result, p)

        # Combine the results by the Chinese remainder theorem
        inverse = extended_euclidean(power, factor_power)[1] % factor_power
        res = (res + (power * inverse * factor_res)) % p_min_1

    return res

def DiffieHellmanCrack(p, g, A, B):
    a = DiscreteLog(p, g, A)
    return pow(B,a,p)

def ElGamalCrack(p, g, A, y1, y2):
    a = DiscreteLog(p, g, A)
    power = pow(y1,a,p)
    inverse = extended_euclidean(power, p)[1] % p
    return (y2 * inverse) % p

g = 33026900258717130131626338278
a = int(np.random.random() * 100003100019100043100057100069)
p = 100003100019100043100057100069

# # g = 3
# # a = 5
# # p = 7

A = pow(g,a,p)
print(f"A = {A}")
start = time.time()
print(DiscreteLog(p, g, A))
print(time.time() - start)
# print(recursed_discrete_log(p, g, A, p))
# # print(pow(g,res,p))


# # if (pow(g,res,p) != A):
# #     print("wrong")

# ## PROFILING ##

# trials = 1000
# tot = 0

# for i in range(trials):
#     # a = np.random.randint(0, p-1)
#     # A = pow(g,a,p)

#     start = time.time()
#     res = DiscreteLog(p, g, A)
#     # res = shank_discrete_log(p, g, A, p)
#     end = time.time()
#     tot += end - start
#     if res[0] != a:
#         print("wrong")
#         print(res, a)
# print(tot/trials)

# # print(DiscreteLog(p, g, A))