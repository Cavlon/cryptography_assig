import numpy as np
import math
import time

def repeated_squaring(generator, target_exponent, modulo):
    power = target_exponent
    current_exponentiation = generator
    result = 1

    # The binary of the exponent shows which powers of 2 should be included in the result
    binary = bin(target_exponent)

    for i in range(len(binary)-2):
        # Only include the required powers
        if power % 2 == 1:
            result = (result * current_exponentiation) % modulo
        
        # Repeatedly square
        current_exponentiation = (current_exponentiation * current_exponentiation) % modulo
        power = power >> 1
    
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
        a = np.random.randint(2, n-1)
        x = repeated_squaring(a, d, n)
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

def rho_function(val):
    return val * val + 1

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
            y = rho_function(y) % n
        
        k = 0
        while k < r and G == 1:
            ys = y
            count = min(m, r-k)
            for i in range(count):
                y = rho_function(y) % n
                q = (q * abs(x-y)) % n
            G = extended_euclidean(q, n)[0]
            k += m

        r *= 2
    
    if G == n:
        ys = rho_function(ys) % n
        G = extended_euclidean(abs(x - ys), n)[0]
        while G == 1:
            ys = rho_function(ys) % n
            G = extended_euclidean(abs(x - ys), n)[0]
    
    if G == n:
        return False
    else:
        return G

def prime_factorise(val, m):
    small_primes = [2, 3, 5, 7, 11, 13, 17]
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
        while primality_test(prime_factor, 2) == False:
            factor = prime_factor
            prime_factor = brent_pollard_rho_factorise(factor, np.random.randint(1, factor), m)

            # Repeat with a different x0 if Pollard Rho fails
            while prime_factor == False:
                prime_factor = brent_pollard_rho_factorise(factor, np.random.randint(1, factor), m)

        # Calculate the power factor
        power = 0
        while n % prime_factor == 0:
            n = n // prime_factor
            power += 1
        
        # Amend the list of power factors
        factors.append((prime_factor, power))
    
    return factors

def shank_discrete_log(p, g, A, cardinality):

    target = A

    # Cardinality is used since the sub-group has lower cardinality but the order is the same
    m = math.ceil(math.sqrt(cardinality))

    # Dictionary holding {value:power of the generator to attain this value}
    big_steps = {1:0}

    val = 1
    # A big step in the sub-group
    step = repeated_squaring(g, m, p)

    # Compute all the big steps
    for i in range(1, m):
        val = (val * step) % p
        big_steps[val] = i*m
    
    # Compute all the small steps
    for i in range(m):
        if target in big_steps:
            return (big_steps[target] - i) % (p-1)
        target = (target * g) % p

def DiffieHellman(p, g, B):
    a = np.random.randint(2, p-2)
    print(a)
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

    # Use Silver-Pholig-Hellman to split it into several smaller discrete log instances
    factors = prime_factorise(p-1, 800)

    # print(factors)

    # The final result of the discrete log
    res = 0

    for factor, order in factors:

        # These change according to the order of the iteration
        power = (p-1)
        factor_power = 1

        # The generators are constant so they are pre-computed
        generator = repeated_squaring(g, power // factor, p)
        generator_inverse = extended_euclidean(g, p)[1] % p

        # print(f"gener = {generator}, inve = {generator_inverse}")

        # The value being logged at each iteration
        beta = A

        # The instance result being built as powers are lifted
        factor_res = 0

        # Iteratively raise the order of the factor until the full power factor is considered
        for i in range(order):
            # Update the current power
            power = power // factor

            # Initialise the beta
            instance_beta = repeated_squaring(beta, power, p)

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
            beta *= repeated_squaring(generator_inverse, coefficient_result, p)

        # Combine the results by the Chinese remainder theorem
        inverse = extended_euclidean(power, factor_power)[1] % factor_power
        res = (res + (power * inverse * factor_res)) % (p-1)

    return res

def DiffieHellmanCrack(p, g, A, B):
    a = DiscreteLog(p, g, A)
    return repeated_squaring(B, a, p)

def ElGamalCrack(p, g, A, y1, y2):
    a = DiscreteLog(p, g, A)
    power = repeated_squaring(y1, a, p)
    inverse = extended_euclidean(power, p)[1] % p
    return (y2 * inverse) % p


def find_generator(p):
    factors = prime_factorise(p-1, 800)
    found = True

    for i in range(1, p):

        for factor, order in factors:
            power = (p-1) // factor
            if repeated_squaring(i, power, p) == 1:
                found = False
                break
        
        if found:
            return i
        else:
            found = True
    return False

# g = 7
# a = 178162
# p = 999959
# A = repeated_squaring(g,a,p)
# print(f"A = {A}")

# print(find_generator(p))

# print(shank_discrete_log(p,g,A,p))
# print(pollard_rho_discrete_log(p, g, A, p//3))
# print(DiscreteLog(p, g, A))

# print(extended_euclidean(6,16))

## PROFILING ##

trials = 10000
tot = 0

for i in range(trials):
    # print(i)
    # a = np.random.randint(1,p)
    # A = repeated_squaring(g,a,p)
    value = np.random.randint(100, 100000000)
    start = time.time()
    # res = shank_discrete_log(p, g, A, p)
    # brent_pollard_rho_discrete_log(p,g,A,p//3)
    # pollard_rho_discrete_log(p,g,A,p//3)
    brent_pollard_rho_factorise(value, np.random.randint(1, value), 800)
    # res = DiscreteLog(p,g,A)
    end = time.time()
    # if a != res:
    #     print("incorrect")
    tot += end - start
print(tot/trials)
# print(pollard_rho_discrete_log(p,g,A,p//3))
# print(pollard_rho_discrete_log(p,49,9,p//3, (p-1)//2))
# print(pollard_rho_discrete_log(7,2,4,p//3, (7-1)//2))
# print(shank_discrete_log(p, 49, 9, (p-1)//2))
# print(DiscreteLog(p, g, A))


# trials = 100
# min_m = 1
# min_time = 100000

# for m in range(500,1001):
#     tot = 0

#     for i in range(trials):
#         start = time.time()
#         brent_prime_factorise(np.random.randint(100, 100000000), m)
#         # prime_factorise(104728)
#         end = time.time()
#         tot += end - start
#     if tot/trials < min_time:
#         min_time = tot/trials
#         min_m = m
# print(min_time)
# print(min_m)