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

def extended_euclidean(a, b, s1 = 1, s2 = 0, t1 = 0, t2 = 1):
    if b == 0:
        return a, s1, t1
    else:
        q = a // b
        return extended_euclidean(b, a % b, s2, s1-(q*s2), t2, t1-(q*t2))

def sieve_of_eratosthenes(limit):
    """Generate a list of prime numbers up to the limit."""
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for start in range(2, int(limit**0.5) + 1):
        if sieve[start]:
            for multiple in range(start * start, limit + 1, start):
                sieve[multiple] = False
    return [num for num, is_prime in enumerate(sieve) if is_prime]

def quadratic_sieve(n):

    B = math.log(n)
    B = math.floor((math.e ** math.sqrt(B * math.log(B)))**(math.sqrt(2) / 4))

    # print(B)

    root = math.ceil(math.sqrt(n))
    primes = sieve_of_eratosthenes(B)
    # print(primes)

    X = []
    Q = []
    Q_decomps = []

    ind = root
    while len(X) < len(primes) + 1:
        ind += 1

        q_orig = (ind * ind) % n
        if q_orig == 0:
            continue

        q = q_orig
        factorisation = []

        for prime in primes:
            power = 0

            while q % prime == 0:
                q = q // prime
                power += 1

            factorisation.append(power)
        
        if q == 1:
            Q_decomps.append(factorisation)
            X.append(ind)
            Q.append(q_orig)

    # print(X)
    # print(Q)
    # print(Q_decomps)

    # for row in Q_decomps:
    #     print(row)
    
    # for row in Q_decomps:
    #     print(row)
    # print()

    matrix = np.array(Q_decomps)
    matrix = matrix % 2
    matrix = matrix.T

    # for row in matrix:
    #     print(row)

    # print()
    
    next_row = 0
    for i in range(len(matrix[0])):
        # Pivoting: Find the row with the largest element in the current column
        max_row = -1
        for k in range(next_row, len(matrix)):
            if matrix[k][i] == 1:
                max_row = k

        if max_row == -1:
            continue

        # Swap rows
        matrix[[next_row, max_row]] = matrix[[max_row, next_row]]

        next_row += 1

        if next_row == len(matrix):
            break
        
        # Make all rows below the current row 0 in the current column
        for k in range(next_row, len(matrix)):
            if matrix[k][i] == 1:
                for j in range(i, len(matrix[0])):
                    matrix[k][j] = (matrix[k][j] + 1) % 2
    
    # for row in matrix:
    #     print(row)

    pivot_columns = []
    free_variables = set(range(len(matrix[0])))

    for i in range(len(matrix)-1, -1, -1):
        pivot_found = False
        for j in range(len(matrix[0])):
            if matrix[i, j] == 1:
                pivot_columns.append(j)
                free_variables.discard(j)
                pivot_found = True
                break
        if not pivot_found:
            pivot_columns.append(-1)
    
    final_x = 0
    final_q = 0

    while final_x == (-final_q % n) or final_x == (final_q % n):
        solution = np.zeros(len(matrix[0]), dtype=int)
        # print(solution)
        while np.all(solution == 0):
            for variable in free_variables:
                solution[variable] = np.random.randint(0, 2)
        
        # print()
        # print(free_variables)
        # print(pivot_columns)
        # print(solution)
        # print()
        
        for i in range(len(matrix)):
            pivot = pivot_columns[i]
            if pivot == -1:
                continue
                
            for j in range(pivot + 1, len(matrix[0])):
                
                if matrix[-i-1][j] == 1:
                    solution[pivot] = (solution[pivot] + solution[j]) % 2
        
        # print()
        # print(solution)
        # print()

        # for row in Q_decomps:
        #     print(row)
        # print()

        # for i in range(len(solution)):
        #     if solution[i] == 1:
        #         print(X[i], Q[i])
        #         print(Q_decomps[i])


        final_x = 1
        final_q = 1
        final_q_decomp = [row for row, mask in zip(Q_decomps, solution) if mask == 1]
        final_q_decomp = np.sum(final_q_decomp, axis=0)
        # print(final_q_decomp)

        for i in range(len(solution)):
            if solution[i] == 1:
                final_x = (final_x * X[i]) % n
        
        for i in range(len(primes)):
                final_q = (final_q * (primes[i]**(final_q_decomp[i] // 2))) % n
    
    # print(final_x, final_q)

    gcd = extended_euclidean(final_x - final_q, n)[0]

    if gcd > 1 and gcd < n:
        return gcd

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

value = 56468475

# print(quadratic_sieve(value))
# print(brent_pollard_rho_factorise(value, np.random.randint(1, value), 800))
# gaussian_elimination(matrix, vector)

# print(gaussian_elimination(matrix, vector))

trials = 10000
tot_time = 0

for i in range(trials):
    # print(i)
    value = np.random.randint(100, 100000000)
    start = time.time()
    # quadratic_sieve(value)
    brent_pollard_rho_factorise(value, np.random.randint(1, value), 800)
    tot_time += time.time() - start
print(tot_time / trials)

# print(quadratic_sieve(20, 10))
