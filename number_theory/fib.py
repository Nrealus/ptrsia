import sys


def fib(n):
    """
    Naive calculation of the n-th Fibonacci number F_n.
    It is not used in the actual problem resolution suggested here.
    """
    if n < 0:
        raise ValueError('Fibonacci number subscript is negative')

    a, b = 0, 1
    if n == 0:
        return a

    for _ in range(n - 1):
        a, b = b, a + b

    return b

def fib_build_sequence_while(predicate) -> [int]:
    """
    Returns a Python list of consequent Fibonacci numbers, from F_0 to F_(k-1), where k is the first integer such as `predicate(k)` is false.
    """
    res = []

    i = 0
    _tmp = 0
    while True:
        if (i == 0):
            _tmp = 0
        elif (i == 1):
            _tmp = 1
        else:
            _tmp = res[i-1]+res[i-2]

        if (predicate(_tmp)):
            res.append(_tmp)
            i+=1
        else:
            break

    return (res,i)

def fib_mod(n, m) -> int:
    """
    Returns F_n mod m, i.e. the remainder of the euclidian division of the n-th Fibonacci number by m.

    The approach is based on recurrence relations giving F_2n and F_(2n+1) using only F_(n) and F(n-1).
    
    We also use the fact that the set Z/mZ is a ring and in particular :
    
    a1≡b1[m] and a2≡b2[m] => a1+a2≡b1+b2[m]

    a1≡b1[m] and a2≡b2[m] => a1*a2≡b1*b2[m]

    a1≡b1[m] => a1^k≡b1^k[m]
    
    Finally, we also use the `fib_build_sequence_while` function.

    Space-complexity : O(m + log(n))
    
    Time-complexity : O(m + log(n))
    
    (See proofs and more details in SOLUTION.md)
    """
    
    # obtain the list of Fibonacci numbers <= m, as well as its.
    (fibs_lower_than_m, length) = fib_build_sequence_while(lambda x : x <= m)
    
    # traces is the list of the remainders of successive divisions by 2 of n.
    # n is halfed until we reach a value lower than `length`, i.e. the number of Fibonacci numbers under m.
    traces = [n%2]
    k = n
    while True:
        if (k < length):
            break
        k = int(k/2)
        traces.append(k%2)

    # initialize r_n = F_k and r_(n-1) = F_(k-1) where k = max{i / F_(i) <= m}
    # r_n is supposed to be the remainder of the division of F_n by m.
    rn = fibs_lower_than_m[k]
    rn_m1 = fibs_lower_than_m[k-1]

    # iterate over traces and compute r_n by using recurrence relations between r_2n, r_(2n+1) (or r_(2n-1)) and r_n, r_(n-1).
    # these are the same recurrence relations than on F_2n, F_(2n+1) (or F_(2n-1)) and F_n, F_(n-1) because of modular congruence properties (Z/mZ being a ring).
    # see more details and proof of the relation in README.md
    l = len(traces)
    for i in range(1,l):
        r2n = (rn + 2*rn_m1)*rn
        r2n_p1 = (rn + rn_m1)**2 + rn**2
        r2n_m1 = rn**2 + rn_m1**2
        
        if traces[l-1-i] == 0:
            rn = r2n % m
            rn_m1 = r2n_m1 % m
        else:
            rn = r2n_p1 % m
            rn_m1 = r2n % m
    
    return rn

#import time
#def main(argv):
#    n = 1000
#    m = 134
#
#    start = time.time_ns()
#    print(fib(n) % m)
#    end = time.time_ns()
#    duration = end - start
#    print(duration)
#    
#    start = time.time_ns()
#    print(fib_mod(n, m))
#    end = time.time_ns()
#    duration = end - start
#    print(duration)
    
def main(argv):
    with open(argv[1], 'r') as input_file:
        for line in input_file:
            if not line.strip():
                continue

            n, m = map(int, line.split())
            print(fib_mod(n, m))

if __name__ == '__main__':
    main(sys.argv)
