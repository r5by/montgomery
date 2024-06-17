## Utils
## @author: Luke Li

import time
import sympy
from functools import reduce
from operator import mul
import random
import json
import signal
from contextlib import contextmanager

# lambdas
random_list = lambda low, high, count: [random.randint(low, high) for _ in range(count)]


def profiler(num_runs=100, enabled=True):
    def decorator(func):
        if not enabled:
            # If profiling is disabled, return the original function unmodified
            return func

        def wrapper(*args, **kwargs):
            total_time = 0
            for _ in range(num_runs):
                start_time = time.time()
                func(*args, **kwargs)
                end_time = time.time()
                total_time += (end_time - start_time)
            average_time = total_time / num_runs
            print(f"Average execution time for {func.__name__}: {average_time} seconds")
            return func(*args, **kwargs)

        return wrapper

    return decorator


def print2D(myArray):
    # mx = len(str(max((sub[0] for sub in myArray))))
    mx = 5  # todo> adjust this accordingly

    for row in myArray:
        print(" ".join(["{:<{mx}}".format(ele, mx=mx) for ele in row]))


def generate_large_primes(count, num_bits=8192):
    """Generate a list of large prime numbers of specified bit length."""
    primes = []
    while len(primes) < count:
        prime = sympy.randprime(2 ** (num_bits - 1), 2 ** num_bits)
        primes.append(prime)
    return primes


def save_to_file(primes, filename="m_primes.data"):
    """Save the list of prime numbers to a file."""
    with open(filename, "w") as file:
        for prime in primes:
            file.write(f"{prime}\n")


def load_primes_from_file(filename="m_primes.data"):
    """Load prime numbers from a file."""
    with open(filename, "r") as file:
        primes = [int(line.strip()) for line in file]
    return primes


def prime_factors(n):
    factors = []
    # Check for number of 2s that divide n
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # n must be odd at this point, thus a skip of 2 (i.e. 3, 5, 7, 9, ..., sqrt(n)) is fine
    for i in range(3, int(n ** 0.5) + 1, 2):
        # While i divides n, add i and divide n
        while n % i == 0:
            factors.append(i)
            n //= i

    # This condition is to check if n is a prime
    # number greater than 2
    if n > 2:
        factors.append(n)

    return factors


def bi_factor(n):
    flist = prime_factors(n)

    if len(flist) == 1:
        return 1, n

    h = len(flist) >> 1
    n1 = reduce(mul, flist[:h])
    return n1, n // n1


def load_config():
    with open('../config.json', 'r') as file:
        return json.load(file)


class TimeoutException(Exception):
    pass


@contextmanager
def timeout(time):
    # Signal handler function
    def raise_timeout(signum, frame):
        raise TimeoutException()

    # Set the signal handler and a timer
    signal.signal(signal.SIGALRM, raise_timeout)
    signal.alarm(time)

    try:
        yield
    finally:
        # Disable the alarm
        signal.alarm(0)
