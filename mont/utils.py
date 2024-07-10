## Utils
## @author: Luke Li

import time
import sympy
from functools import reduce
from operator import mul
import random
import json
import os
import signal
from contextlib import contextmanager

from math import sqrt
from numpy import array
from mpmath import mpf

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
    default_config = {
        "COMP_TYPE": "ter",
        "DEFAULT_MUL_OPT": "real8",
        "DEFAULT_EXP_OPT": "mont-ladder-safe",
        "TOTAL_CORES": 6
    }
    config_filename = 'config.json'

    # Determine the script directory (assumed to be the project root)
    script_directory = os.path.dirname(os.path.abspath(__file__))
    project_root_path = os.path.join(script_directory, '..', )
    project_root_config_path = os.path.join(project_root_path, config_filename)

    # Paths to check for the config file
    paths_to_check = [
        os.path.join(os.getcwd(), config_filename),  # Current Working Directory
        os.path.normpath(project_root_config_path)  # Project Root Directory
    ]

    # Check each path for the config file
    for path in paths_to_check:
        if os.path.exists(path):
            with open(path, 'r') as file:
                return json.load(file)

    # 2. Check if the environment variable MONTGOMERY_CONFIG_PATH is set and file exists
    env_config_path = os.getenv('MONTGOMERY_CONFIG_PATH')
    if env_config_path:
        if os.path.exists(env_config_path):
            with open(env_config_path, 'r') as file:
                return json.load(file)

    # 3. If no configuration found, create a default config.json in the current working directory
    with open(os.path.join(os.getcwd(), config_filename), 'w') as file:
        json.dump(default_config, file, indent=4)
    return default_config


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


#region Find addition chain
# I steal it shamelessly from:
#   https://rosettacode.org/wiki/Addition-chain_exponentiation#Python
'''  Rosetta Code task Addition-chain_exponentiation  '''


class AdditionChains:
    ''' two methods of calculating addition chains '''

    def __init__(self):
        ''' memoization for knuth_path '''
        self.chains, self.idx, self.pos = [[1]], 0, 0
        self.pat, self.lvl = {1: 0}, [[1]]

    def add_chain(self):
        ''' method 1: add chains depth then breadth first until done '''
        newchain = self.chains[self.idx].copy()
        newchain.append(self.chains[self.idx][-1] +
                        self.chains[self.idx][self.pos])
        self.chains.append(newchain)
        if self.pos == len(self.chains[self.idx]) - 1:
            self.idx += 1
            self.pos = 0
        else:
            self.pos += 1
        return newchain

    def find_chain(self, nexp):
        ''' method 1 interface: search for chain ending with n, adding more as needed '''
        assert nexp > 0
        if nexp == 1:
            return [1]
        chn = next((a for a in self.chains if a[-1] == nexp), None)
        if chn is None:
            while True:
                chn = self.add_chain()
                if chn[-1] == nexp:
                    break

        return chn

    def knuth_path(self, ngoal):
        ''' method 2: knuth method, uses memoization to search for a shorter chain '''
        if ngoal < 1:
            return []
        while not ngoal in self.pat:
            new_lvl = []
            for i in self.lvl[0]:
                for j in self.knuth_path(i):
                    if not i + j in self.pat:
                        self.pat[i + j] = i
                        new_lvl.append(i + j)

            self.lvl[0] = new_lvl

        returnpath = self.knuth_path(self.pat[ngoal])
        returnpath.append(ngoal)
        return returnpath


def cpow(xbase, chain):
    ''' raise xbase by an addition exponentiation chain for what becomes x**chain[-1] '''
    pows, products = 0, {0: 1, 1: xbase}
    for i in chain:
        products[i] = products[pows] * products[i - pows]
        pows = i
    return products[chain[-1]]

#endregion
