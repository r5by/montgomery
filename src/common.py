## Common arithmetics
## @author: Luke Li
from utils import *

# Convenient lambdas
nextp2 = lambda x: 1 << (x - 1).bit_length()  # next power of 2 greater than x
ith_bit = lambda n, i: (n >> i) & 1  # get the i-th bit of nunmber n
ith_word = lambda n, i, w: (n >> (i * w)) & ((1 << w) - 1)  # get the i-th word of number n
create_2d = lambda m, n: [[] * n for _ in range(m)]  # create a m-by-n 2D array
is_power_of_two = lambda n: n > 0 and (n & (n - 1)) == 0  # check if n is power of 2


def get_digits_decimal(u, n):
    res = []
    b = 10 ** n
    while u > 0:
        res.append(u % b)
        u //= b

    return res


def get_digit_radix(u, i, b):
    # Convert number u from base 10 to base b
    digits = []
    while u > 0:
        remainder = u % b
        digits.append(remainder)
        u //= b

    # If the requested digit index i is beyond the available digits, return 0 (assuming the higher digits are 0)
    if i >= len(digits):
        return 0

    return digits[i]


@profiler(num_runs=100, enabled=False)
def xgcd(a, b):
    '''
        ax + bx = gcd(a,b), extended gcd without recursion
    :param a:
    :param b:
    :return: x, y and gcd(a,b)
        !!NOTE: x, and y could be negative numbers
    '''
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b != 0:
        q = a // b
        a, b = b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return x0, y0, a


@profiler(num_runs=100, enabled=True)
def bin_gcd(a, b):
    u, v, e = a, b, 1
    while u & 1 == 0 and v & 1 == 0:
        u >>= 1
        v >>= 1
        e <<= 1

    while u != 0:
        while u & 1 == 0:
            u >>= 1
        while v & 1 == 0:
            v >>= 1

        if u >= v:
            u -= v
        else:
            v -= u

    return e * v


@profiler(num_runs=100, enabled=False)
def bin_inv(a, p, b=1):
    '''
        Return b*x % p, for ax = 1 mod p; by default b=1, this shall return a^{-1} % p
    :param a: the given number
    :param b: default to 1
    :param p: the modulo, NOTE p must be Prime!
    :return: (inverse) of a
    '''

    u, v = a, p
    x1, x2 = b, 0

    while u != 1 and v != 1:
        while u & 1 == 0:
            u >>= 1
            x1 = x1 >> 1 if x1 & 1 == 0 else x1 + p >> 1

        while v & 1 == 0:
            v >>= 1
            x2 = x2 >> 1 if x2 & 1 == 0 else x2 + p >> 1

        if u >= v:
            u -= v
            x1 -= x2
        else:
            v -= u
            x2 -= x1

    return x1 % p if u == 1 else x2 % p


def co_prime(a, b):
    _, _, d = xgcd(a, b)
    return d == 1
