## Common arithmetics
## @author: Luke Li
from utils import *

# Convenient lambdas
nextp2 = lambda x: 1 << (x - 1).bit_length()  # next power of 2 greater than x
ith_bit = lambda n, i: (n >> i) & 1  # get the i-th bit of nunmber n
ith_word = lambda n, i, w: (n >> (i * w)) & ((1 << w) - 1)  # get the i-th word of number n
create_2d = lambda m, n: [[] * n for _ in range(m)]  # create a m-by-n 2D array
is_power_of_two = lambda n: n > 0 and (n & (n - 1)) == 0  # check if n is power of 2

concatenate = lambda zlist, j, w: zlist[j] & 1 if w == 1 else ((zlist[j] & 1) << (w - 1)) + (
        zlist[j - 1] & ((1 << (w - 1)) - 1))  #
# concatenate the higher and lower value of Zj and Z{j-1}
num_from_list = lambda zlist, w: sum(bit << (index * w) for index, bit in enumerate(zlist))


def decompose(Z, w, m):
    ''' Z = ZM(major part) + ZR(remainder part), where ZR is the least significant (w-m) bits'''
    mask = (1 << w - m) - 1
    ZR = Z & mask
    ZM = Z & ~mask
    return ZM, ZR


def compress(numbers, bits):
    '''
        Ternary tree CSA of a given number list
        NOTE!! bits must be greater than the sum of all numbers, o.w. will cause the loss of precision error
    '''
    # Process the list until we reduce it to two numbers
    while len(numbers) > 2:
        tmp = []
        # Take every three elements and apply CSA
        i = 0
        while i < len(numbers) - 2:
            # Perform CSA on every three elements
            S, C = csa(numbers[i], numbers[i + 1], numbers[i + 2], bits)
            tmp.append(S)
            tmp.append(C)
            i += 3
        # If there are remaining numbers that couldn't form a complete group of three,
        # just pass them to the next level.
        if i < len(numbers):
            tmp.extend(numbers[i:])
        numbers = tmp

    return csa(numbers[0], numbers[1], 0, bits)


def csa(x, y, z, T):
    '''
        3-2 CSA(Carry-save-adders) for three T-bit integers x, y and z
    :param x:
    :param y:
    :param z:
    :param T:
    :return:
'''

    # Initialize sum (S) and carry (C)
    S = 0
    C = 0

    # Iterate through each bit position
    for i in range(T):
        # Extract the ith bit from x, y, z
        xb = (x >> i) & 1
        yb = (y >> i) & 1
        zb = (z >> i) & 1

        # Sum bit (X XOR Y XOR Z)
        Sb = xb ^ yb ^ zb

        # Carry bit ((X AND Y) OR (Y AND Z) OR (Z AND X))
        Cb = (xb & yb) | (yb & zb) | (zb & xb)

        # Set the ith bit in S and C
        S |= (Sb << i)
        C |= (Cb << i)

    # C needs to be shifted left by one (to account for carry bit positions)
    C <<= 1

    return S, C


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
