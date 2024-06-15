## Common arithmetics
## @author: Luke Li<zhongwei.li@mavs.uta.edu>
from utils import *

# Convenient lambdas
nextp2 = lambda x: 1 << (x - 1).bit_length()  # next power of 2 greater than x
ith_bit = lambda n, i: (n >> i) & 1  # get the i-th bit of nunmber n
ith_word = lambda n, i, w: 0 if i < 0 else (n >> (i * w)) & ((1 << w) - 1)  # get the i-th word of number n,
# if i < 0 returns 0
create_2d = lambda m, n: [[] * n for _ in range(m)]  # create a m-by-n 2D array
is_power_of_two = lambda n: n > 0 and (n & (n - 1)) == 0  # check if n is power of 2
extract_bits = lambda n, i, j: (n >> i) & ((1 << (j - i + 1)) - 1)  # extract i-th to j-th bits of given number n (
# inclusively) for i <= j
update_ith_word = lambda n, i, w, v: (n & ~(((1 << w) - 1) << (i * w))) | (v << (i * w))  # Lambda function to
# update the ith word of w bits in number n with value v

concatenate = lambda a, b, w: a << w | b
num_from_list = lambda zlist, w: sum(bit << (index * w) for index, bit in enumerate(zlist))


# Constants
REG_SIZE = 256  # the register size in bits


def decompose(Z, w, m):
    ''' Z = ZM(major part) + ZR(remainder part), where ZR is the least significant (w-m) bits'''
    mask = (1 << w - m) - 1
    ZR = Z & mask
    ZM = Z & ~mask
    return ZM, ZR


def compress(numbers, _bits=None):
    '''
        Ternary tree CSA of a given number list
        NOTE!! _bits must be greater than the sum of all numbers, o.w. will cause the loss of precision error

        TODO> Find out the correct way of compress (including, compress 2/3/>3 numbers) algorithms that works for real7 and real8
        compress output different results based on input number's orders, for example:

        a, b = compress([1, 2, 3, 4, 5, 6])  # result: a=21, b=0
        a_, b_ = compress([4, 3, 1, 2, 6, 5])  # result: a=5, b=16

        For this reason, the realization 7 & 8 in [2] are not specified clearly. line 9 in real-7 or line 12 in real-8
        explicitly calculate 'c' based on intermediate result that is derived from compress of multiple variables.

        Btw, recursively call:
        a, b = compress([a, b, 0])
        shall produce a final result in form of S(!=0), C(==0).

        Thus, how to define compression alg. on 2 elements is also ambiguous. Therefore raise questions regarding
        line 6 in real-7 or line 8 in real-8 for the compression of qS and qC.
    '''
    # Process the list until we reduce it to two numbers
    while len(numbers) > 2:
        tmp = []
        # Take every three elements and apply CSA
        i = 0
        while i < len(numbers) - 2:
            # Perform CSA on every three elements
            S, C = csa(numbers[i], numbers[i + 1], numbers[i + 2], _bits)
            # S, C = csa(numbers[i], numbers[i + 1], numbers[i + 2])
            tmp.append(S)
            tmp.append(C)
            i += 3
        # If there are remaining numbers that couldn't form a complete group of three,
        # just pass them to the next level.
        if i < len(numbers):
            tmp.extend(numbers[i:])
        numbers = tmp

    # return csa(numbers[0], numbers[1], 0)
    return numbers[0], numbers[1]


def csa(x, y, z, _T=REG_SIZE):
    '''
        3-2 CSA(Carry-save-adders) for three T-bit integers x, y and z;
        Generally speaking, csa of three T-bit integers shall produce two (T+1)-bit integers S(sum) and C(carry)
    :param x:
    :param y:
    :param z:
    :param _T: Max bit length of x,y,z if given, o.w. calc. from x,y,z
    :return:
'''

    # Initialize sum (S) and carry (C)
    S = 0
    C = 0

    T = int(_T) if _T else max(x.bit_length(), y.bit_length(), z.bit_length())
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
