## Common arithmetics
## @author: Luke Li<zhongwei.li@mavs.uta.edu>
from .utils import *

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

# Configuration of compressor
config = load_config()
COMP_TYPE = config.get('COMP_TYPE', 'seq')


#region todo> booth_compressor
# def booth_multiplier(multiplicand, multiplier):
#     # Get the size of the numbers
#     n = max(multiplicand.bit_length(), multiplier.bit_length()) + 1
#     # Create masks for sign extension
#     max_value = (1 << n) - 1
#
#     # Sign extend multiplicand and multiplier
#     if multiplicand & (1 << (n - 1)):
#         multiplicand = multiplicand | (~max_value)
#     if multiplier & (1 << (n - 1)):
#         multiplier = multiplier | (~max_value)
#
#     # Initialize product and carry
#     product = 0
#     carry = 0
#
#     # Booth's encoding for multiplier, extend with zeros
#     multiplier <<= 2
#
#     print(f"Initial values: Multiplicand = {bin(multiplicand)}, Multiplier = {bin(multiplier)}")
#
#     for i in range(n):
#         # Check the two least significant bits of the multiplier and the carry
#         bits = (multiplier & 3) | (carry << 2)
#
#         if bits == 1 or bits == 4:  # 01 or 10
#             product += multiplicand
#         elif bits == 2 or bits == 5:  # 10
#             product -= multiplicand
#
#         # Arithmetic shift right (multiplier and product)
#         carry = multiplier & (1 << 1)
#         multiplier >>= 1
#         multiplier |= (carry << (2 * n - 1))  # Extend sign bit if the second last bit was 1
#         product >>= 1
#         product |= (carry << (2 * n - 1))  # Keep sign extension
#
#         print(f"Step {i + 1}: Product = {bin(product)}, Multiplier = {bin(multiplier)}")
#
#     # Adjust product for final result
#     product >>= 1  # Remove the appended bit
#
#     return product
#
#
# # Example usage:
# multiplicand = 13  # 1101 in binary
# multiplier = -3  # Represented in two's complement
#
# # Simulate Booth's multiplication algorithm
# result = booth_multiplier(multiplicand, multiplier)
# print(f"Final product: {result}, in binary: {bin(result)}")
#endregion

def add_chain_breakdown(n, limit=3):
    ''' Get position of '1's in binary repr. of given integer n to form the addition chain
        Short circuit: we care no more than 3 for addition chain reduction, refer to mont [
    '''
    if n & 1 == 0:
        raise ValueError(f'We care only odd numbers for ')

    cnt = 0
    while n != 0:
        n = n & (n - 1)
        cnt += 1

        if cnt > limit:
            return 0

    return cnt


def decompose(Z, w, m):
    ''' Z = ZM(major part) + ZR(remainder part), where ZR is the least significant (w-m) bits'''
    mask = (1 << w - m) - 1
    ZR = Z & mask
    ZM = Z & ~mask
    return ZM, ZR


def compress_ternary(numbers, _bits=None):
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

    return numbers[0], numbers[1]


def compress_sequential(numbers, _bits=None):
    if len(numbers) < 2:
        raise ValueError(f'Must compress more than 2 numbers')

    if len(numbers) == 2:
        return numbers[0], numbers[1]

    s, c = numbers[0], numbers[1]
    for e in numbers[2:]:
        s, c = csa(s, c, e, _bits)

    return s, c


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


def lcm(a: int, b: int) -> int:
    """Compute the Least Common Multiple of a and b with optimizations."""
    pa, pb = 0, 0

    # Reduce a and b by removing factors of 2
    while a & 1 == 0:
        a >>= 1
        pa += 1
    while b & 1 == 0:
        b >>= 1
        pb += 1

    # Calculate the maximum power of two encountered
    p = max(pa, pb)

    l = _lcm(a, b)

    # Adjust the final LCM for the maximal power of two
    l <<= p

    return l


def _lcm(a: int, b: int) -> int:
    ''' LCM of two odd number a and b'''
    _, _, d = xgcd(a, b)

    if d == 1:
        return a * b

    n = b.bit_length()
    mask = (1 << n) - 1
    # _d, _, _ = xgcd(d, 1 << n)
    _d = -mont_const(d, n)
    return a * ((b * _d) & mask)


def mont_const(x, w):
    '''
        a.k.a Montgomery constant N', for rr'-NN'=1, where r=2^w is the radix that is other than 2 or R
        refer: [3] p.18, Alg. 2.3
    :return:
    '''
    y = 1
    for i in range(2, w + 1):
        if y * x & ((1 << i) - 1) == 1:
            continue

        y += (1 << i - 1)

    return (1 << w) - y  # r - y


def co_prime(a, b):
    _, _, d = xgcd(a, b)
    return d == 1


def get_compressor():
    if COMP_TYPE == 'seq':
        return compress_sequential
    elif COMP_TYPE == 'ter':
        return compress_ternary
    else:
        raise ValueError("Unknown compressor type specified")


# Expose the chosen compress function
compress = get_compressor()
