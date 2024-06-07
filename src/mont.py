## Montgomery domain over F_p, we consider p as a prime other than an arbitrary odd number here
# so that Montgomery domain will be a field in this case.
# In Montgomery Field over F_p, identity is mapped to R%p specifically; and mul operation be like: a_mul_b := REDC(
# a*b), where REDC is the montgomery reduction.
# references:
# [1] ![Montgomery inversion](https://link.springer.com/article/10.1007/s13389-017-0161-x)
# [2] ![High-Radix Design of a Scalable Montgomery Modular Multiplier With Low Latency](
# https://ieeexplore.ieee.org/abstract/document/9328560)
# [3] ![Topics in Computational Number Theory Inspired by Peter L. Montgomery](https://www.cambridge.org/core/books/topics-in-computational-number-theory-inspired-by-peter-l-montgomery/4F7A9AE2CE219D490B7D253558CF6F00)
## @author: Luke Li
from typing import Optional, Union
from common import *
import math


class Montgomery:
    def __init__(self):

        self._N: Optional[int] = None  # modulo

        self.__n: Optional[int] = None  # _R = 2 ** __n; (>>__n) substitutes for (/R) operation
        self._R: Optional[int] = None  # R > N and is a power of 2 (the smallest power of 2 greater than N by default)
        self.__RMASK: Optional[int] = None  # (&__RMASK) substitutes for (%R) operation

        self.__m: Optional[int] = None
        self.__r: Optional[int] = None  # r = 2**m, radix.

        self.mul_opt: Optional[str] = None  # multiplication method configuration (can choose various implementations)
        self.mul = None  # the method

        # multiplication algorithm (realization specific)
        self.n_m: Optional[int] = None  # real-4 specific, R=r**n_m, where r=2**m is the radix

    @classmethod
    def factory(cls, mod: int, mul_opt: str = 'real1') -> 'Montgomery':
        if not mod & 1:
            raise ValueError("The modulus must be an odd number (prime actually).")  # todo> primality test here
            # probably?

        inst = cls()
        inst.N = mod

        return inst.config(mul_opt)

    @property
    def N(self):
        return self._N

    @N.setter
    def N(self, value):
        if not value & 1:
            raise ValueError("The modulus must be an odd number.")  # todo> primality test here probably?
        self._N = value

    @property
    def R(self):
        return self._R

    @R.setter
    def R(self, value):
        if value is not None and not (value & (value - 1)) == 0:
            raise ValueError("R must be a power of 2.")
        self._R = value
        self.__RMASK = self._R - 1
        self.__n = self._R.bit_length() - 1

    @property
    def n(self):
        return self.__n

    @n.setter
    def n(self, value):
        self.__n = value
        self._R = 1 << self.__n
        self.__RMASK = self._R - 1

    @property
    def m(self):
        return self.__m

    @m.setter
    def m(self, value):
        self.__m = value
        self.__r = 1 << self.__m

    @property
    def r(self):
        return self.__r

    @r.setter
    def r(self, value):
        self.__r = value
        self.__m = self.__r.bit_length() - 1

    def build(self, **kwargs: Union[int, str]) -> 'Montgomery':

        # Default properties
        self.n = self.N.bit_length()  # R=2**n
        self.m = self.n  # r = 2**m

        # todo>
        for key, value in kwargs.items():
            if self.mul_opt == 'real1' or self.mul_opt == 'real2':
                if 'R' not in kwargs and 'n' not in kwargs:
                    raise ValueError(f'R (or its power) must be specified for realization 1 or 2.')

                if key != 'R' and key != 'n':
                    continue

                if key == "R":
                    self.R = value
                elif key == "n":
                    self.n = value

                # update radix for mont constant
                self.m = self.n
                self.r = self.R

            elif self.mul_opt == 'real3':
                pass

            elif self.mul_opt == 'real4':
                if 'm' not in kwargs and 'r' not in kwargs:
                    raise ValueError(f' Radix r=2**m must be specified for realization 4 ')

                if key != 'm' and key != 'r':
                    continue

                # set-up radix r=2**m
                if key == "r":
                    self.r = value
                elif key == "m":
                    self.m = value

                # set-up R, we'll use the smallest R here as compared with a wider range of acceptable R's values in
                # the original algorithm described in the paper
                self.n_m = math.ceil(self.n / self.m)
                self.R = self.r ** self.n_m

            elif self.mul_opt == 'real5':
                if 'm' not in kwargs and 'r' not in kwargs:
                    raise ValueError(f' Radix r=2**m must be specified for realization 4 ')

                if key != 'm' and key != 'r':
                    continue

                # set-up radix r=2**m
                if key == "r":
                    self.r = value
                elif key == "m":
                    self.m = value

                # set-up R
                self.d = math.ceil((self.n + self.m + 2) / self.m)
                self.R = self.r ** (self.d - 1) # R = r^{d-1}

        self.__pre_calc()  # Final pre-calculation after all settings
        return self

    def __call__(self, v):
        v_ = self.enter_domain(v)
        return MontgomeryNumber(v_, self)

    def __pre_calc(self):
        self.R2 = self.__pre_calc_R2()  # R^2 % N
        self.N_ = self.__pre_calc_N_()  # N' for rr'-NN'=1

    def enter_domain(self, a):
        return self.multiply(a, self.R2)

    def exit_domain(self, a):
        return self.multiply(a, 1)

    def __pre_calc_R2(self):
        # Use (wn + 1) rounds of modulo addition to get R^2 % N, where R = 2^{wn}
        # refer: <<Topics>> p.19
        ci = self.R % self.N  # c0 = R
        for _ in range(self.n):
            ci = (ci + ci) % self.N

        return ci

    def __pre_calc_N_(self):
        return self.mont_const(self.m)

    def mont_const(self, w):
        '''
            a.k.a Montgomery constant N', for rr'-NN'=1, where r=2^w is the radix that is other than 2 or R
            refer: [3] p.18, Alg. 2.3
        :return:
        '''
        y = 1
        for i in range(2, w + 1):
            if y * self.N & ((1 << i) - 1) == 1:
                continue

            y += (1 << i - 1)

        return (1 << w) - y  # r - y

    def config(self, mul_opt):
        if mul_opt not in ['real1', 'real2', 'real3', 'real4', 'real5']:
            raise ValueError("Unsupported configuration mode for multiplication.")

        if self.mul_opt == mul_opt:
            return self

        self.mul_opt = mul_opt

        #todo>
        if mul_opt == 'real1':
            self.mul = self.real1_FPR2tN
        elif mul_opt == 'real2':
            self.mul = self.real2_FPR2t1
        elif mul_opt == 'real3':
            pass
        elif mul_opt == 'real4':
            self.mul = self.real4_FPR2tm
        elif mul_opt == 'real5':
            self.mul = self.real5_FPR2tm_v1

        return self

    def REDC(self, u):
        '''
            Montgomery reduction (REDC) of number u
            ref: https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
        :param u:
        :return: u * R^{-1} mod N
        '''
        if u < 0 or u > self.R * self.N - 1:
            raise ValueError('Given number if out of montgomery reduction range')

        k = u * self.N_ & self.__RMASK
        t = (u + k * self.N) >> self.__n

        return self.correction(t)

    def correction(self, r):
        if r > 2 * self.N:
            raise ValueError(f'Only number is ready for Montgomery correction step is allowed, check input')

        return r if r < self.N else r - self.N

    def __alm_inv(self, a):
        '''
            Kaliski's almost inversion algorithm
            refer: [1] Alg. 4
        :param a:
        :return:
        '''
        u, v = self.N, a
        r, s = 0, 1
        k = 0

        while v > 0:
            if u & 1 == 0:
                u >>= 1
                s <<= 1
            elif v & 1 == 0:
                v >>= 1
                r <<= 1
            elif u > v:
                u = u - v >> 1
                r += s
                s <<= 1
            else:
                v = v - u >> 1
                s += r
                r <<= 1

            k += 1

        if r >= self.N:
            r -= self.N

        return self.N - r, k

    def __cor_phase_mont(self, r, k):
        '''corretion phase of montgomery inverse
        ref:  [1] Alg. 5
        '''
        n, p = self.__n, self.N

        j = 0
        while j < (n << 1) - k:
            r <<= 1
            if r > p:
                r -= p
            j += 1

        return r

    def __cor_phase_fb(self, r, k):
        '''corretion phase of forward/backward montgomery inversion
        ref:  [1] Alg. 5
        '''
        n, p = self.__n, self.N

        j = 0
        while j < k - n:
            if r & 1:
                r += p

            r >>= 1
            j += 1

        return r

    def _mont_inv(self, a):
        '''
            Calc. the montgomery inverse of a modulo p
            ref:
        :param a: given number in Montgomery representation
        :param R: R = 2^{Wt} >= 2^n, where p < 2^n
        :param p: given odd number
        :return: a*R^{-1} % p
        '''
        t, k = self.__alm_inv(a)
        r = self.__cor_phase_mont(t, k)
        return r

    def real1_FPR2tN(self, a: int, b: int) -> int:
        return self.REDC(a * b)

    def real2_FPR2t1(self, a: int, b: int) -> int:
        '''
            Fixed-precision MM multiplication Radix 2
            refer: [2]. Realization 2
        :param a:
        :param b:
        :return:
        '''
        X, Y, M = a, b, self.N

        Z = 0
        for i in range(self.n):
            Zi = Z + X * ith_bit(Y, i)

            q = Zi & 1  # reason for q = Zi*M_%b, here since b=2, M_=1, b.c. ext(x, 2)=xx, 1 always holds.
            Zi = (Zi + q * M) >> 1
            Z = Zi

        return Z

    def real4_FPR2tm(self, a: int, b: int) -> int:
        '''
            Fixed-Precision MM multiplication with Radix 2**m
            refer: [2] Realization 4
        '''
        X, Y, M, M_, n_m, m = a, b, self.N, self.N_, self.n_m, self.m

        Z = 0
        r = 1 << m
        mask = r - 1  # %r eqs. &mask
        for i in range(n_m):
            Yi = Y & mask  # take Y's ith word
            Zi = Z + X * Yi
            q = (Zi & mask) * M_ & mask
            Zi = (Zi + q * M) >> m  # divides r=2^m
            Z = Zi

            Y >>= m  # prepare Y for next word

        return self.correction(Z)

    def real5_FPR2tm_v1(self, a: int, b: int) -> int:
        '''
            Modified FPR2tm Version 1
            refer: [2] Realization 5
        '''
        X, Y, M, M_, m, d = a, b, self.N, self.N_, self.m, self.d

        Z = 0
        r = 1 << m
        mask = r - 1

        for i in range(d):
            Yi = Y & mask
            q = (Z & mask) * M_ & mask
            Zi = (Z + q * M + (X * Yi << m)) >> m
            Z = Zi

            Y >>= m

        return self.correction(Z)

    def multiply(self, a: int, b: int) -> int:
        # This method delegates to the currently configured multiplication method
        if self.mul is None:
            raise ValueError("Multiplication method not configured (call factory constructor first)")
        return self.mul(a, b)


class MontgomeryNumber:
    def __init__(self, value, mont: Montgomery):
        self.mont = mont
        self.value = value if value < mont.N else value % mont.N

    def __mul__(self, other):
        if isinstance(other, MontgomeryNumber):
            t = self.mont.multiply(self.value, other.value)
            return MontgomeryNumber(t, self.mont)
        elif isinstance(other, int):
            # other = self.mont(other)
            # !! Decision choice !!
            # i * mont(a) should treat i as-is the montgomery repr other than converting it automatically!!
            other = MontgomeryNumber(value=other, mont=self.mont)
            return self.__mul__(other)
        return NotImplemented

    def __rmul__(self, other: int):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, MontgomeryNumber):
            inv = self.mont._mont_inv(other.value)
            return MontgomeryNumber((self * inv).value, self.mont)
        elif isinstance(other, int):
            # Design choice: if we are dividing a '1', we should treat it as the identity in Montgomery Field;
            # o.w. treat it an 'as-is' integer already in Montgomery repr.
            other = self.mont(other) if other == 1 else MontgomeryNumber(other, self.mont)
            return self.__truediv__(other)
        return NotImplemented

    def __rtruediv__(self, other: int):
        t = self.mont(other) if other == 1 else MontgomeryNumber(other, self.mont)
        return t.__truediv__(self)

    def __add__(self, other):
        if isinstance(other, MontgomeryNumber):
            t = self.value + other.value

            # eliminate expensive modulo expression
            if t > self.mont.N:
                t -= self.mont.N

            return MontgomeryNumber(t, self.mont)
        elif isinstance(other, int):
            # Convert int to MontgomeryNumber first
            # other = self.mont(other) # design choice: as-is
            other = MontgomeryNumber(value=other, mont=self.mont)
            return self.__add__(other)
        return NotImplemented

    def __radd__(self, other: int):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, MontgomeryNumber):
            t = self.value - other.value
            if t < 0:
                t += self.mont.N
            return MontgomeryNumber(t, self.mont)
        elif isinstance(other, int):
            # Convert int to MontgomeryNumber first
            other = self.mont(other)
            return self.__sub__(other)
        return NotImplemented

    def __rsub__(self, other: int):
        return self.mont(other).__sub__(self)

    def __eq__(self, other):
        if not isinstance(other, MontgomeryNumber):
            return NotImplemented
        return self.value == other.value and self.mont.N == other.mont.N

    def __repr__(self):
        return f"{self.value} (mod {self.mont.N})"

    def __int__(self):
        # Convert Montgomery representation back to integer using R^{-1}, using int(a)
        return self.mont.exit_domain(self.value)
