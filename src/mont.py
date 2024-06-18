## Montgomery domain over F_p, we consider p as a prime other than an arbitrary odd number here
# so that Montgomery domain will be a field in this case.
# In Montgomery Field over F_p, identity is mapped to R%p specifically; and mul operation be like: a_mul_b := REDC(
# a*b), where REDC is the montgomery reduction.
# references:
# [1] [Montgomery inversion](https://link.springer.com/article/10.1007/s13389-017-0161-x)
# [2] [High-Radix Design of a Scalable Montgomery Modular Multiplier With Low Latency](
# https://ieeexplore.ieee.org/abstract/document/9328560)
# [3] [Topics in Computational Number Theory Inspired by Peter L. Montgomery](https://www.cambridge.org/core/books/topics-in-computational-number-theory-inspired-by-peter-l-montgomery/4F7A9AE2CE219D490B7D253558CF6F00)
## @author: Luke Li<zhongwei.li@mavs.uta.edu>
from typing import Optional, Union
from common import *
import math
from sim import compressor_2m4to2


class Montgomery:
    def __init__(self):

        self._N: Optional[int] = None  # modulo, stated as "M" in [2]

        self.__n: Optional[int] = None  # 2 ** __n > modulo, stated as "N" in [2]
        self._R: Optional[int] = None  # R > N and is a power of 2 (the smallest power of 2 greater than N by default)
        self.__RMASK: Optional[int] = None  # (&__RMASK) substitutes for (%R) operation

        self.__m: Optional[int] = None
        self.__r: Optional[int] = None  # r = 2**m, radix.

        self._w: Optional[int] = None  # a word contains w bits, used for multiple word scenarios

        self.mul_opt: Optional[str] = None  # multiplication method configuration (can choose various implementations)
        self.mul = None  # the method

        # multiplication algorithm (realization specific)
        self.e: Optional[int] = None  # real-3 specific
        self.n_m: Optional[int] = None  # real-4 specific, R=r**n_m, where r=2**m is the radix
        self.d: Optional[int] = None  # real-5 specific

        self.p: Optional[int] = None  # real-8 specific: PE count
        self.k: Optional[int] = None  # real-8 specific

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
        # As the 'M' stated in [2]
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

    @property
    def n(self):
        # As the 'N' stated in [2]
        return self.__n

    @n.setter
    def n(self, value):
        if value is not None and (1 << value) < self.N:
            raise ValueError("The modulo must be smaller than 2**n")
        self.__n = value

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

    @property
    def w(self):
        return self._w

    @w.setter
    def w(self, value):
        self._w = value

    def parse_args(self, **kwargs: Union[int, str]):

        for key, value in kwargs.items():
            if (self.mul_opt == 'real1'
                    or self.mul_opt == 'real2'):
                if 'R' not in kwargs and 'n' not in kwargs:
                    raise ValueError(f'R (or its power) must be specified for realization 1 or 2.')

                if key != 'R' and key != 'n':
                    continue

                if key == "R":
                    self.R = value
                    self.n = value.bit_length() - 1
                elif key == "n":
                    self.n = value
                    self.R = 1 << value

            elif self.mul_opt == 'real3':
                if 'w' not in kwargs:
                    raise ValueError(f'Word length (in bits) must be specified in realization 3.')

                if key != 'w':
                    continue

                self.w = value

            elif self.mul_opt == 'real4':
                if 'm' not in kwargs and 'r' not in kwargs:
                    raise ValueError(f' Radix r=2**m must be specified for realization 4 ')

                if key != 'm' and key != 'r':
                    continue

                # setup radix r=2**m
                if key == "r":
                    self.r = value
                elif key == "m":
                    self.m = value

            elif (self.mul_opt == 'real5'
                  or self.mul_opt == 'real6'
                  or self.mul_opt == 'real7'):
                if 'm' not in kwargs and 'r' not in kwargs:
                    raise ValueError(f' Radix r=2**m must be specified for realization: {self.mul_opt}')

                # setup radix r=2**m
                if key == "r":
                    self.r = value
                elif key == "m":
                    self.m = value

                if key == "w":
                    self.w = value

            elif self.mul_opt == 'real8':
                if 'm' not in kwargs and 'r' not in kwargs:
                    raise ValueError(f' Radix r=2**m must be specified for realization: {self.mul_opt}')

                if 'w' not in kwargs:
                    raise ValueError(f'Word length (in bits) must be specified in realization 8.')

                # set-up radix r=2**m
                if key == "r":
                    self.r = value
                elif key == "m":
                    self.m = value

                # setup word length
                if key == 'w':
                    self.w = value

                # setup PE count, use ideal value only
                # if key == 'p':
                #     self.p = value

            else:
                raise ValueError(f'Unsupported realization!')

    def build(self, **kwargs: Union[int, str]) -> 'Montgomery':

        # Default properties
        self.n = self.N.bit_length()
        self.R = 1 << self.n
        self.m = self.n  # auto: r = 2**m

        if self.mul_opt != 'real1' and not kwargs:
            raise ValueError(f'For multiplication realizations other than the default, current: {self.mul_opt}, '
                             f'requires more input arguments.')

        self.parse_args(**kwargs)

        if (self.mul_opt == 'real1'
                or self.mul_opt == 'real2'):

            # update radix for mont constant
            self.m = self.n
            self.r = self.R

        elif self.mul_opt == 'real3':

            # update default values
            self.R = 1 << self.n + 1  # update R to 2**{N+1}
            self.e = math.ceil((self.n + 2) / self.w) + 1  # e = ceil((N+2)/w) + 1, where R=2**{N+1}

            # update radix for mont constant
            self.m = 1  # auto: r == 2

        elif self.mul_opt == 'real4':

            # setup R, we'll use the smallest R here as compared with a wider range of acceptable R's values in
            # the original algorithm described in the paper
            self.n_m = math.ceil(self.n / self.m)
            self.R = self.r ** self.n_m

        elif self.mul_opt == 'real5':

            # setup R
            self.d = math.ceil((self.n + self.m + 2) / self.m)
            self.R = self.r ** (self.d - 1)  # R = r^{d-1}

        elif (self.mul_opt == 'real6'
              or self.mul_opt == 'real7'):

            # setup R
            self.d = math.ceil((self.n + self.m + 2) / self.m)
            self.R = self.r ** (self.d - 1)  # R = r^{d-1}

            if self.w and self.w < self.m * 2:
                raise ValueError(f'Condition: w (w={self.w}) >= 2m (m={self.m}) is violated for {self.mul_opt}!')

            # w is used for decomposition, default to 2*m
            if self.w is None:
                self.w = 2 * self.m

        elif self.mul_opt == 'real8':

            if self.w < self.m * 2:
                raise ValueError(f'Condition: w (w={self.w}) >= 2m (m={self.m}) is violated for real-8')

            # if self.w >= self.n:
            #     raise ValueError(f'w(={self.w}) cannot be greater than N(={self.n})')

            # setup R
            self.e = math.ceil((self.n + 2 * self.m + 2) / self.w)
            self.p = math.ceil((self.n + self.m + 2) / self.m)  # ref [2] section 3.2.1, ideal PE count
            self.k = math.ceil((self.n + self.m + 2) / (self.m * self.p))
            self.R = self.r ** (self.p * self.k - 1)  # R = r^{kp-1}

        else:
            raise ValueError(f'Unsupported realization!')

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
        # refer: [3] p.19
        ci = self.R % self.N  # c0 = R
        wn = self.R.bit_length() - 1
        for _ in range(wn):
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

    def config(self, mul_opt='real1'):

        if self.mul_opt == mul_opt:
            return self

        self.mul_opt = mul_opt

        if mul_opt == 'real1':
            self.mul = self.real1_FPR2tN
        elif mul_opt == 'real2':
            self.mul = self.real2_FPR2t1
        elif mul_opt == 'real3':
            self.mul = self.real3_MWR2t1
        elif mul_opt == 'real4':
            self.mul = self.real4_FPR2tm
        elif mul_opt == 'real5':
            self.mul = self.real5_FPR2tm_v1
        elif mul_opt == 'real6':
            self.mul = self.real6_FPR2tm_v2
        elif mul_opt == 'real7':
            self.mul = self.real7_FPR2tm_v3
        elif mul_opt == 'real8':
            self.mul = self.real8_MWR2tm
        else:
            raise ValueError("Unsupported configuration mode for multiplication.")

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

        n = self.R.bit_length() - 1
        k = u * self.N_ & self.__RMASK
        t = (u + k * self.N) >> n

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
        n, p = self.R.bit_length() - 1, self.N

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
        n, p = self.R.bit_length() - 1, self.N

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
        X, Y, M, N = a, b, self.N, self.n - 1

        Z = 0
        for i in range(N + 1):  # b.c. self.n = N+1 where R=2**(N+1) for N defined in the paper
            Zi = Z + X * ith_bit(Y, i)

            q = Zi & 1  # reason for q = Zi*M_%b, here since b=2, M_=1, b.c. ext(x, 2)=xx, 1 always holds.
            Zi = (Zi + q * M) >> 1
            Z = Zi

        return Z

    def real3_MWR2t1(self, a: int, b: int) -> int:
        '''
            Scalable Montgomery Modular Multiplication With (Multiplier) Radix 2 (MWR2t1)
            refer: [2] Realization 3
        '''

        X, Y, M, N, w, e = a, b, self.N, self.n, self.w, self.e
        wmask = (1 << w) - 1

        _Z = 0
        for i in range(N + 1):
            Ca, Cb = 0, 0
            Yi = ith_bit(Y, i)

            q = 0
            for j in range(e):
                Xj, Zj = ith_word(X, j, w), ith_word(_Z, j + 1, w)
                # line 5:
                cazj = Zj + Xj * Yi + Ca
                Ca, Zj = cazj >> w, cazj & wmask
                _Z = update_ith_word(_Z, j + 1, w, Zj)

                # line 6-8:
                if j == 0:
                    q = ith_word(_Z, 0 + 1, w) & 1

                # line 9:
                Mj = ith_word(M, j, w)
                tmp = Zj + q * Mj + Cb
                Cb, Zj = tmp >> w, tmp & wmask
                _Z = update_ith_word(_Z, j + 1, w, Zj)

                # line 10:
                v = concatenate(Zj & 1, ith_word(_Z, (j - 1) + 1, w) >> 1, w - 1)
                _Z = update_ith_word(_Z, (j - 1) + 1, w, v)

            _Z = update_ith_word(_Z, (e - 1) + 1, w, 0)
            _Z = update_ith_word(_Z, (-1) + 1, w, 0)

        return _Z >> w

    def real4_FPR2tm(self, a: int, b: int) -> int:
        '''
            Fixed-Precision MM multiplication with Radix 2**m
            refer: [2] Realization 4
        '''
        X, Y, M, M_, n_m, m = a, b, self.N, self.N_, self.n_m, self.m

        Z = 0
        r = 1 << m
        mask = r - 1
        for i in range(n_m):
            Yi = Y & mask
            Zi = Z + X * Yi
            q = (Zi & mask) * M_ & mask
            Zi = (Zi + q * M) >> m
            Z = Zi

            Y >>= m

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

    def real6_FPR2tm_v2(self, a: int, b: int) -> int:
        '''
            Modified FPR2tm Version 2
            refer: [2] Realization 6
        '''
        X, Y, M, M_, m, d, w = a, b, self.N, self.N_, self.m, self.d, self.w
        r = 1 << m
        rmask = r - 1

        ZM2, ZM1, ZR1 = 0, 0, 0

        for i in range(d):
            Yi = ith_word(Y, i, m)
            qi = (((ZM2 >> m) + ZR1) & rmask) * M_ & rmask

            Zi = (ZM2 >> m) + ZR1 + qi * M + (X * Yi << m)

            ZMi, ZRi = decompose(Zi >> m, w, m)
            ZM2 = ZM1
            ZM1 = ZMi
            ZR1 = ZRi

        Z = (ZM2 >> m) + ZM1 + ZR1
        return Z

    def real7_FPR2tm_v3(self, a: int, b: int) -> int:
        '''
            Modified FPR2tm Version 3
            refer: [2] Realization 7
        '''
        X, Y, M, M_, m, d, N, w = a, b, self.N, self.N_, self.m, self.d, self.n, self.w
        r = 1 << m
        rmask = r - 1

        # line 1:
        ZSM2, ZSM1, ZSR1 = 0, 0, 0

        # line 2:
        ZCM2, ZCM1, ZCR1 = 0, 0, 0

        # line 3:
        c1 = 0

        for i in range(d):
            # line 5, 6
            TSi, TCi = compress([ZSM2 >> m, ZCM2 >> m, ZSR1, ZCR1, c1])
            qSi, qCi = compress([(TSi & rmask) * M_, (TCi & rmask) * M_])

            qi = (qSi + qCi) & rmask

            # line 8:
            Yi = ith_word(Y, i, m)
            ZSi, ZCi = compress([TSi, TCi, qi * M, (X * Yi) << m])

            # line 9:
            ci = 1 if (ZCi & rmask) != 0 else 0

            # line 10:
            ZSMi, ZSRi = decompose(ZSi >> m, w, m)

            # line 11:
            ZCMi, ZCRi = decompose(ZCi >> m, w, m)

            # update rolling array
            ZSM2 = ZSM1
            ZSM1 = ZSMi
            ZSR1 = ZSRi

            ZCM2 = ZCM1
            ZCM1 = ZCMi
            ZCR1 = ZCRi

            c1 = ci

        # line 13
        Z = (ZSM2 >> m) + ZSM1 + ZSR1 + (ZCM2 >> m) + ZCM1 + ZCR1 + c1

        return self.correction(Z)

    def real8_MWR2tm(self, a: int, b: int) -> int:
        X, Y, M, M_, m, k, p, e, w = a, b, self.N, self.N_, self.m, self.k, self.p, self.e, self.w
        r = 1 << m
        rmask = r - 1
        wmask = (1 << w) - 1

        X_ = X << m  # X'=Xr

        _ZSM_, _ZSM, _ZSR = 0, 0, 0
        _ZCM_, _ZCM, _ZCR = 0, 0, 0
        c = 0

        for i in range(k * p):  # outer loop

            FBS = FBC = 0
            Yi = ith_word(Y, i, m)

            # print(f'Outer-loop:{i}\n')
            # print('FBS:     ', hex(FBS))
            # print('FBC:     ', hex(FBC))
            # print('ZSM\':   ', hex(_ZSM_ >> w))
            # print('ZCM\':   ', hex(_ZCM_ >> w))
            # print('ZSM:     ', hex(_ZSM >> w))
            # print('ZCM:     ', hex(_ZCM >> w))
            # print('ZSR:     ', hex(_ZSR >> w - m))
            # print('ZCR:     ', hex(_ZCR >> w - m))
            # print('Yi:      ', hex(Yi))
            # print('M\':     ', hex(M_))
            # print('c:       ', hex(c))
            # print("\n")

            q = 0
            for j in range(e):  # inner loop

                # line 6:
                ZSRj, ZCRj = ith_word(_ZSR, j + 1, w - m), ith_word(_ZCR, j + 1, w - m)  # ZSR_j, ZCR_j
                ZSM_j, ZCM_j = ith_word(_ZSM_, j + 1, w), ith_word(_ZCM_, j + 1, w)  # ZSM'_j, ZCM'_j
                TS, TC = compress([ZSRj, ZCRj, ZSM_j >> m, ZCM_j >> m, c if j == 0 else 0])

                # line 7-10:
                if j == 0:
                    qS, qC = compress([(TS & rmask) * M_, (TC & rmask) * M_])
                    q = (qS + qC) & rmask

                # line 11:
                X_j, Mj = ith_word(X_, j, w), ith_word(M, j, w)  # X'_j, M_j
                # OS, OC = compress([TS, TC, X_j * Yi, q * Mj, FBS, FBC])  # naive approach fails!
                OS, OC = compressor_2m4to2(TS, TC, X_j, Yi, q, Mj, FBS, FBC, w, m)

                # line 12:
                c = 1 if (j == 0 and (OC & rmask) != 0) else 0

                # line 13-16:
                _ZSM_ = update_ith_word(_ZSM_, j - 1 + 1, w, ith_word(_ZSM, j - 1 + 1, w))  # ZSM'_{j-1}
                _ZCM_ = update_ith_word(_ZCM_, j - 1 + 1, w, ith_word(_ZCM, j - 1 + 1, w))  # ZCM'_{j-1}
                _ZSM = update_ith_word(_ZSM, j - 1 + 1, w, (OS & rmask) << w - m)  # ZSM_{j-1}
                _ZCM = update_ith_word(_ZCM, j - 1 + 1, w, (OC & rmask) << w - m)  # ZCM_{j-1}

                # line 17, 18:
                _ZSR = update_ith_word(_ZSR, j + 1, w - m, extract_bits(OS, m, w - 1))  # ZSR_{j}
                _ZCR = update_ith_word(_ZCR, j + 1, w - m, extract_bits(OC, m, w - 1))  # ZCR_{j}

                # line 19, 20:
                FBS, FBC = extract_bits(OS, w, w + m + 1), extract_bits(OC, w, w + m + 1)

            # line 22:
            _ZSM = update_ith_word(_ZSM, e - 1 + 1, w, 0)
            _ZSM = update_ith_word(_ZSM, -1 + 1, w, 0)

            # line 23:
            _ZCM = update_ith_word(_ZCM, e - 1 + 1, w, 0)
            _ZCM = update_ith_word(_ZCM, -1 + 1, w, 0)

            # print('Outter loop results: \n')
            # print('ZSM\':   ', hex(_ZSM_ >> w))
            # print('ZCM\':   ', hex(_ZCM_ >> w))
            # print('ZSM:     ', hex(_ZSM >> w))
            # print('ZCM:     ', hex(_ZCM >> w))
            # print('ZSR:     ', hex(_ZSR >> w - m))
            # print('ZCR:     ', hex(_ZCR >> w - m))
            # print("\n")

        # line 25, 26
        DO1 = DO2 = 0
        Carry = c

        Z = 0
        ZSM, ZCM = _ZSM >> w, _ZCM >> w
        ZSM_, ZCM_ = _ZSM_ >> w, _ZCM_ >> w
        ZSR, ZCR = _ZSR >> w - m, _ZCR >> w - m
        for i in range(e):
            # line 28
            ZSMi, ZCMi = ith_word(ZSM, i, w), ith_word(ZCM, i, w)
            ZSRi, ZCRi = ith_word(ZSR, i, w - m), ith_word(ZCR, i, w - m)
            ZSM_i, ZCM_i = ith_word(ZSM_, i, w), ith_word(ZCM_, i, w)
            PS, PC = compress([ZSMi, ZCMi,
                               ZSRi, ZCRi,
                               ZSM_i >> m, ZCM_i >> m,
                               DO1, DO2])

            # line 29
            CarryZi = (PS & wmask) + (PC & wmask) + Carry
            Carry = CarryZi >> w
            Zi = CarryZi & wmask

            # line 30
            DO1, DO2 = (PS >> w) & 1, (PC >> w) & 1

            # update Z
            Z += Zi << (i * w)

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


# todo> remove: quick tb
if __name__ == '__main__':
    # region sample usage
    # M = Montgomery.factory(mod=31, mul_opt='real5').build(m=4)
    p = 61
    # M = Montgomery.factory(mod=p, mul_opt='real3').build(w=8)

    M = Montgomery.factory(mod=p, mul_opt='real8').build(m=2, w=5)
    # M = Montgomery.factory(mod=p, mul_opt='real7').build(m=2, w=5)
    # M = Montgomery.factory(mod=p, mul_opt='real6').build(m=2, w=5)
    x = 46

    _x = M(x)

    print('done')

    # x_, y_ = 7, 21  # x', y' as elements in prime field F_{31}
    #
    # # enter the Montgomery domain over F_{31}
    # x, y = M(7), M(21)
    # a = x + y
    # b = x - y
    # c = x * y
    # d = 1 / x
    # e = x / y
    #
    # # exit the Montomery domain
    # _x, _y = int(x), int(y)
    # assert _x == x_
    # assert _y == y_
    #endregion

#     mont = Montgomery.factory(mod=257, mul_opt='real8').build(m=8, w=16)
#
#     x = 128
#     x_mont = mont(x)
#
#     x_mont_exp = x * mont.R % mont.N
#     assert x_mont.value == x_mont_exp
