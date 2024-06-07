## Montgomery domain over F_p, we consider p as a prime other than an arbitrary odd number here
# so that Montgomery domain will be a field in this case.
## @author: Luke Li

class Montgomery:
    def __init__(self, mod):
        if not mod & 1:
            raise ValueError("The modulus must be an odd number.")
        self.N = mod
        self.__pow = mod.bit_length()
        self.R = 1 << self.__pow  # Set R as the smallest power of 2 greater than N by default
        self.__MASK = self.R - 1

        self.config_mode = 'real1'  # Default to 'real1'
        self.w = self.__pow
        self.r = 1 << self.w  # r = 2^w, by default r==R
        self.__pre_calc()  # do some pre-calculation

    def __call__(self, v):
        v_ = self.__enter_domain(v)
        return MontgomeryNumber(v_, self)

    def __pre_calc(self):
        self.R2 = self.__pre_calc_R2()  # R^2 % N
        self.N_ = self.__pre_calc_N_()  # N' for rr'-NN'=1

    def __enter_domain(self, a):
        return self.REDC(a * self.R2)

    def __pre_calc_R2(self):
        # Use (wn + 1) rounds of modulo addition to get R^2 % N, where R = 2^{wn}
        # refer: <<Topics>> p.19
        ci = self.R % self.N  # c0 = R
        for _ in range(self.__pow):
            ci = (ci + ci) % self.N

        return ci

    def __pre_calc_N_(self):
        return self.mont_const(self.w)

    def mont_const(self, w):
        '''
            a.k.a Montgomery constant N', for rr'-NN'=1, where r=2^w
            refer: <<Topics>> p.18, Alg. 2.3
        :return:
        '''
        y = 1
        for i in range(2, w + 1):
            if y * self.N & ((1 << i) - 1) == 1:
                continue

            y += (1 << i - 1)

        return (1 << w) - y  # r - y

    def config(self, mode):
        if mode not in ['real1', 'real2']:
            raise ValueError("Unsupported configuration mode.")

        if self.config_mode == mode:
            return

        self.config_mode = mode
        self.__pre_calc()

    def REDC(self, u):
        '''
            Montgomery reduction (REDC) of number u
            ref: https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
        :param u:
        :return: u * R^{-1} mod N
        '''
        if u < 0 or u > self.R * self.N - 1:
            raise ValueError('Given number if out of montgomery reduction range')

        k = u * self.N_ & self.__MASK
        t = (u + k * self.N) >> self.__pow

        return self.correction(t)

    def correction(self, r):
        if r > 2 * self.N:
            raise ValueError(f'Only number is ready for Montgomery correction step is allowed, check input')

        return r if r < self.N else r - self.N

    def __alm_inv(self, a):
        '''
            Kaliski's almost inversion algorithm
            refer: <<Montgomery inversion>> Alg. 4
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
        ref:  <<Montgomery inversion>> Alg. 5
        '''
        n, p = self.__pow, self.N

        j = 0
        while j < (n << 1) - k:
           r <<= 1
           if r > p:
               r -= p
           j += 1

        return r

    def __cor_phase_fb(self, r, k):
        '''corretion phase of forward/backward montgomery inversion
        ref:  <<Montgomery inversion>> Alg. 5
        '''
        n, p = self.__pow, self.N

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


class MontgomeryNumber:
    def __init__(self, value, mont: Montgomery):
        self.mont = mont
        self.value = value if value < mont.N else value % mont.N

    def __mul__(self, other):
        if isinstance(other, MontgomeryNumber):
            t = self.value * other.value  # todo> reals
            return MontgomeryNumber(self.mont.REDC(t), self.mont)
        elif isinstance(other, int):
            # other = self.mont(other)
            # !! Decision choice !!
            # i * mont(a) should treat i as the montgomery repr other than converting it automatically!!
            other = MontgomeryNumber(value=other, mont=self.mont)
            return self.__mul__(other)
        return NotImplemented

    def __rmul__(self, other: int):
        return self.__mul__(other)

    def __truediv__(self, other):
        #todo>
        if isinstance(other, MontgomeryNumber):
            inv = self.mont._mont_inv(other.value)
            return MontgomeryNumber((self * inv).value, self.mont)
        elif isinstance(other, int):
            # Design choice: if we are dividing a '1', we should treat it as the identity in Montgomery Field;
            # o.w. treat it an 'as-is' integer already in Montgomery repr.
            other = self.mont(other) if other == 1 else MontgomeryNumber(other, self.mont)
            return self.__truediv__(other)
        return NotImplemented

    def __rtruediv__(self, other):
        if not isinstance(other, int):
            return NotImplemented

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
        return self.mont.REDC(self.value)
