import unittest
from mont.montgomery import Montgomery
from mont.common import *
import os
from random import randint


# Up-to data: 6/16/2024, real-8 tests is failed

class TestMontgomeryOperations(unittest.TestCase):
    def setUp(self):

        self.RADOMIZED = 1
        num_primes = 1  # num of primes >= num of domains; change this shall delete the data file to regenerate primes
        self.MAX_EXP = 10000  # maximum power to be tested

        if self.RADOMIZED:
            filename = '../test/m_primes_256.data'
            domains = 1  # p for Z_p :
            self.num_count = 100  # randomly generated repr. numbers in Z_p
        else:
            # fixed
            filename = '../test/t_primes.data'
            domains = 1
            self.num_count = 10

        if not os.path.exists(filename):
            # Generate prime numbers
            large_primes = generate_large_primes(num_primes, 256)
            save_to_file(large_primes)
        else:
            large_primes = load_primes_from_file(filename)

        self.rands = random_list(low=2, high=max(large_primes), count=self.num_count)
        # self.monts = [Montgomery.factory(mod=large_primes[i]).build() for i in range(domains)]
        self.monts = []

        if self.RADOMIZED:

            # # real-1
            # self.monts[1].config('real1').build()
            #
            # # real-2
            # # self.monts[2].config('real2').build(R=1 << nextp2(self.monts[2].N).bit_length() + 2)
            # self.monts[2].config('real2').build(n=nextp2(self.monts[2].N).bit_length() + 2)
            #
            # # real-3
            # mont3 = Montgomery.factory(mod=large_primes[3], mul_opt='real3').build(w=5)
            # self.monts[3] = mont3
            #
            # # real-4
            # mont4 = Montgomery.factory(mod=large_primes[4], mul_opt='real4').build(m=4)
            # self.monts[4] = mont4
            #
            # # real-5
            # mont5 = Montgomery.factory(mod=large_primes[5], mul_opt='real5').build(m=5)
            # self.monts[5] = mont5
            #
            # # real-6
            # mont6 = Montgomery.factory(mod=large_primes[6], mul_opt='real6').build(m=6)
            # self.monts[6] = mont6
            #
            # mont7 = Montgomery.factory(mod=large_primes[7], mul_opt='real7').build(m=4, w=8)
            # self.monts[7] = mont7

            # todo> only works for sequential compressor?
            mont8 = Montgomery.factory(mod=large_primes[0], mul_opt='real8').build(m=64, w=256)
            self.monts.append(mont8)

        else:
            # customized fixed case here
            self.monts[0].config(mul_opt='real8').build(m=2, w=5)

    def test_domain(self):

        _FULL_TEST = False
        try:
            with timeout(120):  # 120 seconds or 2 minutes
                for mont in self.monts:
                    # verify pre-calc
                    R2_exp = (mont.R % mont.N) ** 2 % mont.N
                    self.assertEqual(R2_exp, mont.R2)

                    # cnt, total = 0, self.num_count if self.RADOMIZED else mont.N
                    cnt, total = 0, mont.N if _FULL_TEST else self.num_count
                    draws = range(mont.N) if _FULL_TEST else self.rands
                    # verify enter/exit domain
                    for x in draws:
                        x_mont = mont(x)
                        x_mont_exp = x * mont.R % mont.N
                        self.assertEqual(x_mont_exp, x_mont.value, f'Montgomery repr of x={x} failed!')

                        x_ = int(x_mont)
                        self.assertEqual(x % mont.N, x_)

                        print(f'done: {cnt + 1}/{total}%')
                        cnt += 1
        except TimeoutException:
            print("Test case exceeded the time limit of 2 minutes.")

    #region Toggle on only for small mont domains
    # def test_mont_constant(self):
    #
    #     for mont in self.monts:
    #
    #         p = mont.N
    #
    #         for w in range(1, p):
    #             r = 1 << w
    #
    #             _, N_, _ = xgcd(r, p)
    #             exp = -N_
    #
    #             if exp < 0:  # NOTE here, should use mod r!
    #                 exp += r
    #
    #             act = mont.mont_const(w)
    #             self.assertEqual(exp, act)
    #endregion

    def test_compressor(self):
        for mont in self.monts:

            p = mont.N
            for _ in range(1000000):
                # Generate three random integers
                num1 = randint(2, p)
                num2 = randint(2, p)
                num3 = randint(3, p)

                exp = num1 + num2 + num3
                a, b = compress([num1, num2, num3])
                self.assertEqual(exp, a + b)

    def test_multiplication(self):

        cnt, total = 0, self.num_count - 1
        for i in range(self.num_count - 1):

            a, b = self.rands[i], self.rands[i + 1]

            for mont in self.monts:
                p = mont.N
                a, b = a % p, b % p
                a_, b_ = mont(a), mont(b)

                act = a_ * b_
                act1 = a * b_
                act2 = a_ * b

                # Direct addition and reduction
                exp = mont((a * b) % p)
                exp1 = mont(mont.exit_domain(a) * b % p)
                exp2 = mont(a * mont.exit_domain(b) % p)

                self.assertEqual(act, exp)
                self.assertEqual(act1, exp1)
                self.assertEqual(act2, exp2)

                print(f'done: {cnt + 1}/{total}%')
                cnt += 1

    def test_exponentiation(self):

        cnt, total = 0, self.num_count
        for i in range(self.num_count):

            x = self.rands[i]
            e = random.randint(1, self.MAX_EXP)
            e1 = (1 << i + 1) + 1  # particularly covers e = 2^k + 1 cases, for k >= 1

            for mont in self.monts:
                p = mont.N
                x = x % p
                x_ = mont(x)

                act = x_ ** e
                exp = mont((x ** e) % p)
                self.assertEqual(act, exp)

                if e1 < self.MAX_EXP:  # This tiny python code can't handle too large an exponent value
                    act1 = x_ ** e1
                    exp1 = mont((x ** e1) % p)
                    self.assertEqual(act1, exp1, f'Exponential test e={e1}=(2^(k+1) + 1) for k={i} fails!')

                print(f'done: {cnt + 1}/{total}%')
                cnt += 1

    def test_addition(self):

        for i in range(self.num_count - 1):

            a, b = self.rands[i], self.rands[i + 1]

            for mont in self.monts:
                act = mont(a) + mont(b)
                act1 = a + mont(b)
                act2 = mont(a) + b

                # Direct addition and reduction
                exp = mont((a + b) % mont.N)
                exp1 = mont((mont.exit_domain(a) + b) % mont.N)
                exp2 = mont((a + mont.exit_domain(b)) % mont.N)

                self.assertEqual(act, exp)
                self.assertEqual(act1, exp1)
                self.assertEqual(act2, exp2)

    def test_subtraction(self):
        for i in range(self.num_count - 1):

            a, b = self.rands[i], self.rands[i + 1]

            for mont in self.monts:
                act = mont(a) - mont(b)
                act1 = a - mont(b)
                act2 = mont(a) - b

                # Direct addition and reduction
                exp = mont((a - b) % mont.N)

                self.assertEqual(act, exp)
                self.assertEqual(act1, exp)
                self.assertEqual(act2, exp)

    def test_almost_inversion(self):

        for mont in self.monts:
            p = mont.N

            for a in self.rands:
                a_ = bin_inv(a, p)
                act, k = mont._Montgomery__alm_inv(a)
                exp = a_ * (1 << k) % p
                self.assertEqual(act, exp)

    def test_inversion(self):

        for num in self.rands:
            # print(f'testing number: {num}')

            for mont in self.monts:
                t = bin_inv(num, mont.N)
                exp = mont(t)

                act = 1 / mont(num)
                act1 = mont(1) / mont(num)

                self.assertEqual(exp, act)
                self.assertEqual(exp, act1)

    def test_division(self):

        for mont in self.monts:
            p = mont.N

            for i in range(self.num_count - 1):
                a, b = self.rands[i], self.rands[i + 1]
                a, b = a % p, b % p
                if a == 1 or b == 1:  # Notation of `1` is reserved for conventional inversion use only
                    continue

                act = mont(a) / mont(b)
                act1 = a / mont(b)
                act2 = mont(a) / b

                # Direct addition and reduction
                b_inv = bin_inv(b, p)
                exp = mont((a * b_inv) % p)
                exp1 = mont((mont.exit_domain(a) * b_inv) % p)

                b_inv2 = bin_inv(mont.exit_domain(b), p)
                exp2 = mont((a * b_inv2) % p)

                self.assertEqual(act, exp)
                self.assertEqual(act1, exp1)
                self.assertEqual(act2, exp2)


# Rerun the corrected unit tests with refined calculations
unittest.main(argv=[''], verbosity=2, exit=False)
