import unittest

from mont.montgomery import Montgomery
from mont.common import *
import os


class TestMontgomeryOperations(unittest.TestCase):
    def setUp(self):
        self.FULL_TEST = False  # run all test cases? caution for large prime domains!!
        self.RADOMIZED = 1

        # NOTE: change this number should delete the m_primes_xxx.data first.
        num_primes = 10  # num of primes >= num of domains; change this shall delete the data file to regenerate primes

        self.MAX_EXP = 10000  # maximum power to be tested

        if self.RADOMIZED:
            filename = './m_primes_256.data'
            self.num_count = 3  # randomly generated repr. numbers in Z_p
        else:
            # fixed
            filename = './t_primes.data'
            self.num_count = 10

        if not os.path.exists(filename):
            # Generate prime numbers
            large_primes = generate_large_primes(num_primes, 256)
            save_to_file(large_primes, filename=filename)
        else:
            large_primes = load_primes_from_file(filename)

        self.rands = random_list(low=2, high=max(large_primes), count=self.num_count)
        self.monts = []
        mont0 = Montgomery.factory(mod=large_primes[0], mul_opt='real0').build()
        self.monts.append(mont0)

        if self.RADOMIZED:
            # real-1
            mont1 = Montgomery.factory(mod=large_primes[1], mul_opt='real1').build()
            self.monts.append(mont1)

            # real-2
            mont2 = Montgomery.factory(mod=large_primes[2]).build()
            mont2.config(mul_opt='real2').build(n=nextp2(mont2.N).bit_length() + 2)
            self.monts.append(mont2)

            # real-3
            mont3 = Montgomery.factory(mod=large_primes[3], mul_opt='real3').build(w=5)
            # self.monts[3] = mont3
            self.monts.append(mont3)

            # real-4
            mont4 = Montgomery.factory(mod=large_primes[4], mul_opt='real4').build(m=4)
            self.monts.append(mont4)

            # real-5
            mont5 = Montgomery.factory(mod=large_primes[5], mul_opt='real5').build(m=5)
            self.monts.append(mont5)

            # real-6
            mont6 = Montgomery.factory(mod=large_primes[6], mul_opt='real6').build(m=6)
            self.monts.append(mont6)

            mont7 = Montgomery.factory(mod=large_primes[7], mul_opt='real7').build(m=4, w=8)
            self.monts.append(mont7)

            mont8 = Montgomery.factory(mod=large_primes[0], mul_opt='real8').build(m=64, w=256)
            self.monts.append(mont8)

        else:
            # customized fixed case here
            self.monts[0].config(mul_opt='real8').build(m=2, w=5)


    def test_domain(self):

        _FULL_TEST = self.FULL_TEST
        try:
            with timeout(120):  # 2 minutes

                cnt = 0
                total = sum(mont.N for mont in self.monts) if _FULL_TEST else self.num_count * len(self.monts)
                for mont in self.monts:
                    # verify pre-calc
                    R2_exp = (mont.R % mont.N) ** 2 % mont.N
                    self.assertEqual(R2_exp, mont.R2)


                    draws = range(mont.N) if _FULL_TEST else self.rands
                    # verify enter/exit domain
                    for x in draws:
                        x_mont = mont(x)
                        x_mont_exp = x * mont.R % mont.N
                        self.assertEqual(x_mont_exp, x_mont.value, f'Montgomery repr of x={x} failed!')

                        x_ = int(x_mont)
                        self.assertEqual(x % mont.N, x_)

                        print(f'{cnt + 1}/{total}% : Test domain on {mont} succeeded!')
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

    # def test_compressor(self):
    #     for mont in self.monts:
    #
    #         p = mont.N
    #         for _ in range(100000):
    #             # Generate three random integers
    #             num1 = randint(2, p)
    #             num2 = randint(2, p)
    #             num3 = randint(3, p)
    #
    #             exp = num1 + num2 + num3
    #             a, b = compress([num1, num2, num3])
    #             self.assertEqual(exp, a + b)

    def test_multiplication(self):

        cnt, total = 0, (self.num_count - 1) * len(self.monts)
        for i in range(self.num_count - 1):

            a, b = self.rands[i], self.rands[i + 1]

            for mont in self.monts:
                p = mont.N
                a, b = a % p, b % p
                a_, b_ = mont(a), mont(b)

                act = a_ * b_
                act1 = a * b_
                act2 = a_ * b

                act3 = a_ + a_
                act4 = b_ + b_

                # Direct addition and reduction
                exp = mont((a * b) % p)
                exp1 = mont(mont.exit_domain(a) * b % p)
                exp2 = mont(a * mont.exit_domain(b) % p)

                exp3 = mont((2*a) % p)
                exp4 = mont((b*2) % p)

                self.assertEqual(act, exp)
                self.assertEqual(act1, exp1)
                self.assertEqual(act2, exp2)

                self.assertEqual(act3, exp3)
                self.assertEqual(act4, exp4)

                print(f'{cnt + 1}/{total}% : Test multiplication on {mont} succeeded!')
                cnt += 1

    def test_exponentiation(self):
        _FULL_TEST = self.FULL_TEST

        test_per_mont = self.num_count if _FULL_TEST else 1
        cnt, total = 0, test_per_mont * len(self.monts)
        for i in range(test_per_mont):

            x = self.rands[i]
            e = random.randint(1, self.MAX_EXP)
            e1 = (1 << i + 1) + 1  # particularly covers e = 2^k + 1 cases, for k >= 1

            # corner cases
            corners = [1, 2]

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

                for ec in corners:
                    act2 = x_ ** ec
                    exp2 = mont((x**ec) % p)
                    self.assertEqual(act2, exp2)

                print(f'{cnt + 1}/{total}% : Test exponentiation on {mont} succeeded!')
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
        print(f'Test additions all succeeded!')

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

        print(f'Test subtraction all succeeded!')

    def test_almost_inversion(self):

        for mont in self.monts:
            p = mont.N

            for a in self.rands:
                a_ = bin_inv(a, p)
                act, k = mont._Montgomery__alm_inv(a)
                exp = a_ * (1 << k) % p
                self.assertEqual(act, exp)

        print(f'Test (almost) inversion all succeeded!')

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

        print(f'Test inversion all succeeded!')

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

        print(f'Test division all succeeded!')


# Rerun the corrected unit tests with refined calculations
unittest.main(argv=[''], verbosity=2, exit=False)
