import unittest
from src.mont import Montgomery
from common import *


class TestMontgomeryOperations(unittest.TestCase):
    def setUp(self):

        RADOMIZED = 1

        if RADOMIZED:
            filename = 'm_primes.data'
            num_primes = 5  # p for Z_p : change this shall delete the data file to regenerate primes
            self.num_count = 100  # randomly generated repr. numbers in Z_p
        else:
            # fixed
            filename = 't_primes.data'
            num_primes = 1
            self.num_count = 10

        if not os.path.exists(filename):
            # Generate prime numbers
            large_primes = generate_large_primes(num_primes, 1024)
            save_to_file(large_primes)
        else:
            large_primes = load_primes_from_file(filename)

        self.monts = [Montgomery(mod=large_primes[i]) for i in range(num_primes)]
        self.rands = random_list(low=2, high=max(large_primes), count=self.num_count)

    def test_domain(self):

        for mont in self.monts:
            # verify pre-calc
            R2_exp = (mont.R % mont.N) ** 2 % mont.N
            self.assertEqual(R2_exp, mont.R2)

            # verify enter/exit domain
            for x in self.rands:
                x_mont = mont(x)
                x_mont_exp = x * mont.R % mont.N
                self.assertEqual(x_mont.value, x_mont_exp)

                x_ = int(x_mont)
                self.assertEqual(x_, x % mont.N)

    def test_multiplication(self):

        for i in range(self.num_count - 1):

            a, b = self.rands[i], self.rands[i + 1]

            for mont in self.monts:
                p = mont.N
                a, b = a % p, b % p

                act = mont(a) * mont(b)
                act1 = a * mont(b)
                act2 = mont(a) * b

                # Direct addition and reduction
                exp = mont((a * b) % p)
                exp1 = mont(mont.REDC(a) * b % p)
                exp2 = mont(a * mont.REDC(b) % p)


                self.assertEqual(act, exp)
                self.assertEqual(act1, exp1)
                self.assertEqual(act2, exp2)

    def test_addition(self):

        for i in range(self.num_count - 1):

            a, b = self.rands[i], self.rands[i + 1]

            for mont in self.monts:
                act = mont(a) + mont(b)
                act1 = a + mont(b)
                act2 = mont(a) + b

                # Direct addition and reduction
                exp = mont((a + b) % mont.N)
                exp1 = mont((mont.REDC(a) + b) % mont.N)
                exp2 = mont((a + mont.REDC(b)) % mont.N)

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
        # todo> remove
        # p = 61
        # a = 1
        # mont = Montgomery(mod=p)
        # a_mont = mont(a)
        # a_, _, _ = xgcd(a, p)
        # t = bin_inv(a, p)
        # exp = mont(t)
        # act = 1 / a_mont

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
                exp1 = mont((mont.REDC(a) * b_inv) % p)

                b_inv2 = bin_inv(mont.REDC(b), p)
                exp2 = mont((a * b_inv2) % p)

                self.assertEqual(act, exp)
                self.assertEqual(act1, exp1)
                self.assertEqual(act2, exp2)



# Rerun the corrected unit tests with refined calculations
unittest.main(argv=[''], verbosity=2, exit=False)
