import unittest
from mont.common import *


class TestCommon(unittest.TestCase):

    def test_random_lcm(self):
        """Generates random pairs of integers with random bit lengths and tests the LCM implementation."""
        n = 10000  # number of random integers
        min_bit_length = 64  # minimum bit length of integers
        max_bit_length = 2048  # maximum bit length of integers
        random_integers = [random.getrandbits(random.randint(min_bit_length, max_bit_length)) for _ in range(n)]

        # Generate pairs and test each
        for i in range(0, n, 2):  # Assuming n is even for simplicity
            a = random_integers[i]
            b = random_integers[i + 1]
            _, _, d = xgcd(a, b)
            expected_lcm = abs(a * b) // d
            calculated_lcm = lcm(a, b)
            self.assertEqual(calculated_lcm, expected_lcm, f"Failed for pair ({a}, {b})")

    def test_addition_chain(self):
        acs = AdditionChains()

        e = 45
        chain = acs.find_chain(e)

        # e_rsa = 65537
        # chain_rsa = acs.find_chain(e_rsa)
        steps = [x.bit_length() for x in chain]

        print(steps)


unittest.main(argv=[''], verbosity=2, exit=False)
