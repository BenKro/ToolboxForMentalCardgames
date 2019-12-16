# ! /usr/bin/env python
# -*- coding: utf-8 -*-
# ===================================================================
# random_generator.py
#
# 27.06.2019, Benjamin Kromer <benjamin.kromer1@web.de>
#
# @desc: Random numbers and permutations used for elliptic curve
#        cryptography and zero-knowledge proofs.
# ===================================================================
import random
import secrets
import hashlib


class RandomGenerator:
    """Class for different random values and permutations

    Attributes:
        order (int): order of the elliptic curve subgroup
    """

    def __init__(self, order):
        """
        Args:
            order (int): order of the elliptic curve subgroup
        """
        self.order = order

    def get_random_value(self):
        """Get one random value in range 1 to order- 1

        Returns:
            int: random value in range 1 to order - 1
        """
        return 1 + secrets.randbelow(self.order-1)

    @staticmethod
    def get_random_value_range(x, y):
        """Get one random value in range x to y-1

        Args:
            x (int): lower bound
            y (int): upper bound

        Returns:
            int: value from [x,y)
        """
        return x + secrets.randbelow(y-x)

    def get_random_array(self, size):
        """Get a list with size random values between 1 and order-1

        Args:
            size (int): number of random values

        Returns:
            List[int]: list with random value
        """
        var = []
        for i in range(size):
            var.append(self.get_random_value())
        return var

    def get_random_permutation(self, size, array=None):
        """Permute an array randomly with fisher-yates algorithm,
        if no array is given as argument, array = [0,1,...,size]

        Args:
            size (int): number of random values
            array (List): array to be permuted

        Returns:
            List: permuted array
        """
        if array is None:
            array = list(range(0, size))

        for i in range(size-1):
            j = self.get_random_value_range(i, size)
            var0 = array[i]
            array[i] = array[j]
            array[j] = var0

        return array

    @staticmethod
    def get_random_permutation_seed(size, seed, array=None):
        """Permute an array randomly by seed with fisher-yates algorithm,
        if no array is given as argument, array = [0,1,...,size]

        Args:
            size (int): number of random values
            seed (int): seed for PRNG
            array (List): array to be permuted

        Returns:
            List: permuted array
        """
        random.seed(seed)

        if array is None:
            array = list(range(0, size))

        for i in range(size-1):
            j = random.randint(i, size-1)
            var0 = array[i]
            array[i] = array[j]
            array[j] = var0

        return array

    @staticmethod
    def get_random_from_hash(i, j, number, *args):
        """Get random values in range i to j-1 from PRNG, while hash from
        *args is used as seed.

        Args:
            i (int): lower bound
            j (int): upper bound
            number (int): number of random values
            *args (int): input for hash function

        Returns:
            [List[int], int]: list with random values and seed from hash
        """
        var0 = hashlib.sha3_256()
        for arg in args:
            var1 = arg.bit_length()
            var2 = var1 // 8 + (var1 % 8 > 0)
            var3 = int.to_bytes(arg, var2, "big")
            var0.update(var3)

        seed = var0.hexdigest()
        random.seed(seed)
        var1 = []
        for k in range(number):
            var1.append(random.randint(i, j))
        return var1, seed
