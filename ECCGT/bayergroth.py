# ! /usr/bin/env python
# -*- coding: utf-8 -*-
# ===================================================================
# bayergroth.py
#
# 10.09.2019, Benjamin Kromer <benjamin.kromer1@web.de>
#
# @desc: Zero-Knowledge Argument for Correctness of a Shuffle by Stephanie
#        Bayer and Jens Groth to prove that: ciphers_out[i] = ciphers_in[pi[
#        i]] + Enc_pk(O, rho[i]).
# ===================================================================
from copy import deepcopy
from ECCGT.eccwrapper import Fastecdsa as ECCobj
from ECCGT.eccwrapper import ShortPoint


class BayGroProver:
    """Prover in Zero-Knowledge Argument for Correctness of a Shuffle such
    that ciphers_out[i] = ciphers_in[pi[i]] + Enc_pk(O, rho[i])

    """
    def __init__(self, N, m, n, curve, pi, gen, pk, ciphers_in, ciphers_out,
                 rho):
        """
        Args:
            N (int): number of cards N=m*n
            m (int): rows
            n (int): columns
            curve (ECCObj): curve object
            pi (List[int]): permutation
            gen (List[ShortPoint]): generators for pedersen commitment
            pk (ShortPoint): public key
            ciphers_in (List[List[ShortPoint]]): ciphers before the
                permutation and re-masking
            ciphers_out (List[List[ShortPoint]]): ciphers after the
                permutation and re-masking
            rho (List[int]): random parameters for re-masking
        """
        self.m = m
        self.n = n
        self.N = N
        assert (m * n == N)

        self.curve = curve
        self.pubKey = deepcopy(pk)

        self.ciphers_in = deepcopy(ciphers_in)
        self.ciphers_out = deepcopy(ciphers_out)
        self.rho = rho

        self.generators_ck = gen
        self.pedersen = Pedersen(gen, n, self.curve)
        self.order = self.curve.order

        # R1 (Shuffle)
        self.a = [pi[i]+1 for i in range(len(pi))]
        self.A, self.r_A, self.c_A = (None,)*3

        # R3 (Shuffle)
        self.x2 = None
        self.B, self.r_B, self.c_B = (None,)*3

        # R1 (Hadamard/Zero)
        self.y4, self.z4 = (None,)*2
        self.D = None
        self.F, self.r_F, self.c_F = (None,)*3
        self.G, self.r_G, self.c_G = (None,)*3
        self.z, self.r_z, self.c_z = (None,)*3

        # R3 (Hadamard/Zero)
        self.x6, self.y6, self.x6_array = (None,)*3
        self.H, self.r_H = (None,) * 2
        self.H_m, self.r_H_m, self.c_H_m = (None,) * 3
        self.F_0, self.r_F_0, self.c_F_0 = (None,)*3
        self.P, self.r_P, self.c_P = (None,) * 3

        # R5 (Hadamard/Zero)
        self.x8 = None
        self.f, self.r_f = (None,)*2
        self.h, self.r_h = (None,) * 2
        self.r_p = None

        # R1 (Single Value Product)
        self.g = None
        self.alpha = None
        self.delta, self.s_delta, self.c_delta = (None,) * 3
        self.gamma, self.r_gamma, self.c_gamma = (None,) * 3
        self.s_Delta, self.c_Delta = (None,)*2

        # R3 (Single Value Product)
        self.gamma_tilde, self.r_gamma_tilde = (None,) * 2
        self.alpha_tilde, self.r_alpha_tilde = (None,) * 2

        # R1 (Multi-Exponent)
        self.b_0, self.r_B0, self.c_B0 = (None,) * 3
        self.beta, self.r_beta, self.c_beta = (None,) * 3
        self.E, self.tau_k = (None,)*2

        # R3 (Multi-Exponent)
        self.b, self.r_b = (None,) * 2
        self.beta_tilde, self.r_beta_tilde = (None,) * 2
        self.tau = None

    def r1_shuffle(self):
        """Form array a to Matrix A and create Pedersen commitment c_A with
        m random values r_A

        Returns:
            List[ShortPoint]: commitment to A
        """
        self.r_A = self.curve.rand_gen.get_random_array(self.m)

        # copy array a of size N into matrix A of size mxn
        self.A = [self.a[self.n*i: self.n*(i+1)] for i in range(self.m)]

        # pedersen commitment of A
        self.c_A = self.pedersen.commit_matrix_vector(self.A, self.r_A)

        return self.c_A

    def r3_shuffle(self, x2):
        """Create b from a, form it to Matrix B and create Pedersen
        commitment c_B with m random values r_A

        Args:
            x2 (int): challenge from verfifier

        Returns:
            List[ShortPoint]: commitment to B
        """
        # save challenge
        self.x2 = x2

        self.r_B = self.curve.rand_gen.get_random_array(self.m)

        # set array bi as x^ai
        pi = [exp_mod(self.x2, self.a[i], self.curve.order)
              for i in range(len(self.a))]

        # copy array b of size N into matrix B of size mxn
        self.B = [pi[self.n * i: self.n * (i + 1)] for i in range(self.m)]

        # pedersen commitment of B
        self.c_B = self.pedersen.commit_matrix_vector(self.B, self.r_B)

        return self.c_B

    def r1_hadamard_zero(self):
        """Product argument: calculate D, F, G and commit to G

        Returns:
            List[ShortPoint], ShortPoint: commitment to G and z
        """
        # calculate D = y4 * A + B and F = D - z4
        self.D = []
        self.F = []
        for i in range(self.m):
            tmp0 = []
            tmp1 = []
            for j in range(self.n):
                var0 = mul_mod(self.A[i][j], self.y4, self.curve.order)
                var1 = add_mod(var0, self.B[i][j], self.curve.order)
                var2 = sub_mod(var1, self.z4, self.curve.order)
                tmp0.append(var1)
                tmp1.append(var2)
            self.D.append(tmp0)
            self.F.append(tmp1)

        # calculate G as Hadamard of F and create Pedersen commitment c_G
        # with random r_A of size m
        self.G = hadamard(self.F, self.m, self.n, self.curve.order)
        self.r_G = self.curve.rand_gen.get_random_array(self.m)
        self.c_G = self.pedersen.commit_matrix_vector(self.G, self.r_G)

        # fill vector of size n with z4 and commit to it
        self.z = [self.z4] * self.n
        self.r_z = 0
        self.c_z = self.pedersen.commit_vector_value(self.z, self.r_z)

        return self.c_G, self.c_z

    def r3_hadamard_zero(self):
        """Calculate the commitment c_F(0), c_H(m) and c_P. P(k) is the sum
        over the bilinear map of F(i) and H(j). F and H are modified:
        F: {F_0, F(2), F(3), ... , F(m), -1}
        H: {H(1), H(2), ... , H(m-1), H, H_m}

        Returns:
            ShortPoint, ShortPoint, List[ShortPoint]: commitment to F(0),
            H(m) and P
        """
        self.H = []
        self.r_H = []
        self.x6_array = [self.x6]
        # calculate vectors H(1), ... , H(m-1) as H(i) = G(i)*x^i
        # calculate values t_H(i) = r_G*x^i
        for i in range(self.m-1):
            tmp0 = []
            self.x6_array.append(mul_mod(self.x6_array[i], self.x6,
                                         self.curve.order))
            for j in range(self.n):
                var1 = mul_mod(self.x6_array[i], self.G[i][j],
                               self.curve.order)
                tmp0.append(var1)
            self.r_H.append(mul_mod(self.x6_array[i], self.r_G[i],
                                    self.curve.order))
            self.H.append(tmp0)

        tmp0 = [0] * self.n
        var0 = 0
        # calculate vector H = sum(x^i * G(i+1)) for i = 1, ... , m-1
        # calculate value t_H = sum(x^i * r_G(i+1)) for i = 1, ... , m-1
        for i in range(self.m - 1):
            for j in range(self.n):
                var1 = mul_mod(self.x6_array[i], self.G[i + 1][j],
                               self.curve.order)
                tmp0[j] = add_mod(tmp0[j], var1, self.curve.order)
            var2 = mul_mod(self.x6_array[i], self.r_G[i + 1], self.curve.order)
            var0 = add_mod(var0, var2, self.curve.order)
        self.H.append(tmp0)
        self.r_H.append(var0)

        # Pick H_m and t_H_m randomly and calculate commitment
        self.H_m = self.curve.rand_gen.get_random_array(self.n)
        self.r_H_m = self.curve.rand_gen.get_random_value()
        self.H.append(self.H_m)
        self.r_H.append(self.r_H_m)
        self.c_H_m = self.pedersen.commit_vector_value(self.H_m, self.r_H_m)

        # Pick F_0 and r_F_0 randomly and calculate commitment
        self.F_0 = self.curve.rand_gen.get_random_array(self.n)
        self.r_F_0 = self.curve.rand_gen.get_random_value()
        self.c_F_0 = self.pedersen.commit_vector_value(self.F_0, self.r_F_0)

        # Calculate -1
        var0 = neg_mod(1, self.curve.order)
        tmp0 = [var0] * self.n

        # Modify F: {F_0, F(2), F(3), ... , F(m), -1}
        self.F.pop(0)
        self.F.insert(0, self.F_0)
        self.F.insert(self.m, tmp0)

        # Calculate c_P as sum over bilinear map and the commitment for it
        var_l = 2 * self.m+1
        self.P = [0] * var_l
        self.r_P = self.curve.rand_gen.get_random_array(var_l)
        for k in range(var_l):
            var0 = 0
            for i in range(self.m+1):
                j = (self.m - k) + i
                if 0 <= j <= self.m:
                    var1 = bilinearmap(self.F[i], self.H[j], self.y6,
                                       self.curve.order)
                    var0 = add_mod(var0, var1, self.curve.order)
            self.P[k] = var0

        self.r_P[self.m + 1] = 0
        self.c_P = self.pedersen.commit_vector_vector(self.P, self.r_P)

        return self.c_F_0, self.c_H_m, self.c_P

    def r5_hadamard_zero(self):
        """Calculate second Zero Argument commitments

        Returns:
           List[int], int, List[int], int, int
        """
        x8_v = [self.x8]
        for i in range(2 * self.m):
            x8_v.append(mul_mod(x8_v[i], self.x8, self.curve.order))

        # f = sum(F(i) * x^i) for i = 0,...,m
        self.f = deepcopy(self.F[0][:])
        for i in range(self.m):
            for j in range(self.n):
                var0 = mul_mod(self.F[i + 1][j], x8_v[i], self.curve.order)
                var1 = add_mod(self.f[j], var0, self.curve.order)
                self.f[j] = var1

        # r_f = r_F_0 + sum((r_A(i)* y4 + r_B(i) - r_z)*x8^i)
        #         + r_-1 * x_8^m for i = 1, ... ,m-1
        self.r_f = self.r_F_0 + x8_v[self.m-1]
        for i in range(1, self.m):
            var0 = mul_mod(self.r_A[i], self.y4, self.curve.order)
            var1 = sub_mod(self.r_B[i], self.r_z, self.curve.order)
            var2 = add_mod(var0, var1, self.curve.order)
            var3 = mul_mod(var2, x8_v[i - 1], self.curve.order)
            self.r_f = add_mod(self.r_f, var3, self.curve.order)

        # h     = sum(x^(m-j)*H(j))
        # t_h   = sum(x^(m-j)*t_H(j)) for j = 0, ... ,m
        self.h = deepcopy(self.H[self.m])
        self.r_h = self.r_H[self.m]

        for i in range(self.m):
            for j in range(self.n):
                var1 = mul_mod(self.H[i][j], x8_v[self.m - i - 1],
                               self.curve.order)
                self.h[j] = add_mod(self.h[j], var1, self.curve.order)
            var2 = mul_mod(self.r_H[i], x8_v[self.m - i - 1], self.curve.order)
            self.r_h = add_mod(self.r_h, var2, self.curve.order)

        # t_p = sum(x^k * t_P(k)) for k = 0, ... ,2*m
        self.r_p = self.r_P[0]
        for i in range(self.m * 2):
            var0 = mul_mod(x8_v[i], self.r_P[i + 1], self.curve.order)
            self.r_p = add_mod(self.r_p, var0, self.curve.order)

        return self.f, self.r_f, self.h, self.r_h, self.r_p

    def r1_single_value_product(self):
        """Single value argument
        Returns:
            ShortPoint, ShortPoint, ShortPoint
        """
        self.g = self.G[self.m-1]
        self.alpha = [self.g[0]]
        for i in range(1, self.n):
            self.alpha.append(
                mul_mod(self.alpha[i - 1], self.g[i], self.curve.order))

        self.gamma = self.curve.rand_gen.get_random_array(self.n)
        self.r_gamma = self.curve.rand_gen.get_random_value()

        self.c_gamma = self.pedersen.commit_vector_value(
            self.gamma, self.r_gamma)

        self.delta = self.curve.rand_gen.get_random_array(self.n)
        self.delta[0] = self.gamma[0]
        self.delta[self.n - 1] = 0

        self.s_delta = self.curve.rand_gen.get_random_value()
        self.s_Delta = self.curve.rand_gen.get_random_value()

        var0 = []
        for i in range(self.n - 1):
            var1 = mul_mod(self.delta[i], self.gamma[i + 1], self.curve.order)
            var0.append(neg_mod(var1, self.curve.order))

        self.c_delta = self.pedersen.commit_vector_value(var0, self.s_delta)

        var0 = []
        for i in range(self.n - 1):
            var1 = mul_mod(self.g[i+1], self.delta[i], self.curve.order)
            var2 = mul_mod(self.alpha[i], self.gamma[i + 1], self.curve.order)
            var3 = sub_mod(self.delta[i+1], var1, self.curve.order)
            var0.append(sub_mod(var3, var2, self.curve.order))

        self.c_Delta = self.pedersen.commit_vector_value(var0, self.s_Delta)

        return self.c_gamma, self.c_delta, self.c_Delta

    def r3_single_value_product(self):
        """Single value argument
        Returns:
            List[int], List[int], int, int
        """
        self.gamma_tilde = []
        self.alpha_tilde = []
        for i in range(self.n):
            var0 = mul_mod(self.x6, self.g[i], self.curve.order)
            self.gamma_tilde.append(add_mod(
                var0, self.gamma[i], self.curve.order))
            var0 = mul_mod(self.x6, self.alpha[i], self.curve.order)
            self.alpha_tilde.append(add_mod(
                var0, self.delta[i], self.curve.order))

        var0 = mul_mod(self.x6, self.r_G[self.m - 1], self.curve.order)
        self.r_gamma_tilde = add_mod(var0, self.r_gamma, self.curve.order)

        var0 = mul_mod(self.x6, self.s_Delta, self.curve.order)
        self.r_alpha_tilde = add_mod(var0, self.s_delta, self.curve.order)

        return self.gamma_tilde, self.alpha_tilde, \
            self.r_gamma_tilde, self.r_alpha_tilde

    def r1_multi_exponent(self):
        """Calculate the diagonal sums E_k, B0, Beta and commit to them

        Returns:
            ShortPoint, List[ShortPoint], List[List[ShortPoint]]
        """
        # Commitment B_0
        self.b_0 = self.curve.rand_gen.get_random_array(self.n)
        self.r_B0 = self.curve.rand_gen.get_random_value()
        self.c_B0 = self.pedersen.commit_vector_value(self.b_0, self.r_B0)

        # Commitment Beta
        self.beta = self.curve.rand_gen.get_random_array(2 * self.m)
        self.beta[self.m] = 0
        self.r_beta = self.curve.rand_gen.get_random_array(2 * self.m)
        self.r_beta[self.m] = 0
        self.c_beta = self.pedersen.commit_vector_vector(
            self.beta, self.r_beta)

        # Commitment E_k
        self.tau_k = self.curve.rand_gen.get_random_array(2 * self.m)

        # Calculate tau(m) = -sum(rho(i)*b(i)) for i = 1, ... , N
        var0 = 0
        for i in range(self.N):
            var1 = mul_mod(self.rho[i], self.B[int(i / self.n)][i % self.n],
                           self.curve.order)
            var0 = add_mod(var0, var1, self.curve.order)
        self.tau_k[self.m] = neg_mod(var0, self.curve.order)

        # Form reencrypted and shuffled cards into matrix of size mxn
        c_shuffled = [self.ciphers_out[self.n * i: self.n * (i + 1)]
                      for i in range(self.m)]

        # Calculate Ek = Enc(G^c_k; tau_k) + sum(C'_i*b_j)
        # for i = 1,...,m, j = 0,...,m, j = (k-m)+i
        b_v = deepcopy(self.B)
        b_v.insert(0, self.b_0)

        # Calculate diagonal sums
        c_prod = []
        for k in range(2 * self.m):
            var0 = None
            for i in range(self.m):
                j = k - self.m + i + 1
                if 0 <= j <= self.m:
                    var1 = mul_cipher(b_v[j][0], c_shuffled[i][0], self.curve)
                    for i2 in range(1, self.n):
                        var2 = mul_cipher(b_v[j][i2], c_shuffled[i][i2],
                                          self.curve)
                        var1 = add_cipher(var2, var1, self.curve)
                    if var0 is None:
                        var0 = var1
                    else:
                        var0 = add_cipher(var0, var1, self.curve)
            c_prod.append(var0)

        # Add encryption E_pk(G*c_k; tau_k)
        kha = [self.curve.multiplication(self.tau_k[x], self.pubKey) for x
               in range(2*self.m)]
        enc = []
        for i in range(2 * self.m):
            enc_a = self.curve.multiplication(self.tau_k[i],
                                              self.curve.generator)
            var0 = self.curve.multiplication(
                self.beta[i], self.curve.generator)
            enc_b = self.curve.addition(var0, kha[i])
            enc.append([enc_a, enc_b])

        self.E = [add_cipher(enc[i], c_prod[i], self.curve)
                  for i in range(2 * self.m)]

        # calculate E(m) separately with c(m) = 0
        enc_a = self.curve.multiplication(self.tau_k[self.m],
                                          self.curve.generator)
        enc_b = self.curve.multiplication(self.tau_k[self.m], self.pubKey)
        enc = [enc_a, enc_b]

        self.E[self.m] = add_cipher(enc, c_prod[self.m], self.curve)

        return self.c_B0, self.c_beta, self.E

    def r3_multi_exponent(self):
        """Calculate commitment as sums over random values from last round
        using challenge x6_array = (x, x^2, ... , x^m)^T

        Returns:
            List[int], int, int, int, int
        """
        x6_array = [exp_mod(self.x6, i + 1, self.curve.order)
                    for i in range(self.m)]

        # Commitment vector b and value r_b
        self.b = []
        for j in range(self.n):
            var0 = self.b_0[j]
            for i in range(self.m):
                var1 = mul_mod(self.B[i][j], x6_array[i], self.curve.order)
                var0 = add_mod(var0, var1, self.curve.order)
            self.b.append(var0)

        self.r_b = deepcopy(self.r_B0)
        for i in range(self.m):
            var0 = mul_mod(self.r_B[i], x6_array[i], self.curve.order)
            self.r_b = add_mod(self.r_b, var0, self.curve.order)

        # Commitment beta_tilde and r_beta_tilde
        self.beta_tilde = self.beta[0]
        self.r_beta_tilde = self.r_beta[0]
        for i in range(2*self.m - 1):
            var0 = exp_mod(self.x6, i + 1, self.curve.order)
            var1 = mul_mod(var0, self.beta[i + 1], self.curve.order)
            var2 = mul_mod(var0, self.r_beta[i + 1], self.curve.order)
            self.beta_tilde = add_mod(self.beta_tilde, var1, self.curve.order)
            self.r_beta_tilde = add_mod(
                self.r_beta_tilde, var2, self.curve.order)

        # Commitment tau to check E
        self.tau = self.tau_k[0]
        for i in range(1, 2*self.m):
            var0 = exp_mod(self.x6, i, self.curve.order)
            var1 = mul_mod(self.tau_k[i], var0, self.curve.order)
            self.tau = add_mod(self.tau, var1, self.curve.order)

        return self.b, self.r_b, self.beta_tilde, self.r_beta_tilde, self.tau

    def round1(self):
        """Round 1

        Returns:
            List[ShortPoint]: commitment to A
        """
        return self.r1_shuffle()

    def round3(self, x2):
        """Round 3

        Args:
            x2 (int): challenge from verifier

        Returns:
            List[ShortPoint]: commitment to A
        """
        return self.r3_shuffle(x2)

    def round5(self, y4, z4):
        """Round 5

        Args:
            y4 (int): challenge from verifier
            z4 (int): challenge from verifier

        Returns:
            return values from the functions as stated bellow
        """
        self.y4 = y4
        self.z4 = z4

        var0, var1 = self.r1_hadamard_zero()
        var2, var3, var4 = self.r1_multi_exponent()
        var5, var6, var7 = self.r1_single_value_product()

        return var0, var1, var2, var3, var4, var5, var6, var7

    def round7(self, x6, y6):
        """Round 7

        Args:
            x6 (int): challenge from verifier
            y6 (int): challenge from verifier

        Returns:
            return values from the functions as stated bellow
        """
        self.x6 = x6
        self.y6 = y6
        var0, var1, var2 = self.r3_hadamard_zero()
        var3, var4, var5, var6, var7 = self.r3_multi_exponent()
        var8, var9, var10, var11 = self.r3_single_value_product()
        return var0, var1, var2, var3, var4, var5, var6, var7, var8, var9, \
            var10, var11

    def round9(self, x8):
        """Round 9

        Args:
            x8 (int): challenge from verifier

        Returns:
            return values from the functions as stated bellow
        """
        self.x8 = x8
        return self.r5_hadamard_zero()

    def nizk_prover(self):
        """Proof non-interactive Zero-Knowledge Argument for Correctness of a
        Shuffle, challenge is generated by hash

        Returns:
            Proofs from all zero-knowledge arguments
        """
        # R1
        self.r1_shuffle()

        # R2
        # G, ck, pk, C, C' c_A, N, m, n, ID
        # seed = H(G, ck, c_A)
        # x2 = PRG(1, seed)
        list_ck = curve_points_to_list(self.generators_ck)
        var1 = curve_points_to_list(self.c_A)

        x2, seed = \
            self.curve.rand_gen.get_random_from_hash(1, self.curve.order, 1,
                                                     self.curve.generator.x,
                                                     self.curve.generator.y,
                                                     *list_ck, *var1)

        # R3
        self.r3_shuffle(x2[0])

        # R4
        # seed = H(G, ck, c_B, x2)
        # y4 = PRG(1, seed)
        # z4 = PRG(2, seed)
        var1 = curve_points_to_list(self.c_B)

        values, seed = \
            self.curve.rand_gen.get_random_from_hash(1, self.curve.order, 2,
                                                     self.curve.generator.x,
                                                     self.curve.generator.y,
                                                     *list_ck, *var1, self.x2)
        self.y4 = values[0]
        self.z4 = values[1]

        # R5
        self.r1_hadamard_zero()
        self.r1_multi_exponent()
        self.r1_single_value_product()

        # R6
        # seed = H(G, ck, c_G, c_z, c_B0, c_C, E,
        #          c_u, c_delta, c_Delta, y4, z4)
        # x6 = PRG(1, seed)
        # y6 = PRG(2, seed)
        var1 = curve_points_to_list(self.c_G)
        var2 = curve_points_to_list(self.c_beta)
        var3 = ciphers_to_list(self.E)

        values, seed = self.curve.rand_gen.get_random_from_hash(
            1, self.curve.order, 2, self.curve.generator.x,
            self.curve.generator.y, *list_ck, *var1, self.c_z.x, self.c_z.y,
            self.c_B0.x, self.c_B0.y, *var2, *var3, self.c_gamma.x,
            self.c_gamma.y, self.c_delta.x, self.c_delta.y, self.c_Delta.x,
            self.c_Delta.y, self.y4, self.z4)
        self.x6 = values[0]
        self.y6 = values[1]

        # R7
        self.r3_hadamard_zero()
        self.r3_multi_exponent()
        self.r3_single_value_product()

        # R8
        # seed = H(G, ck, c_F_0, c_H_m, c_P, b, s, c, t, tau,
        # u_tilde, q_tilde, r_u_tilde, r_q_tilde, x6, y6)
        # x8 = PRG(1, seed)
        var1 = curve_points_to_list(self.c_P)
        values, seed = self.curve.rand_gen.get_random_from_hash(
            1, self.curve.order, 1, self.curve.generator.x,
            self.curve.generator.y, *list_ck, self.c_F_0.x, self.c_F_0.y,
            self.c_H_m.x, self.c_H_m.y, *var1, *self.b, self.r_b,
            self.beta_tilde, self.r_beta_tilde, self.tau, *self.gamma_tilde,
            *self.alpha_tilde, self.r_gamma_tilde, self.r_alpha_tilde,
            self.x6, self.y6)

        self.x8 = values[0]

        # R9
        self.r5_hadamard_zero()
        return [[self.c_A], [self.c_B], [self.c_G, self.c_z, self.c_B0,
                                         self.c_beta, self.E, self.c_gamma,
                                         self.c_delta, self.c_Delta],
                [self.c_F_0, self.c_H_m, self.c_P, self.b, self.r_b,
                 self.beta_tilde, self.r_beta_tilde, self.tau,
                 self.gamma_tilde, self.alpha_tilde, self.r_gamma_tilde,
                 self.r_alpha_tilde], [self.f, self.r_f, self.h, self.r_h,
                                       self.r_p]]


class BayGroVerifier:
    """Verifier in Zero-Knowledge Argument for Correctness of a Shuffle such
    that ciphers_out[i] = ciphers_in[pi[i]] + Enc_pk(O, rho[i])
    """
    def __init__(self, N, m, n, curve, gen, pk, ciphers_out):
        """
        Args:
            N (int): number of cards N=m*n
            m (int): rows
            n (int): columns
            curve (ECCObj): curve object
            gen (List[ShortPoint]): generators for pedersen commitment
            pk (ShortPoint): public key
            ciphers_out (List[List[ShortPoint]]): ciphers after the
                permutation and re-masking
        """
        self.m = m
        self.n = n
        self.N = N
        assert (m * n == N)

        self.curve = curve

        self.generators_ck = gen
        self.pedersen = Pedersen(gen, n, self.curve)
        self.order = self.curve.order
        self.pubKey = deepcopy(pk)
        self.ciphers_out = ciphers_out

        # R2
        self.x2 = None
        self.c_A = None

        # R4
        self.y4, self.z4 = (None,)*2
        self.c_B = None

        # R2 (Hadamard/Zero)
        self.x6, self.y6, self.x6_array = (None,)*3
        self.c_G, self.c_z = (None,)*2

        # R4 (Hadamard/Zero)
        self.x8 = None
        self.c_F_0 = None
        self.c_H_m = None
        self.c_P = None

        # R6/Verify (Hadamard/Zero)
        self.f, self.r_f = (None,)*2
        self.h, self.r_h = (None,) * 2
        self.r_p = None
        self.c_F = None
        self.c_H = None

        # R2 (Single Value Product)
        self.c_gamma, self.c_delta, self.c_Delta = (None,) * 3

        # R4 (Single Value Product)
        self.gamma_tilde, self.r_gamma_tilde = (None,) * 2
        self.alpha_tilde, self.r_alpha_tilde = (None,) * 2

        # R2 (Multi-Exponent)
        self.c_B0, self.c_beta, self.E = (None,) * 3

        # R4/Verify (Multi-Exponent)
        self.b, self.r_b = (None,) * 2
        self.beta_tilde, self.r_beta_tilde = (None,) * 2
        self.tau = None

    def r2_shuffle(self, c_A):
        """Save commitment c_A and send random challenge x2

        Args:
            c_A (List[ShortPoint]): commitment to A

        Returns:
            int: random challenge x2
        """
        self.c_A = c_A
        self.x2 = self.curve.rand_gen.get_random_value()

        return self.x2

    def r4_shuffle(self, c_B):
        """Save commitment c_B and send random challenge y4 and z4

        Args:
            c_B (List[ShortPoint]): commitment to B

        Returns:
            int, int: random challenge y4 and z4
        """
        self.c_B = c_B
        self.y4 = self.curve.rand_gen.get_random_value()
        self.z4 = self.curve.rand_gen.get_random_value()

        return self.y4, self.z4

    def r2_hadamard_zero(self, c_G, c_z):
        """Save commitment c_G and c_z

        Args:
            c_G (List[ShortPoint]): commitment to G
            c_z (ShortPoint): commitment to z
        """
        self.c_G = c_G
        self.c_z = c_z

    def r4_hadamard_zero(self, c_F_0, c_H_m, c_P):
        """Save commitment c_F_0, c_H_m, c_P

        Args:
            c_F_0 (ShortPoint): commitment to F[0]
            c_H_m (ShortPoint): commitment to H[m]
            c_P (List[ShortPoint]): commitment to P
        """
        self.c_F_0 = c_F_0
        self.c_H_m = c_H_m
        self.c_P = c_P

    def r6_verify_hadamard_zero(self, f, r_f, h, t_h, t_p):
        """Save Zero commitments f, r_f, h, t_h, t_p.
        Verify sum(c_F*x^i) = com_ck(f;r_f) for i  = 0, ... ,m
        Verify sum(c_H*x^(m-j)) = com_ck(h;t_h) for j = 0, ...,m
        Verify sum(c_P*x^k) = com_ck(bilinearmap(f,h);t_p) for k = 0, ... ,2m

        Args:
            f (List[int]): commitment from prover
            r_f (int): commitment from prover
            h (List[int]): commitment from prover
            t_h (int): commitment from prover
            t_p (int): commitment from prover

        Returns:
            bool: True if verification successful, False else
        """
        self.f = f
        self.r_f = r_f
        self.h = h
        self.r_h = t_h
        self.r_p = t_p

        verification = True

        # Calculate x8_array = {x8, x8^2, ... , x8^2m}
        x8_array = [self.x8]
        for i in range(2*self.m-1):
            x8_array.append(mul_mod(x8_array[i], self.x8, self.curve.order))

        # Verify sum(c_F*x^i) = com_ck(f;r_f) for i  = 0, ... ,m
        self.c_F = [self.c_F_0]
        for i in range(self.m-1):
            var0 = self.curve.multiplication(self.y4, self.c_A[i + 1])
            var1 = self.curve.subtraction(self.c_B[i + 1], self.c_z)
            self.c_F.append(self.curve.addition(var0, var1))

        var0 = neg_mod(1, self.curve.order)
        var1 = [var0] * self.n

        var0 = self.pedersen.commit_vector_value(var1, 1)

        self.c_F.append(var0)

        ver_c_f = self.c_F[0]
        for i in range(self.m):
            var0 = self.curve.multiplication(x8_array[i], self.c_F[i + 1])
            ver_c_f = self.curve.addition(ver_c_f, var0)

        ped_c_f = self.pedersen.commit_vector_value(self.f, self.r_f)

        if ver_c_f != ped_c_f:
            verification = False

        # Verify sum(c_H*x^(m-j)) = com_ck(h;t_h) for j = 0, ...,m
        # c_Hi = c_Gi*x6^i
        self.c_H = [self.curve.multiplication(self.x6_array[i + 1],
                    self.c_G[i]) for i in range(self.m-1)]

        # c_H = sum(c_G(i+1)*x^i
        if self.m > 1:
            var0 = self.curve.multiplication(self.x6_array[1], self.c_G[1])
            for i in range(1, self.m-1):
                var1 = self.curve.multiplication(self.x6_array[i + 1],
                                                 self.c_G[i+1])
                var0 = self.curve.addition(var0, var1)

            self.c_H.append(var0)
        self.c_H.append(self.c_H_m)

        if self.m > 1:
            ver_c_h = self.c_H[self.m]
            for i in range(self.m):
                var0 = self.curve.multiplication(x8_array[self.m - i - 1],
                                                 self.c_H[i])
                ver_c_h = self.curve.addition(ver_c_h, var0)
        else:
            ver_c_h = self.c_H[0]

        ped_c_h = self.pedersen.commit_vector_value(self.h, self.r_h)

        if ver_c_h != ped_c_h:
            verification = False

        # Verify sum(c_P*x^k) = com_ck(bilinearmap(f,h);t_p) for k = 0, ... ,2m
        ver_c_p = self.c_P[0]
        for i in range(2*self.m):
            var1 = self.curve.multiplication(x8_array[i], self.c_P[i + 1])
            ver_c_p = self.curve.addition(ver_c_p, var1)

        var0 = bilinearmap(self.f, self.h, self.y6, self.curve.order)
        ped_c_p = self.pedersen.commit_vector_vector([var0], [self.r_p])

        if ver_c_p != ped_c_p[0]:
            verification = False

        return verification

    def r2_single_value_product(self, c_gamma, c_delta, c_Delta):
        """Save commitment c_gamma, c_delta, c_Delta
        Args:
            c_gamma (ShortPoint): commitment from prover
            c_delta (ShortPoint): commitment from prover
            c_Delta (ShortPoint): commitment from prover
        """
        self.c_gamma = c_gamma
        self.c_delta = c_delta
        self.c_Delta = c_Delta

    def r4_verify_single_value_product(self, gamma_tilde, alpha_tilde,
                                       r_gamma_tilde, r_alpha_tilde):
        """Save commitment gamma_tilde, alpha_tilde, r_gamma_tilde,
        r_alpha_tilde and verify single value argument.

        Args:
            gamma_tilde (List[int]): commitment from prover
            alpha_tilde (List[int]): commitment from prover
            r_gamma_tilde (int): commitment from prover
            r_alpha_tilde (int): commitment from prover

        Returns:
            bool: True if verification successful, False else
        """
        self.gamma_tilde = gamma_tilde
        self.alpha_tilde = alpha_tilde
        self.r_gamma_tilde = r_gamma_tilde
        self.r_alpha_tilde = r_alpha_tilde

        verification = True

        if gamma_tilde[0] != alpha_tilde[0]:
            verification = False

        var_x = 1
        g = 1
        for i in range(1, self.N+1):
            var_x = mul_mod(var_x, self.x2, self.curve.order)
            var_y = mul_mod(i, self.y4, self.curve.order)
            var0 = add_mod(var_x, var_y, self.curve.order)
            var1 = sub_mod(var0, self.z4, self.curve.order)
            g = mul_mod(var1, g, self.curve.order)

        g = mul_mod(g, self.x6, self.curve.order)

        if alpha_tilde[self.n - 1] != g:
            verification = False

        ped = self.pedersen.commit_vector_value(
            self.gamma_tilde, self.r_gamma_tilde)

        var0 = self.curve.multiplication(self.x6, self.c_G[self.m-1])
        ver = self.curve.addition(var0, self.c_gamma)

        if ped != ver:
            verification = False

        var0 = []
        for i in range(self.n-1):
            var1 = mul_mod(self.x6, self.alpha_tilde[i + 1], self.curve.order)
            var2 = mul_mod(self.alpha_tilde[i], self.gamma_tilde[i + 1],
                           self.curve.order)
            var0.append(sub_mod(var1, var2, self.curve.order))

        ped = self.pedersen.commit_vector_value(var0, self.r_alpha_tilde)

        var0 = self.curve.multiplication(self.x6, self.c_Delta)
        ver = self.curve.addition(var0, self.c_delta)

        if ped != ver:
            verification = False

        return verification

    def r2_multi_exponent(self, c_B0, c_C, E):
        """Save Multi-Exponent commitments c_B0, c_C, E

        Args:
            c_B0 (ShortPoint): commitment from prover
            c_C (List[ShortPoint]): commitment from prover
            E (List[List[ShortPoint]]): commitment from prover
        """
        self.c_B0 = c_B0
        self.c_beta = c_C
        self.E = E

    def r4_verify_multi_exponent(self, b, r_b, beta_tilde, r_beta_tilde, tau):
        """Save Multi-Exponent commitments b, s, c, t, tau.
        Verify c_B0+c_B*x6_array = com_ck(b,s).
        For k = 1, ... ,2m-1:
        Verify c_C0+sum(c_Ck*x^k) = com_ck(c,t)
        Verify E0+sum(E_k*x^k) = Enc_pk(G^c, tau)+sum(C_i*(x^m-i*b))

        Args:
            b (List[int]): commitment from prover
            r_b (int): commitment from prover
            beta_tilde (int): commitment from prover
            r_beta_tilde (int): commitment from prover
            tau (int): commitment from prover

        Returns:
            bool: True if verification successful, False else
        """
        self.b = b
        self.r_b = r_b
        self.beta_tilde = beta_tilde
        self.r_beta_tilde = r_beta_tilde
        self.tau = tau

        verification = True

        # Verify c_B0+c_B*x6_array = com_ck(b,s).
        ped = self.pedersen.commit_vector_value(self.b, self.r_b)

        ver = self.c_B0
        for i in range(self.m):
            var0 = self.curve.multiplication(self.x6_array[i + 1], self.c_B[i])
            ver = self.curve.addition(ver, var0)
        if ped != ver:
            verification = False

        # Verify c_C0+sum(c_Ck*x^k) = com_ck(c,t)
        ped = self.pedersen.commit_vector_vector(
            [self.beta_tilde], [self.r_beta_tilde])

        ver = self.c_beta[0]
        for i in range(2*self.m - 1):
            var1 = self.curve.multiplication(self.x6_array[i + 1],
                                             self.c_beta[i + 1])
            ver = self.curve.addition(ver, var1)

        if ped[0] != ver:
            verification = False

        # Verify E0+sum(E_k*x^k) = Enc_pk(G^c, tau)+sum(C_i*(x^m-i*b))
        c_shuffled = [self.ciphers_out[self.n*i: self.n*(i+1)]
                      for i in range(self.m)]

        # E0+sum(E_k*x^k)
        e1 = self.E[0]
        for i in range(1, 2 * self.m):
            var1 = mul_cipher(self.x6_array[i], self.E[i], self.curve)
            e1 = add_cipher(var1, e1, self.curve)

        # Enc_pk(G^c, tau)
        enc_a = self.curve.multiplication(self.tau, self.curve.generator)
        kha = self.curve.multiplication(self.tau, self.pubKey)
        var0 = self.curve.multiplication(self.beta_tilde, self.curve.generator)
        enc_b = self.curve.addition(var0, kha)

        e2_1 = [enc_a, enc_b]

        # sum(C_i * (x ^ m - i * b))
        e2_2 = []
        for i in range(self.m):
            tmp = None
            for j in range(self.n):
                var1 = mul_mod(self.b[j], self.x6_array[self.m - i - 1],
                               self.curve.order)
                var2 = mul_cipher(var1, c_shuffled[i][j], self.curve)

                if tmp is None:
                    tmp = var2
                else:
                    tmp = add_cipher(tmp, var2, self.curve)
            e2_2.append(tmp)

        e2 = e2_2[0]
        for i in range(1, self.m):
            e2 = add_cipher(e2, e2_2[i], self.curve)

        e2 = add_cipher(e2, e2_1, self.curve)

        if e1 != e2:
            verification = False

        return verification

    def round2(self, c_A):
        return self.r2_shuffle(c_A)

    def round4(self, c_B):
        return self.r4_shuffle(c_B)

    def round6(self, c_G, c_z, c_B0, c_C, E, c_u, c_delta, c_Delta):
        self.x6 = self.curve.rand_gen.get_random_value()
        self.y6 = self.curve.rand_gen.get_random_value()

        self.x6_array = [exp_mod(self.x6, i, self.curve.order)
                         for i in range(2*self.m + 1)]

        self.r2_hadamard_zero(c_G, c_z)
        self.r2_multi_exponent(c_B0, c_C, E)
        self.r2_single_value_product(c_u, c_delta, c_Delta)
        return self.x6, self.y6

    def round8(self, c_F_0, c_H_m, c_P, b, s, c, t, tau, u_tilde, q_tilde,
               r_u_tilde, r_q_tilde):
        self.x8 = self.curve.rand_gen.get_random_value()

        self.r4_hadamard_zero(c_F_0, c_H_m, c_P)
        verification1 = self.r4_verify_multi_exponent(b, s, c, t, tau)
        verification2 = self.r4_verify_single_value_product(u_tilde, q_tilde,
                                                            r_u_tilde,
                                                            r_q_tilde)
        return self.x8, verification1, verification2

    def round10(self, f, r_f, h, t_h, t_p):
        return self.r6_verify_hadamard_zero(f, r_f, h, t_h, t_p)

    def nizk_verifier(self, nizk_proof):
        """Verification non-interactive Zero-Knowledge Argument for
        Correctness of a Shuffle, challenge is generated by hash

        Args:
            nizk_proof: nizk proof

        Returns:
            bool, bool, bool: True if verification successful, False else
            for all three arguments: Multi-Exponent, Single Value, Hadamard

        """
        self.c_A = nizk_proof[0][0]

        self.c_B = nizk_proof[1][0]

        self.c_G = nizk_proof[2][0]
        self.c_z = nizk_proof[2][1]
        self.c_B0 = nizk_proof[2][2]
        self.c_beta = nizk_proof[2][3]
        self.E = nizk_proof[2][4]
        self.c_gamma = nizk_proof[2][5]
        self.c_delta = nizk_proof[2][6]
        self.c_Delta = nizk_proof[2][7]

        self.c_F_0 = nizk_proof[3][0]
        self.c_H_m = nizk_proof[3][1]
        self.c_P = nizk_proof[3][2]
        self.b = nizk_proof[3][3]
        self.r_b = nizk_proof[3][4]
        self.beta_tilde = nizk_proof[3][5]
        self.r_beta_tilde = nizk_proof[3][6]
        self.tau = nizk_proof[3][7]
        self.gamma_tilde = nizk_proof[3][8]
        self.alpha_tilde = nizk_proof[3][9]
        self.r_gamma_tilde = nizk_proof[3][10]
        self.r_alpha_tilde = nizk_proof[3][11]

        self.f = nizk_proof[4][0]
        self.r_f = nizk_proof[4][1]
        self.h = nizk_proof[4][2]
        self.r_h = nizk_proof[4][3]
        self.r_p = nizk_proof[4][4]

        # R2
        # G, ck, pk, C, C' c_A, N, m, n, ID
        # seed = H(G, ck, c_A)
        # x2 = PRG(1, seed)
        list_ck = curve_points_to_list(self.generators_ck)
        var1 = curve_points_to_list(self.c_A)

        x2, seed = \
            self.curve.rand_gen.get_random_from_hash(1, self.curve.order, 1,
                                                     self.curve.generator.x,
                                                     self.curve.generator.y,
                                                     *list_ck, *var1)
        self.x2 = x2[0]
        # R4
        # seed = H(G, ck, c_B, x2)
        # y4 = PRG(1, seed)
        # z4 = PRG(2, seed)
        var1 = curve_points_to_list(self.c_B)

        values, seed = \
            self.curve.rand_gen.get_random_from_hash(1, self.curve.order, 2,
                                                     self.curve.generator.x,
                                                     self.curve.generator.y,
                                                     *list_ck, *var1, self.x2)
        self.y4 = values[0]
        self.z4 = values[1]

        # R6
        # seed = H(G, ck, c_G, c_z, c_B0, c_C, E),
        #          c_u, c_delta, c_Delta, y4, z4)
        # x6 = PRG(1, seed)
        # y6 = PRG(2, seed)
        var1 = curve_points_to_list(self.c_G)
        var2 = curve_points_to_list(self.c_beta)
        var3 = ciphers_to_list(self.E)

        values, seed = self.curve.rand_gen.get_random_from_hash(
            1, self.curve.order, 2, self.curve.generator.x,
            self.curve.generator.y, *list_ck, *var1, self.c_z.x, self.c_z.y,
            self.c_B0.x, self.c_B0.y, *var2, *var3, self.c_gamma.x,
            self.c_gamma.y, self.c_delta.x, self.c_delta.y, self.c_Delta.x,
            self.c_Delta.y, self.y4, self.z4)
        self.x6 = values[0]
        self.y6 = values[1]

        self.x6_array = [exp_mod(self.x6, i, self.curve.order)
                         for i in range(2*self.m + 1)]

        # R8
        # seed = H(G, ck, c_F_0, c_H_m, c_P, b, s, c, t, tau,
        # u_tilde, q_tilde, r_u_tilde, r_q_tilde, x6, y6)
        # x8 = PRG(1, seed)
        var1 = curve_points_to_list(self.c_P)
        values, seed = self.curve.rand_gen.get_random_from_hash(
            1, self.curve.order, 1, self.curve.generator.x,
            self.curve.generator.y, *list_ck, self.c_F_0.x, self.c_F_0.y,
            self.c_H_m.x, self.c_H_m.y, *var1, *self.b, self.r_b,
            self.beta_tilde, self.r_beta_tilde, self.tau, *self.gamma_tilde,
            *self.alpha_tilde, self.r_gamma_tilde, self.r_alpha_tilde,
            self.x6, self.y6)

        self.x8 = values[0]

        verification1 = self.r4_verify_multi_exponent(
            self.b, self.r_b, self.beta_tilde, self.r_beta_tilde, self.tau)
        verification2 = self.r4_verify_single_value_product(self.gamma_tilde,
                                                            self.alpha_tilde,
                                                            self.r_gamma_tilde,
                                                            self.r_alpha_tilde)

        verification3 = self.r6_verify_hadamard_zero(self.f, self.r_f, self.h,
                                                     self.r_h, self.r_p)

        return verification1, verification2, verification3


class Pedersen:
    """Pedersen commitment: com_ck(a_1,...,a_n;r) = g_1*a_1+...+g_n*a_n+h*r
    with randomness r

    Attributes:
        Curve (ECCobj): elliptic curve
        gen_G_Curve (List[ShortPoint]): generators g_1,...g_n
        gen_H_Curve (ShortPoint): generator h
    """
    def __init__(self, gen, n, Curve):
        """
        Args:
            gen (List[ShortPoint]): generators for pedersen commitment
            n (int): columns of shuffle argument
            Curve (ECCobj): elliptic curve
        """
        self.Curve = Curve

        self.gen_G_Curve = gen[0:n]
        self.gen_H_Curve = gen[n]

    def commit_vector_value(self, a_v, r):
        """Commit to n values in a_v with randomness r

        Args:
            a_v (List[int]): elements to for commitment
            r (int): randomness

        Returns:
            ShortPoint: commitment

        """
        var0: ShortPoint = self.Curve.multiplication(a_v[0],
                                                     self.gen_G_Curve[0])
        for i in range(1, len(a_v)):
            var1 = self.Curve.multiplication(a_v[i], self.gen_G_Curve[i])
            var0 = self.Curve.addition(var0, var1)

        if r != 0:
            var2 = self.Curve.multiplication(r, self.gen_H_Curve)
            var0 = self.Curve.addition(var0, var2)

        return var0

    def commit_matrix_vector(self, A_v, r_v):
        """Commit to m*n values with randomness r_v

        Args:
            A_v (List[List[int]]): m vectors with element for commitment
            r_v (List[int]): m randomness values

        Returns:
            List[ShortPoint]: m commitments
        """
        assert(len(A_v) <= len(r_v))
        var0 = []
        for i in range(len(A_v)):
            var0.append(self.commit_vector_value(A_v[i], r_v[i]))

        return var0

    def commit_vector_vector(self, p, t_p):
        """Commitment: com_ck(p[i],0,...,0; t_p[i]
        Args:
            p (List[int]): values for commitment
            t_p (List[int]): randomness for commitment

        Returns:
            List[ShortPoint]: len(p) == len(t_p) commitments
        """
        assert(len(p) == len(t_p))
        var0 = []
        for i in range(len(p)):
            var2 = self.Curve.multiplication(p[i], self.gen_G_Curve[0])
            var1 = self.Curve.multiplication(t_p[i], self.gen_H_Curve)
            var2 = self.Curve.addition(var1, var2)

            var0.append(var2)
        return var0


def mul_cipher(a, b, curve):
    """Multiply cipher with integer
    Args:
        a (int): integer
        b ([ShortPoint, ShortPoint]): cipher
        curve (ECCobj): elliptic curve

    Returns:
        [ShortPoint, ShortPoint]: multiplied cipher
    """
    var0 = curve.multiplication(a, b[0])
    var1 = curve.multiplication(a, b[1])

    return [var0, var1]


def add_cipher(a, b, curve: ECCobj):
    """Add two ciphers
    Args:
        a ([ShortPoint, ShortPoint]): cipher
        b ([ShortPoint, ShortPoint]): cipher
        curve (ECCobj): elliptic curve

    Returns:
        [ShortPoint, ShortPoint]: added cipher
    """
    var0 = curve.addition(b[0], a[0])
    var1 = curve.addition(b[1], a[1])

    return [var0, var1]


def add_mod(a: int, b: int, order: int) -> int:
    var = (a+b) % order
    return var


def mul_mod(a: int, b: int, order: int) -> int:
    var = (a * b) % order
    return var


def exp_mod(a: int, b: int, order: int) -> int:
    var = (a**b) % order
    return var


def sub_mod(a: int, b: int, order: int) -> int:
    var = (a-b) % order
    return var


def neg_mod(a: int, order: int) -> int:
    var = (-a) % order
    return var


def hadamard(A, m, n, order):
    """Calculate Hadamard product: b_0 = A[0], b_1 = A[0]*A[1],...,
    b_m = A[0]*A[1]*...*A[m]

    Args:
        A (List[List[int]]): Matrix with n*m elements
        m (int): columns
        n (int): rows
        order (int): order of elliptic curve subgroup

    Returns:
        List[List[int]]: Hadamard product
    """
    var0 = [A[0]]
    for i in range(1, m):
        tmp = []
        for j in range(n):
            var1 = mul_mod(A[i][j], var0[i - 1][j], order)
            tmp.append(var1)
        var0.append(tmp)

    return var0


def bilinearmap(f, h, y, order):
    """Bilinear map: a = sum_{j=1}^{n}(f_j*h_j*y^j)

    Args:
        f (List[int]): list 1
        h (List[int]): list 2
        y (int): challenge
        order (int): order of elliptic curve subgroup

    Returns:
        int: solution from bilinear map
    """
    n = len(f)
    assert (n == len(h))

    var0 = 0

    for j in range(n):
        var1 = mul_mod(f[j], h[j], order)
        var2 = exp_mod(y, j + 1, order)
        var3 = mul_mod(var1, var2, order)

        var0 = add_mod(var0, var3, order)

    return var0


def ciphers_to_list(ciphers):
    """Fill list with elements from ciphers

    Args:
        ciphers (List[[ShortPoint, ShortPoint]]): ciphers

    Returns:
        List[int]
    """
    var = []
    for var0 in ciphers:
        for var1 in var0:
            var.append(var1.x)
            var.append(var1.y)
    return var


def curve_points_to_list(curve_points):
    """Fill list with elements from ShortPoints

    Args:
        curve_points ((List[ShortPoint])): point on elliptic curve

    Returns:
        List(int): list with x and y values from ShortPoint
    """
    var = []
    for var0 in curve_points:
        var.append(var0.x)
        var.append(var0.y)
    return var
