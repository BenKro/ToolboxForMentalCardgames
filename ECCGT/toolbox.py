# ! /usr/bin/env python
# -*- coding: utf-8 -*-
# ===================================================================
# toolbox.py
#
# 10.09.2019, Benjamin Kromer <benjamin.kromer1@web.de>
#
# @desc: Toolbox for mental card games using the elliptic curve encryption
#        scheme. Functions for masking, re-masking and unmasking cards are
#        implemented, such as for key generation, player order generation,
#        shuffle cards and open/draw cards.
# ===================================================================
from ECCGT.eccwrapper import ShortPoint
import ECCGT.proofs as proofs
from ECCGT.bayergroth import BayGroProver, BayGroVerifier


class Toolbox:
    """Toolbox for Mental Card Games

    Attributes:
        N (int): number of cards
        curve: elliptic curve
        secret_key (int): own secret key
        public_key_share_own (ShortPoint): own public key share
        keygen_proof (List[int]): proof that public_key_share_own =
            generator*secret_key
        public_key (ShortPoint): common public key
        keygen_successful (List[bool]): keygen successful flags
        n (int): number of columns for shuffle proof
        m (int): number of rows for shuffle proof
        bay_gro_secret_keys (List[int]): secret keys for commitment generators
        bay_gro_generator_shares_proofs (List[List[int]]): proofs for
            generator shares
        bay_gro_generator_shares_data (List[List[ShortPoint]]): generator
            shares
        bay_gro_generators (List[ShortPoint]): common generators for shuffle
            proof
        bay_gro_generator_creation_successful (List[List[bool]]): generator
            creation successful flags

    """
    def __init__(self, curve, N):
        """
        Args:
            curve: elliptic curve
            N (int): number of cards
        """
        self.N = N
        self.curve = curve

        # keygen
        self.secret_key = None
        self.public_key_share_own = None
        self.keygen_proof = None
        self.keygen_data = None
        self.public_key = None
        self.keygen_successful = None

        # bayer groth shuffle reference string
        self.n = None
        self.m = None
        self.bay_gro_secret_keys = None
        self.bay_gro_generator_shares_proofs = None
        self.bay_gro_generator_shares_data = None
        self.bay_gro_generators = None
        self.bay_gro_generator_creation_successful = None

    # init --------------------------------------------------------------------
    def init_cards_to_curve(self, list_of_cards=None):
        """Force card values to curve: card_list[i] = list_of_cards[i]*G,
        if no array is given as argument, array = [1,2,...,N]

        Args:
            list_of_cards (List[int]): list of integer representing the real
                cards

        Returns:
            card_list (List[ShortPoint]): forced card values
        """
        card_list = []

        if list_of_cards is None:
            list_of_cards = list(range(1, self.N + 1))

        for x in list_of_cards:
            card_list.append(self.curve.multiplication(
                x, self.curve.generator))

        return card_list

    def init_enc_card_points(self, cards_raw):
        """Mask forced card values with randomness k = 1 and public private key

        Args:
            cards_raw (List[ShortPoint]): forced card values

        Returns:
            cards_masked (List[ShortPoint, ShortPoint]): masked cards
        """
        k = [1]*self.N
        cards_masked = self.enc_cards(
            self.public_key, cards_raw, self.N, k)

        return cards_masked

    # keygen ------------------------------------------------------------------
    def key_generate(self):
        """Generate secret key and public key share and a proof that
        public_key_share_own = generator*secret_key
        """
        self.secret_key = self.curve.rand_gen.get_random_value()
        self.public_key_share_own = self.curve.multiplication(
            self.secret_key, self.curve.generator)
        var0 = proofs.PokProver(self.curve, self.curve.generator,
                                self.public_key_share_own,
                                self.secret_key)
        self.keygen_proof = var0.pok_nizk()

    def key_combine(self, *args):
        """Check proofs from other player, if all proofs are correct, combine
        all key shares to one public key pk = sum(pk_i)

        Args:
            *args ([List[ShortPoint], List[int]: public key shares and proofs
                from other players
        """
        self.keygen_successful = []
        for arg in args:
            var0 = proofs.PokVerifier(self.curve, self.curve.generator, arg[0])
            self.keygen_successful.append(var0.pok_nizk(arg[1][0], arg[1][1]))
        if all(self.keygen_successful):
            var0 = self.public_key_share_own
            for arg in args:
                var0 = self.curve.addition(var0, arg[0])
            self.public_key = var0

    # plaord ------------------------------------------------------------------
    def order_generate(self):
        """Generate random order_share_secret and calculate order_share_public
         = G*order_share_secret for random order

        Returns:
            order_share_secret (int), order_share_public (ShortPoint)
        """
        order_share_secret = self.curve.rand_gen.get_random_value()
        order_share_public = self.curve.multiplication(
            order_share_secret, self.curve.generator)

        return order_share_secret, order_share_public

    def order_combine(self, order_share_secret, *args):
        """Check for all args if a_i = G*k_i and calculates seed a = sum(a_i)
        mod q for i = 1,...,N

        Args:
            order_share_secret (int): own order share secret
            *args ([ShortPoint, int]): public and secret key shares
                from other players

        Returns:
            order_seed (int):
        """
        order_seed = order_share_secret
        for arg in args:
            if arg[0] == self.curve.multiplication(arg[1],
                                                   self.curve.generator):
                order_seed = (order_seed + arg[1]) % self.curve.order
            else:
                return None
        return order_seed

    def order_finalize(self, player_no, order_seed):
        """Use seed for PRNG to generate permutation using Fisher–Yates
        shuffle

        Args:
            player_no (int): number of players
            order_seed (int): order seed calculated by all players

        Returns:
            player_order (List[int]): player order (P1,...PN)
        """
        player_order = self.curve.rand_gen.get_random_permutation_seed(
            player_no, order_seed)

        return player_order

    @staticmethod
    def order_get_own(player_ids, id1, player_order):
        """Get order index by id

        Args:
            player_ids (List[str]): list of player IDs
            id1 (str): ID
            player_order (List[int]): player order (P1,...PN)

        Returns:
            int: index in player order (P1,...PN)
        """
        try:
            var0 = player_ids.index(id1)
            return player_order[var0]

        except Exception:
            return None

    # ck-init -----------------------------------------------------------------
    def shuffle_baygro_ck_generate(self, m, n):
        """Init m and n and generate n+1 generator shares for shuffle proof,
        G_i = bay_gro_secret_keys[i]*curve.generator for random
        bay_gro_secret_keys

        Args:
            n (int): number of columns for shuffle proof
            m (int): number of rows for shuffle proof
        """
        self.m = m
        self.n = n
        self.bay_gro_secret_keys = self.curve.rand_gen.get_random_array(n+1)
        self.bay_gro_generator_shares_proofs = []
        self.bay_gro_generator_shares_data = []
        for i in range(n+1):
            var0 = self.curve.multiplication(self.bay_gro_secret_keys[i],
                                             self.curve.generator)
            var1 = proofs.PokProver(self.curve,
                                    self.curve.generator,
                                    var0,
                                    self.bay_gro_secret_keys[i])
            self.bay_gro_generator_shares_proofs.append(var1.pok_nizk())
            self.bay_gro_generator_shares_data.append(var0)

    def shuffle_baygro_ck_combine(self, *args):
        """Verify and combine all generator shares from other players with own
        generator shares for commitment-keys

        Args:
            *args (List[List[ShortPoint], List[List[int]]): generator shares
            and proofs from other players
        """
        self.bay_gro_generator_creation_successful = []
        for arg in args:
            var0 = []
            for i in range(self.n + 1):
                var1 = proofs.PokVerifier(self.curve,  self.curve.generator,
                                          arg[0][i])
                var0.append(var1.pok_nizk(arg[1][i][0], arg[1][i][1]))
            self.bay_gro_generator_creation_successful.append(var0)

        if all(all(test) for test in
               self.bay_gro_generator_creation_successful):
            self.bay_gro_generators = []
            for i in range(self.n + 1):
                var0 = self.bay_gro_generator_shares_data[i]
                for arg in args:
                    var0 = self.curve.addition(var0, arg[0][i])
                self.bay_gro_generators.append(var0)

    # shuffle -----------------------------------------------------------------
    def shuffle_mix_remask_cards(self, enc_cards, size):
        """Shuffle cards due to a random permutation and remask them with
        random masking value

        Args:
            enc_cards (List[List[ShortPoint]]): masked cards
            size (int): number of masked cards

        Returns:
            List[List[List[ShortPoint]], List[int], List[int]]: shuffled
            and remasked cards, masking value, permutation
        """
        permuted_cards, pi = self.shuffle_permute_cards(enc_cards, size)
        shuffled_cards, rho = self.re_enc_cards(self.public_key,
                                                permuted_cards, size)

        return shuffled_cards, rho, pi

    def shuffle_permute_cards(self, enc_cards, size):
        """Generate random permutation and permute cards

        Args:
            enc_cards (List[List[ShortPoint]]): list of masked cards
            size (int): number of masked cards

        Returns:
            List[List[List[ShortPoint]], List[int]]: permuted cards,
            permutation
        """
        permutation = self.curve.rand_gen.get_random_permutation(size)

        permuted_cards = [enc_cards[permutation[x]] for x in range(size)]

        return permuted_cards, permutation

    def shuffle_baygro_prove(self, ciphers_in, ciphers_out, rho, pi):
        """Generate non-interactive shuffle proof

        Args:
            ciphers_in (List[List[ShortPoint]]): masked cards
            ciphers_out (List[List[ShortPoint]]): shuffled and re-masked cards
            rho (List[int]): random masking values used for re-masking
            pi (List[int]): random permutation used for remasking

        Returns:
            List: shuffle proof (see bayergroth.py)
        """
        prover = BayGroProver(self.N, self.m, self.n, self.curve, pi,
                              self.bay_gro_generators, self.public_key,
                              ciphers_in, ciphers_out, rho)
        return prover.nizk_prover()

    def shuffle_baygro_verify(self, ciphers_out, proof):
        """Verify non-interactive shuffle proof

        Args:
            ciphers_out (List[List[ShortPoint]]): shuffled and re-masked cards
            proof (List): shuffle proof

        Returns:
            bool: True if verification successful, False else
        """
        verifier = BayGroVerifier(self.N, self.m, self.n, self.curve,
                                  self.bay_gro_generators, self.public_key,
                                  ciphers_out)
        if all(verifier.nizk_verifier(proof)):
            return True
        else:
            return False

    def shuffle_shuffle_and_proof(self, ciphers_in):
        """Shuffle and re-mask cards and generate shuffle proof

        Args:
            ciphers_in (List[List[ShortPoint]]): masked cards

        Returns:
            List[List[ShortPoint]], List: shuffled and re-masked cards,
            shuffle proof
        """
        ciphers_out, rho, pi = self.shuffle_mix_remask_cards(
            ciphers_in, self.N)
        proof = self.shuffle_baygro_prove(ciphers_in, ciphers_out, rho, pi)
        return ciphers_out, proof

    # enc ---------------------------------------------------------------------
    def enc_cards(self, pubkey, forcedcard, size, k=None):
        """Mask curve points. enc_i = (c_i1, c_i2) = Enc(forced_card[i]; k[i])
        = (k[i]*curve.generator; forced_card[i] + k[i]*pubkey)

        Args:
            pubkey (ShortPoint): public key
            forcedcard (List[ShortPoint]): cards represented as curve points
            size (int): number of cards in forcedcard
            k (List[int]): list of with size elements of random masking
                values, if no array is given as argument, k is chosen randomly
        Returns:
            List[List[ShortPoint]]: masked curve points
        """
        # random integer for masking
        if k is None:
            k = self.curve.rand_gen.get_random_array(size)

        kha = [self.curve.multiplication(k[x], pubkey) for x in range(0, size)]
        enc = []
        for x in range(0, size):  # mask cards
            enc_a = self.curve.multiplication(k[x], self.curve.generator)
            enc_b = self.curve.addition(forcedcard[x], kha[x])
            enc.append([enc_a, enc_b])

        return enc

    # reenc -------------------------------------------------------------------
    def re_enc_cards(self, pubkey, enc_cards, size=0):
        """Re-mask card ciphers with random masking value,
        re_enc[i] = Enc(ZeroPoint; k[i]) + enc_cards[i]

        Args:
            pubkey (ShortPoint): public key
            enc_cards (List[List[ShortPoint]]): masked curve points
            size (int): number of cards in enc_cards

        Returns:
            List[List[List[ShortPoint]], List[int]]: remasked cards,
            random masking value
        """
        kha = []
        re_enc = []

        k = self.curve.rand_gen.get_random_array(size)

        for x in range(0, size):
            kha.append(self.curve.multiplication(k[x], pubkey))

        for x in range(0, size):  # re-mask cards c -> c''
            re_enc_a = self.curve.addition(enc_cards[x][0],
                                           self.curve.multiplication(k[x],
                                           self.curve.generator))
            re_enc_b = self.curve.addition(enc_cards[x][1], kha[x])
            re_enc.append([re_enc_a, re_enc_b])

        return re_enc, k

    # dec ---------------------------------------------------------------------
    def dec_generate(self, c):
        """Generate decryption share and proof for unmasking a cipher,
        d = secret_key*c[0]

        Args:
            c (List[ShortPoint, ShortPoint]): card cipher which should be
                unmasked

        Returns:
            List[List[ShortPoint], List[int]]: public key share, cipher
            share, decryption share and DLEQ(G,public_key_share_own, c[0],d).
        """
        d = self.curve.multiplication(self.secret_key, c[0])
        dleq = proofs.PeqProver(
            self.curve, self.curve.generator, self.public_key_share_own,
            c[0], d, self.secret_key)
        proof = dleq.peq_nizk()
        return [self.public_key_share_own, c[0], d], proof

    def dec_combine(self, c, *args):
        """Verify decryption share proof from all players and combine
        decryption shares to unmask card,
        Dec(c) = c_2 - sum(decryption shares)

        Args:
            c (List[ShortPoint, ShortPoint]): card cipher which should be
                unmasked
            *args (List[List[ShortPoint], List[int]]): public key share,
                cipher share, decryption share and DLEQ(G,public_key_share_own,
                c[0],d).

        Returns:
            ShortPoint: unmasked card/ elliptic curve point which represents
            one card.
        """
        var0 = self.curve.multiplication(self.secret_key, c[0])
        for arg in args:
            dleq = proofs.PeqVerifier(self.curve, self.curve.generator,
                                      arg[0][0], arg[0][1], arg[0][2])
            if dleq.peq_nizk(arg[1][0], arg[1][1]) is False:
                return None

            var0 = self.curve.addition(var0, arg[0][2])

        return self.curve.subtraction(c[1], var0)

    def dec_with_seckey(self, seckey, enc_cards, size=0):
        """Unmask cipher with secret key,
        Dec(enc_cards[i]) = enc[i][1] - enc[i][0]*seckey

        Args:
            seckey (int): secret key which public key is used for masking
            enc_cards (List[List[ShortPoint]]): card ciphers which should be
                unmasked
            size (int): number of card ciphers which should be unmasked

        Returns:
            List[ShortPoint]: unmasked ciphers
        """
        dec = []

        for x in range(0, size):
            dec_a = self.curve.multiplication(seckey, enc_cards[x][0])
            dec_b = enc_cards[x][1]
            dec.append(self.curve.subtraction(dec_b, dec_a))

        return dec

    @staticmethod
    def dec_index_raw_card(cards_raw, card):
        """Get index of card in cards_raw

        Args:
            cards_raw (List[ShortPoint]): List of elliptic curve points with
                all card points
            card (ShortPoint): elliptic curve point representing one card

        Returns:
            int: index of card in cards_raw
        """
        if card is None:
            return False

        try:
            return cards_raw.index(card)

        except Exception:
            return None

    # -------------------------------------------------------------------------
