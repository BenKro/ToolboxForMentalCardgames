from ECCGT.eccwrapper import Fastecdsa as ECCobj
import fastecdsa.curve as curvelib
from ECCGT.gameset import GameSet
from ECCGT.toolbox import Toolbox
from copy import deepcopy

curve = ECCobj(curvelib.brainpoolP160r1)

player_ids = ["A", "B", "C"]
m = 4
n = 8
N = 32
M = 3

# Instantiate Toolbox objects -------------------------------------------------
P_Toolbox = [Toolbox(curve, N) for i in range(M)]

# Instantiate GameSet objects -------------------------------------------------
P_GameSet = [GameSet(M, N, player_ids, player_ids[i]) for i in range(M)]

# Order Key Share Generation --------------------------------------------------
for i in range(M):
    P_GameSet[i].order_share_secret, P_GameSet[i].order_share_public = \
        P_Toolbox[i].order_generate()

# Order Key Share Combination -------------------------------------------------
empt_list = []
for i in range(M):
    empt_list = []
    for j in range(M):
        if j != i:
            empt_list.append([P_GameSet[j].order_share_public,
                              P_GameSet[j].order_share_secret])
    P_GameSet[i].order_seed = \
        P_Toolbox[i].order_combine(P_GameSet[i].order_share_secret, *empt_list)

# Order Key Share Finalize ----------------------------------------------------
for i in range(M):
    P_GameSet[i].player_order = P_Toolbox[i].order_finalize(
        M, P_GameSet[i].order_seed)
    P_GameSet[i].player_own_order = P_Toolbox[i].order_get_own(
        player_ids, P_GameSet[i].own_id, P_GameSet[i].player_order)

# Key Share Generation --------------------------------------------------------
for i in range(M):
    P_Toolbox[i].key_generate()

# Key Share Combination -------------------------------------------------------
for i in range(M):
    empt_list = []
    for j in range(M):
        if j != i:
            empt_list.append(
                [P_Toolbox[j].public_key_share_own, P_Toolbox[j].keygen_proof])
    P_Toolbox[i].key_combine(*empt_list)

# Force Cards To Curve --------------------------------------------------------
for i in range(M):
    P_GameSet[i].cards_raw = P_Toolbox[i].init_cards_to_curve()

# Mask Cards ------------------------------------------------------------------
for i in range(M):
    P_GameSet[i].cards_masked = P_Toolbox[i].init_enc_card_points(
        P_GameSet[i].cards_raw)

# CK Share Generation ---------------------------------------------------------
for i in range(M):
    P_Toolbox[i].shuffle_baygro_ck_generate(m, n)

# CK Share Combination --------------------------------------------------------
for i in range(M):
    empt_list = []
    for j in range(M):
        if j != i:
            empt_list.append(
                [P_Toolbox[j].bay_gro_generator_shares_data,
                 P_Toolbox[j].bay_gro_generator_shares_proofs])
    P_Toolbox[i].shuffle_baygro_ck_combine(*empt_list)

# Shuffle Prove and Verify ----------------------------------------------------
ciphers_out, proof = P_Toolbox[0].shuffle_shuffle_and_proof(
    P_GameSet[0].cards_masked)

x = None
for i in range(M - 1):
    assert P_Toolbox[i+1].shuffle_baygro_verify(ciphers_out, proof)

for i in range(1, M):
    ciphers_out, proof = P_Toolbox[i].shuffle_shuffle_and_proof(ciphers_out)
    for j in range(M):
        if j != i:
            assert P_Toolbox[j].shuffle_baygro_verify(ciphers_out, proof)

for i in range(M):
    P_GameSet[i].cards_shuffled = ciphers_out

# Check Shuffle ---------------------------------------------------------------
for i in range(1, M):
    assert P_GameSet[0].cards_shuffled == P_GameSet[i].cards_shuffled

# Draw 10 Cards per Player ----------------------------------------------------
for i in range(M):
    for k in range(10):
        empt_list = []
        for j in range(M):
            if j != i:
                empt_list.append(P_Toolbox[j].dec_generate(P_GameSet[
                    j].cards_shuffled[k+i*10]))
                P_GameSet[j].cards_owner[k+i*10] = player_ids[i]

        card = P_Toolbox[i].dec_combine(
            P_GameSet[i].cards_shuffled[k+i*10], *empt_list)

        P_GameSet[i].cards_unmasked[k+i*10] = P_Toolbox[
            i].dec_index_raw_card(P_GameSet[i].cards_raw, card)
        P_GameSet[i].cards_owner[k+i*10] = player_ids[i]

# Draw the 2 Skat Cards per Player --------------------------------------------
for k in range(30, 32):
    empt_list = []
    for j in range(M):
        if j != 0:
            empt_list.append(P_Toolbox[j].dec_generate(P_GameSet[
                j].cards_shuffled[k]))
            P_GameSet[j].cards_owner[k] = "Skat"

    card = P_Toolbox[0].dec_combine(
        P_GameSet[0].cards_shuffled[k], *empt_list)

    P_GameSet[0].cards_unmasked[k] = P_Toolbox[
        0].dec_index_raw_card(P_GameSet[0].cards_raw, card)
    P_GameSet[0].cards_owner[k] = "Skat"

# Open All Cards ----------------------------------------------------------
for k in range(N):
    empt_list = []
    for i in range(M):
        empt_list.append(P_Toolbox[i].dec_generate(P_GameSet[
            i].cards_shuffled[k]))

    for i in range(M):
        testlist = deepcopy(empt_list)
        testlist.pop(i)
        card = P_Toolbox[i].dec_combine(
            P_GameSet[i].cards_shuffled[k], *testlist)

        var = P_Toolbox[i].dec_index_raw_card(P_GameSet[i].cards_raw, card)

# Print Result ----------------------------------------------------------------
player_ids2 = ["Alice", "Bob", "Charlie"]
print("-----------------------------------------")
for i in range(M):
    print(player_ids2[i] + "-CARDS:")
    print(P_GameSet[i].cards_unmasked)
    print(player_ids2[i] + "-OWNER:")
    print(P_GameSet[i].cards_owner)
    print("-----------------------------------------")

# Löschen der Hilfsvariablen --------------------------------------------------
del card, testlist, var, i, j, k
del ciphers_out, empt_list, x, player_ids2, proof, player_ids
pass
