from ECCGT.eccwrapper import Fastecdsa as ECCobj
import fastecdsa.curve as curvelib
from ECCGT.gameset import GameSet
from ECCGT.toolbox import Toolbox
import time
from copy import deepcopy

curve = ECCobj(curvelib.secp256k1)

player_ids = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
m = 4
n = 13
N = 52
P = 10

# Instantiate Toolbox objects -------------------------------------------------
start = time.time()
P_Toolbox = [Toolbox(curve, N) for i in range(P)]
end = time.time()
print("Init Toolbox", (end - start))

# Instantiate GameSet objects -------------------------------------------------
start = time.time()
P_GameSet = [GameSet(P, N, player_ids, player_ids[i]) for i in range(P)]
end = time.time()
print("Init GameSet", (end - start))

# Order Key Share Generation --------------------------------------------------
for i in range(P):
    start = time.time()
    P_GameSet[i].order_share_secret, P_GameSet[i].order_share_public = \
        P_Toolbox[i].order_generate()
    end = time.time()
    print("Order-Get Key", (end - start))

# Order Key Share Combination -------------------------------------------------
empt_list = []
for i in range(P):
    empt_list = []
    for j in range(P):
        if j != i:
            empt_list.append([P_GameSet[j].order_share_public,
                              P_GameSet[j].order_share_secret])
    start = time.time()
    P_GameSet[i].order_seed = \
        P_Toolbox[i].order_combine(P_GameSet[i].order_share_secret, *empt_list)
    end = time.time()
    print("Order-Combine Key", (end - start))

# Order Key Share Finalize ----------------------------------------------------
for i in range(P):
    start = time.time()
    P_GameSet[i].player_order = P_Toolbox[i].order_finalize(
        P, P_GameSet[i].order_seed)
    P_GameSet[i].player_own_order = P_Toolbox[i].order_get_own(
        player_ids, P_GameSet[i].own_id, P_GameSet[i].player_order)
    end = time.time()
    print("Order-Finalize Key", (end - start))

# Key Share Generation --------------------------------------------------------
for i in range(P):
    start = time.time()
    P_Toolbox[i].key_generate()
    end = time.time()
    print("KeyGen-Generate Key", (end - start))

# Key Share Combination -------------------------------------------------------
for i in range(P):
    empt_list = []
    for j in range(P):
        if j != i:
            empt_list.append(
                [P_Toolbox[j].public_key_share_own, P_Toolbox[j].keygen_proof])
    start = time.time()
    P_Toolbox[i].key_combine(*empt_list)
    end = time.time()
    print("KeyGen-Combine/Finalize Key", (end - start))

# CK Share Generation ---------------------------------------------------------
for i in range(P):
    start = time.time()
    P_Toolbox[i].shuffle_baygro_ck_generate(m, n)
    end = time.time()
    print("CK-Generate Key", (end - start))

# CK Share Combination --------------------------------------------------------
for i in range(P):
    empt_list = []
    for j in range(P):
        if j != i:
            empt_list.append(
                [P_Toolbox[j].bay_gro_generator_shares_data,
                 P_Toolbox[j].bay_gro_generator_shares_proofs])
    start = time.time()
    P_Toolbox[i].shuffle_baygro_ck_combine(*empt_list)
    end = time.time()
    print("CK-Combine/Finalize Key", (end - start))

# Force Cards To Curve --------------------------------------------------------
for i in range(P):
    start = time.time()
    P_GameSet[i].cards_raw = P_Toolbox[i].init_cards_to_curve()
    end = time.time()
    print("Init-Cards To Curve", (end - start))

# Mask Cards ------------------------------------------------------------------
for i in range(P):
    start = time.time()
    P_GameSet[i].cards_masked = P_Toolbox[i].init_enc_card_points(
        P_GameSet[i].cards_raw)
    end = time.time()
    print("Init-Mask Cards", (end - start))

# Shuffle Prove and Verify ----------------------------------------------------
start = time.time()
ciphers_out, proof = P_Toolbox[0].shuffle_shuffle_and_proof(
    P_GameSet[0].cards_masked)
end = time.time()
print("Shuffle-Proof", (end - start))
x = None
for i in range(P-1):
    start = time.time()
    assert P_Toolbox[i+1].shuffle_baygro_verify(ciphers_out, proof)
    end = time.time()
    print("Shuffle-Verify", (end - start))
for i in range(1, P):
    start = time.time()
    ciphers_out, proof = P_Toolbox[i].shuffle_shuffle_and_proof(ciphers_out)
    end = time.time()
    print("Shuffle-Proof", (end - start))
    for j in range(P):
        if j != i:
            start = time.time()
            assert P_Toolbox[j].shuffle_baygro_verify(ciphers_out, proof)
            end = time.time()
            print("Shuffle-Verify", (end - start))

for i in range(P):
    P_GameSet[i].cards_shuffled = ciphers_out

# Check Shuffle ---------------------------------------------------------------
for i in range(1, P):
    assert P_GameSet[0].cards_shuffled == P_GameSet[i].cards_shuffled

# Draw 2 Cards per Player -----------------------------------------------------
for i in range(P):
    for k in range(2):
        empt_list = []
        for j in range(P):
            if j != i:
                start = time.time()
                empt_list.append(P_Toolbox[j].dec_generate(P_GameSet[
                    j].cards_shuffled[k+i*2]))
                end = time.time()
                print("Draw-Generate Dec", (end - start))
                P_GameSet[j].cards_owner[k+i*2] = player_ids[i]
        start = time.time()
        card = P_Toolbox[i].dec_combine(
            P_GameSet[i].cards_shuffled[k+i*2], *empt_list)
        end = time.time()
        print("Draw-Combine Dec", (end - start))
        start = time.time()
        P_GameSet[i].cards_unmasked[k+i*2] = P_Toolbox[
            i].dec_index_raw_card(P_GameSet[i].cards_raw, card)
        end = time.time()
        print("Draw-Index Dec", (end - start))
        P_GameSet[i].cards_owner[k+i*2] = player_ids[i]

# Open 5 Table Cards ----------------------------------------------------------
for k in range(P*2, P*2+5):
    empt_list = []
    for i in range(P):
        start = time.time()
        empt_list.append(P_Toolbox[i].dec_generate(P_GameSet[
            i].cards_shuffled[k]))
        end = time.time()
        print("Open-Generate Dec", (end - start))

    for i in range(P):
        testlist = deepcopy(empt_list)
        testlist.pop(i)
        start = time.time()
        card = P_Toolbox[i].dec_combine(
            P_GameSet[i].cards_shuffled[k], *testlist)
        end = time.time()
        print("Open-Combine Dec", (end - start))
        start = time.time()
        P_GameSet[i].cards_unmasked[k] = P_Toolbox[i].dec_index_raw_card(
            P_GameSet[i].cards_raw, card)
        end = time.time()
        print("Open-Index Dec", (end - start))
        P_GameSet[i].cards_owner[k] = "Open"

# Print Result ----------------------------------------------------------------
for i in range(P):
    print(P_GameSet[i].cards_unmasked)
    print(P_GameSet[i].cards_owner)

pass
