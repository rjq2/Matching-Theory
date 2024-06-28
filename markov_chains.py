import numpy as np

def comb(n, m):
    numbers = np.arange(m + 1)
    permutations = []

    if n == 2:
        for i in numbers:
            for j in numbers:
                permut = [i, j]
                if sum(np.array(permut) == 0) <= n and sum(np.array(permut) == 1) <= 1 and sum(np.array(permut) == 2) <= 1 and sum(np.array(permut) == 3) <= 1 and sum(np.array(permut) == 4) <= 1 and sum(np.array(permut) == 5) <= 1 and sum(np.array(permut) == 6) <= 1:
                    permutations.append(permut)
    
    elif n == 3:
        for i in numbers:
            for j in numbers:
                for k in numbers:
                    permut = [i, j, k]
                    if sum(np.array(permut) == 0) <= n and sum(np.array(permut) == 1) <= 1 and sum(np.array(permut) == 2) <= 1 and sum(np.array(permut) == 3) <= 1 and sum(np.array(permut) == 4) <= 1 and sum(np.array(permut) == 5) <= 1 and sum(np.array(permut) == 6) <= 1:
                        permutations.append(permut)
    
    elif n == 4:
        for i in numbers:
            for j in numbers:
                for k in numbers:
                    for l in numbers:
                        permut = [i, j, k, l]
                        if sum(np.array(permut) == 0) <= n and sum(np.array(permut) == 1) <= 1 and sum(np.array(permut) == 2) <= 1 and sum(np.array(permut) == 3) <= 1 and sum(np.array(permut) == 4) <= 1 and sum(np.array(permut) == 5) <= 1 and sum(np.array(permut) == 6) <= 1:
                            permutations.append(permut)
    
    elif n == 5:
        for i in numbers:
            for j in numbers:
                for k in numbers:
                    for l in numbers:
                        for o in numbers:
                            permut = [i, j, k, l, o]
                            if sum(np.array(permut) == 0) <= n and sum(np.array(permut) == 1) <= 1 and sum(np.array(permut) == 2) <= 1 and sum(np.array(permut) == 3) <= 1 and sum(np.array(permut) == 4) <= 1 and sum(np.array(permut) == 5) <= 1 and sum(np.array(permut) == 6) <= 1:
                                permutations.append(permut)
    
    elif n == 6:
        for i in numbers:
            for j in numbers:
                for k in numbers:
                    for l in numbers:
                        for o in numbers:
                            for p in numbers:
                                permut = [i, j, k, l, o, p]
                                if sum(np.array(permut) == 0) <= n and sum(np.array(permut) == 1) <= 1 and sum(np.array(permut) == 2) <= 1 and sum(np.array(permut) == 3) <= 1 and sum(np.array(permut) == 4) <= 1 and sum(np.array(permut) == 5) <= 1 and sum(np.array(permut) == 6) <= 1:
                                    permutations.append(permut)

    return np.array(permutations)

from itertools import permutations

def possible_pref(n, rep, n_f):
    v = np.arange(1, n_f + 1)
    
    # Generate all permutations of v
    possible_pref = np.array(list(permutations(v))).T
    
    m = possible_pref.shape[1]
    
    results_iid = np.zeros((rep, n), dtype=int)
    pref_domain = np.zeros((n_f, n, rep), dtype=int)
    
    counter = 0
    for i in range(rep):
        results_iid[i, :] = np.random.randint(1, m + 1, size=n)
        pref_domain[:, :, counter] = possible_pref[:, results_iid[i, :]-1]
        counter += 1
    
    return pref_domain

def possible_pref_f(n, rep, n_f):
    v = np.arange(1, n + 1)
    
    # Generate all permutations of v
    possible_pref = np.array(list(permutations(v))).T
    
    m = possible_pref.shape[1]
    
    results_iid = np.zeros((rep, n_f), dtype=int)
    pref_domain = np.zeros((n, n_f, rep), dtype=int)
    
    counter = 0
    for i in range(rep):
        results_iid[i, :] = np.random.randint(1, m + 1, size=n_f)
        pref_domain[:, :, counter] = possible_pref[:, results_iid[i, :]-1]
        counter += 1
    
    return pref_domain


def markov_transition_biro(pref_1, pref_2, MatchingM, actual_match_M):
    n = pref_1.shape[1]
    n_f = pref_2.shape[1]
    MatchingF = np.zeros(n_f, dtype=int)

    # Create MatchingF
    for j in range(n_f):
        if j + 1 not in MatchingM:
            MatchingF[j] = 0
        else:
            column = np.where(MatchingM == j + 1)[0]
            if column.size > 0:
                MatchingF[j] = column[0] + 1
            else:
                MatchingF[j] = 0

    # Create b_m
    b_m = np.zeros((n, n_f), dtype=int)
    for i in range(n):
        for j in range(n_f):
            if MatchingM[i] == 0:
                b_m[i, j] = 1
            else:
                if MatchingM[i] != j + 1 and pref_1[j, i] > pref_1[MatchingM[i] - 1, i]:
                    b_m[i, j] = 1
                else:
                    b_m[i, j] = 0

    # Create b_f
    b_f = np.zeros((n_f, n), dtype=int)
    for i in range(n_f):
        for j in range(n):
            if MatchingF[i] == 0:
                b_f[i, j] = 1
            else:
                if MatchingF[i] != j + 1 and pref_2[j, i] > pref_2[MatchingF[i] - 1, i]:
                    b_f[i, j] = 1
                else:
                    b_f[i, j] = 0

    blocking_pair = b_m & b_f.T
    blocking_pair_1 = blocking_pair.copy()

    new_match_M_cell = np.zeros((1, n), dtype=int)

    if not np.any(blocking_pair_1 == 1):
        if np.array_equal(MatchingM, actual_match_M):
            p = 1
        else:
            p = 0
    else:
        index = 0
        for i in range(n):
            for j in range(n_f):
                if blocking_pair_1[i, j] == 1:
                    new_match_M = MatchingM.copy()
                    column = np.where(MatchingM == j + 1)[0]
                    if column.size > 0:
                        new_match_M[column[0]] = 0
                    new_match_M[i] = j + 1
                    if index < new_match_M_cell.shape[0]:
                        new_match_M_cell[index, :] = new_match_M
                    else:
                        new_match_M_cell = np.vstack([new_match_M_cell, new_match_M])
                    index += 1

        coincidence = np.any(np.all(new_match_M_cell == actual_match_M, axis=1))
        n_coincidence = np.sum(coincidence)

        p = n_coincidence / new_match_M_cell.shape[0]

    blocking_pair_4 = blocking_pair.copy()



    bloq_f_1 = blocking_pair.T
    bloq_f_2 = pref_2.T

    bloq_m_1 = blocking_pair
    bloq_m_2 = pref_1.T

    for i in range(n):
        if np.sum(blocking_pair_4[i, :] == 1) > 1:
            blocking_indices = np.where(bloq_m_1[i, :] == 1)[0]  
            m_payoff = bloq_m_2[i, blocking_indices]
            f_payoff = bloq_f_2[blocking_indices, i]
            total = m_payoff + f_payoff
            a = blocking_indices
            be = a[np.where(total == np.max(total))[0]]
            if len(be) >= 2:
                e = len(be)
                blocking_pair_4[i, :] = 0
                blocking_pair_4[i, be[:e]] = 1
            else:
                a = blocking_indices
                be = a[np.where(total == np.max(total))[0]]
                blocking_pair_4[i, :] = 0
                blocking_pair_4[i, be] = 1

    new_match_M_cell = np.zeros((1, n), dtype=int)

    if not np.any(blocking_pair_4 == 1):
        if np.array_equal(MatchingM, actual_match_M):
            p4 = 1
        else:
            p4 = 0
    else:
        index = 0
        for i in range(n):
            for j in range(n_f):
                if blocking_pair_4[i, j] == 1:
                    new_match_M = MatchingM.copy()
                    column = np.where(MatchingM == j + 1)[0]
                    if column.size > 0:
                        new_match_M[column[0]] = 0
                    new_match_M[i] = j + 1
                    if index < new_match_M_cell.shape[0]:
                        new_match_M_cell[index, :] = new_match_M
                    else:
                        new_match_M_cell = np.vstack([new_match_M_cell, new_match_M])
                    index += 1

        coincidence = np.any(np.all(new_match_M_cell == actual_match_M, axis=1))
        n_coincidence = np.sum(coincidence)

        p4 = n_coincidence / new_match_M_cell.shape[0]


    blocking_pair_5 = blocking_pair.T.copy()
    
    bloq_f_1 = blocking_pair.T
    bloq_f_2 = pref_2.T

    bloq_m_1 = blocking_pair
    bloq_m_2 = pref_1
    
    for i in range(n_f):
        if np.sum(blocking_pair_5[i, :] == 1) > 1:
            blocking_indices = np.where(bloq_f_1[i, :] == 1)[0]
            if np.max(blocking_indices) < bloq_m_2.shape[1]:  # Check index bounds for bloq_m_2
                m_payoff = bloq_m_2[i, blocking_indices]
                if np.max(blocking_indices) < bloq_f_2.shape[0]:  # Check index bounds for bloq_f_2
                    f_payoff = bloq_f_2[i, blocking_indices]
                    total = m_payoff + f_payoff
                    a = np.where(bloq_f_1[i, :] == 1)[0]
                    be = a[np.where(total == np.max(total))[0]]
                    if len(be) >= 2:
                        e = len(be)
                        blocking_pair_5[i, :] = 0
                        blocking_pair_5[i, be[:e]] = 1
                    else:
                        a = np.where(bloq_f_1[i, :] == 1)[0]
                        be = a[np.where(total == np.max(total))[0]]
                        blocking_pair_5[i, :] = 0
                        blocking_pair_5[i, be] = 1

    blocking_pair_5 = blocking_pair_5.T
    new_match_M_cell = np.zeros((1, n), dtype=int)

    if not np.any(blocking_pair_5 == 1):
        if np.array_equal(MatchingM, actual_match_M):
            p5 = 1
        else:
            p5 = 0
    else:
        index = 0
        for i in range(n):
            for j in range(n_f):
                if blocking_pair_5[i, j] == 1:
                    new_match_M = MatchingM.copy()
                    column = np.where(MatchingM == j + 1)[0]
                    if column.size > 0:
                        new_match_M[column[0]] = 0
                    new_match_M[i] = j + 1
                    if index < new_match_M_cell.shape[0]:
                        new_match_M_cell[index, :] = new_match_M
                    else:
                        new_match_M_cell = np.vstack([new_match_M_cell, new_match_M])
                    index += 1

        coincidence = np.any(np.all(new_match_M_cell == actual_match_M, axis=1))
        n_coincidence = np.sum(coincidence)

        p5 = n_coincidence / new_match_M_cell.shape[0]

    blocking_pair_2 = blocking_pair.copy()

    for i in range(n):
        if np.sum(blocking_pair_2[i, :] == 1) > 1:
            a = np.where(blocking_pair_2[i, :] == 1)[0]
            if n_f == 2:
                blocking_pair_2[i, :] = 0
                blocking_pair_2[i, np.argmax(pref_1[:, i])] = 1
            elif n_f == 3:
                if len(a) == 2:
                    if pref_1[a[0], i] > pref_1[a[1], i]:
                        blocking_pair_2[i, a[0]] = 1
                        blocking_pair_2[i, a[1]] = 0
                    else:
                        blocking_pair_2[i, a[1]] = 1
                        blocking_pair_2[i, a[0]] = 0
                elif len(a) == 3:
                    blocking_pair_2[i, :] = 0
                    blocking_pair_2[i, np.argmax(pref_1[:, i])] = 1
            elif n_f == 4:
                if len(a) == 2:
                    if pref_1[a[0], i] > pref_1[a[1], i]:
                        blocking_pair_2[i, a[0]] = 1
                        blocking_pair_2[i, a[1]] = 0
                    else:
                        blocking_pair_2[i, a[1]] = 1
                        blocking_pair_2[i, a[0]] = 0
                elif len(a) == 3:
                    if pref_1[a[0], i] > pref_1[a[1], i] and pref_1[a[0], i] > pref_1[a[2], i]:
                        blocking_pair_2[i, :] = 0
                        blocking_pair_2[i, a[0]] = 1
                    elif pref_1[a[1], i] > pref_1[a[0], i] and pref_1[a[1], i] > pref_1[a[2], i]:
                        blocking_pair_2[i, :] = 0
                        blocking_pair_2[i, a[1]] = 1
                    elif pref_1[a[2], i] > pref_1[a[0], i] and pref_1[a[2], i] > pref_1[a[1],i] :
                        blocking_pair_2[i, :] = 0
                        blocking_pair_2[i, a[2]] = 1
                elif len(a) == 4:
                    blocking_pair_2[i, :] = 0
                    blocking_pair_2[i, np.argmax(pref_1[:, i])] = 1

    new_match_M_cell = np.zeros((1, n), dtype=int)

    if not np.any(blocking_pair_2 == 1):
        if np.array_equal(MatchingM, actual_match_M):
            p2 = 1
        else:
            p2 = 0
    else:
        index = 0
        for i in range(n):
            for j in range(n_f):
                if blocking_pair_2[i, j] == 1:
                    new_match_M = MatchingM.copy()
                    column = np.where(MatchingM == j + 1)[0]
                    if column.size > 0:
                        new_match_M[column[0]] = 0
                    new_match_M[i] = j + 1
                    if index < new_match_M_cell.shape[0]:
                        new_match_M_cell[index, :] = new_match_M
                    else:
                        new_match_M_cell = np.vstack([new_match_M_cell, new_match_M])
                    index += 1

        coincidence = np.any(np.all(new_match_M_cell == actual_match_M, axis=1))
        n_coincidence = np.sum(coincidence)

        p2 = n_coincidence / new_match_M_cell.shape[0]

    blocking_pair_3 = blocking_pair.copy()
    pref_2 = pref_2.T

    for i in range(n_f):
        if np.sum(blocking_pair_3[:, i] == 1) > 1:
            a = np.where(blocking_pair_3[:, i] == 1)[0]
            if n == 2:
                blocking_pair_3[:, i] = 0
                blocking_pair_3[np.argmax(pref_2[:, i]), i] = 1
            elif n == 3:
                if len(a) == 2:
                    if pref_2[i, a[0]] > pref_2[i, a[1]]:
                        blocking_pair_3[a[0], i] = 1
                        blocking_pair_3[a[1], i] = 0
                    else:
                        blocking_pair_3[a[1], i] = 1
                        blocking_pair_3[a[0], i] = 0
                elif len(a) == 3:
                    blocking_pair_3[:, i] = 0
                    blocking_pair_3[np.argmax(pref_2[i, :]), i] = 1
            elif n == 4:
                if len(a) == 2:
                    if pref_2[i, a[0]] > pref_2[i, a[1]]:
                        blocking_pair_3[a[0], i] = 1
                        blocking_pair_3[a[1], i] = 0
                    else:
                        blocking_pair_3[a[1], i] = 1
                        blocking_pair_3[a[0], i] = 0
                elif len(a) == 3:
                    if pref_2[i, a[0]] > pref_2[i, a[1]] and pref_2[i, a[0]] > pref_2[i, a[2]]:
                        blocking_pair_3[:, i] = 0
                        blocking_pair_3[a[0], i] = 1
                    elif pref_2[i, a[1]] > pref_2[i, a[0]] and pref_2[i, a[1]] > pref_2[i, a[2]]:
                        blocking_pair_3[:, i] = 0
                        blocking_pair_3[a[1], i] = 1
                    elif pref_2[i, a[2]] > pref_2[i, a[0]] and pref_2[i, a[2]] > pref_2[i, a[1]]:
                        blocking_pair_3[:, i] = 0
                        blocking_pair_3[a[2], i] = 1
                elif len(a) == 4:
                    blocking_pair_3[:, i] = 0
                    blocking_pair_3[np.argmax(pref_2[i, :]), i] = 1

    new_match_M_cell = np.zeros((1, n), dtype=int)
    pref_2 = pref_2.T
    if not np.any(blocking_pair_3 == 1):
        if np.array_equal(MatchingM, actual_match_M):
            p3 = 1
        else:
            p3 = 0
    else:
        index = 0
        for i in range(n):
            for j in range(n_f):
                if blocking_pair_3[i, j] == 1:
                    new_match_M = MatchingM.copy()
                    column = np.where(MatchingM == j + 1)[0]
                    if column.size > 0:
                        new_match_M[column[0]] = 0
                    new_match_M[i] = j + 1
                    if index < new_match_M_cell.shape[0]:
                        new_match_M_cell[index, :] = new_match_M
                    else:
                        new_match_M_cell = np.vstack([new_match_M_cell, new_match_M])
                    index += 1

        coincidence = np.any(np.all(new_match_M_cell == actual_match_M, axis=1))
        n_coincidence = np.sum(coincidence)

        p3 = n_coincidence / new_match_M_cell.shape[0]

    return p, p2, p3, p4, p5

def markov_absorb(M, all_matchings_M):
    # Filas con uno en la diagonal de M
    filas_con_uno = np.diag(M) == 1
    
    # Estables
    stables = all_matchings_M[filas_con_uno, :]
    
    # Obtener filas con valor '1' y eliminar esas filas de la matriz original
    filas_con_uno_valores = M[filas_con_uno, :]
    M_2 = M[~filas_con_uno, :]
    Q = M_2[:, ~filas_con_uno]
    R = np.zeros((len(Q), len(M)))
    
    for i in range(len(M)):
        if np.any(M[i, i] == 1):
            R[:, i] = M_2[:, i]
    
    # Columnas sin ceros
    columnas_sin_ceros = ~np.all(R == 0, axis=0)
    R = R[:, columnas_sin_ceros]
    
    zero = np.zeros((len(R[0]), len(Q)))
    identity = np.eye(len(R[0]))
    
    # Forma Canónica
    P = np.block([[Q, R], [zero, identity]])
    
    # Matriz N del ppt
    N = np.linalg.inv((np.eye(len(Q)) - Q))
    
    # Tiempos esperados de absorción
    t = np.dot(N, np.ones(len(N)))
    t = np.mean(t)
    t_2 = np.dot(N, np.ones(len(N)))
    
    # Probabilidad de Absorción
    B = np.dot(N, R)
    
    # Número de estables
    variability = np.any(stables != stables[0, :], axis=0)
    n_h_stables = np.sum(variability)
    
    return t, n_h_stables, B, t_2

def markov_sim(pref_1, pref_2, all_matchings_M):
    m = len(all_matchings_M)
    M = np.zeros((m, m))
    M_m = np.zeros((m, m))
    M_f = np.zeros((m, m))
    M_e = np.zeros((m, m))
    M_e_f = np.zeros((m, m))
    
    for i in range(m):
        for j in range(m):
            M[i,j], M_m[i,j], M_f[i,j], M_e[i,j], M_e_f[i,j] = markov_transition_biro(pref_1, pref_2, all_matchings_M[i,:], all_matchings_M[j,:])
    
    filas_con_uno = np.diag(M) == 1
    stables = all_matchings_M[filas_con_uno, :]
    
    t = np.zeros(5)
    n_h_stables = np.zeros(1)
    
    t[0], n_h_stables, _, _ = markov_absorb(M, all_matchings_M)
    t[1], _ , _, _ = markov_absorb(M_m, all_matchings_M)
    t[2], _ , _, _,= markov_absorb(M_f, all_matchings_M)
    t[3], _ , _, _, = markov_absorb(M_e, all_matchings_M)
    t[4], _ , _, _, = markov_absorb(M_e_f, all_matchings_M)
    return t, n_h_stables, M, M_m ,M_f, M_e, M_e_f