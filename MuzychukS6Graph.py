def MuzychukS6Graph(n, d, phi='random', sigma='random'):
    '''
    Muzychuk strongly regular graph on n^d * ((n^d-1)/(n-1)+1) vertices, denoted
    S6 in [Mu07].
    n is a prime power.
    n is even or d is odd.
    phi defaults to random and is either:
        - 'random': phi_i are generated at random.
        - 'fixed': this will generate the same phi_i every time.
        - A dictionary describing the phi_i functions in Muzychuk's description:
          phi[(i, line)] should  be in {0,..., (n^d-1)/(n-1) - 1}. Note line is
          a tuple in ascending order. Also note each phi_i must be injective.
    sigma is similar, either:
        - 'random': sigma_ij are generated at random (from phi if phi is given).
        - 'fixed': this will generate the same sigma_ij every time the same
          phi_i are used.
        - A dictionary describing the sigma_ij: sigma[(i, j, n)] = m where i, j
          in some line in L, n in phi[(i, line)], m in phi[(j, line)]. Note a
          value must be given for all keys of this form. Also note sigma_ij
          must be (sigma_ji)^(-1).

    EXAMPLE:
        sage: MuzychukS6Graph(3, 3).is_strongly_regular(parameters=True)
        (378, 116, 34, 36)
    '''
    ### TO DO: optimise
    ###        add option to return phi, sigma? generate phi, sigma from seed? (int say?)

    from sage.combinat.designs.block_design import ProjectiveGeometryDesign
    from random import randrange
    from time import time

    assert is_even(n * (d-1)), 'n must be even or d must be odd'
    assert is_prime_power(n), 'n must be a prime power'
    assert phi or not sigma, 'sigma may only be given if phi is'
    t = time()

    #build L, L_i and the design
    m = int((n^d-1)/(n-1) + 1) #from m = p + 1, p = (n^d-1) / (n-1)
    L = graphs.CompleteGraph(m)
    L.delete_edges([(2*x, 2*x + 1) for x in range(m/2)])
    L_i = [0]*m
    for x in range(m):
        L_i[x] = [edge for edge in L.edges(labels=False) if x in edge]
    Design = ProjectiveGeometryDesign(d, d-1, GF(n, 'a'), point_coordinates=False)
    projBlocks = Design.blocks()
    atInf = projBlocks[-1]
    Blocks = [[x for x in block if x not in atInf] for block in projBlocks[:-1]]
    print 'finished preamble at %f (+%f)' % (time() - t, time() - t)
    t1 = time()

    #sort the hyperplanes into parallel classes
    ParClasses = [Blocks]
    while ParClasses[0]:
        nextHyp = ParClasses[0].pop()
        for C in ParClasses[1:]:
            listC = sum(C,[])
            for x in nextHyp:
                if x in listC:
                    break
            else:
                C.append(nextHyp)
                break
        else:
            ParClasses.append([nextHyp])
    del ParClasses[0]
    print 'finished ParClasses at %f (+%f)' % (time() - t, time() - t1)
    t1 = time()

    #build E^C_j
    E = {}
    v = n^d
    k = n ^ (d-1)
    ones = ones_matrix(n^d)
    for C in ParClasses:
        EC = matrix(QQ, n^d)
        for i in range(n^d):
            for j in range(n^d):
                for line in C:
                    if i in line and j in line:
                        EC[i, j] = 1/k
        EC -= ones/v
        E[tuple(C[0])] = EC
    print 'finished E at %f (+%f)' % (time() - t, time() - t1)
    t1 = time()

    #handle phi
    if phi == 'random':
        phi = {}
        for x in range(m):
            temp = range(len(ParClasses))
            for line in L_i[x]:
                rand = randrange(0, len(temp))
                phi[(x, line)] = temp.pop(rand)
    elif phi == 'fixed':
        phi = {}
        for x in range(m):
            val = 0
            for line in L_i[x]:
                phi[(x, line)] = val
                val+=1
    else:
        assert type(phi) == dict, 'phi must be a dictionary or\
        \'fixed\': alternatively, remove this argument and it will be\
        generated randomly'
        assert set(phi.keys()) == \
        set([(x, line) for x in range(m) for line in L_i[x]]), \
        'each phi_i must have domain L_i'
        for x in range(m):
            assert m - 2 == len(set([val
                for (key, val) in phi.items() if key[0] == x])), \
            'each phi_i must be injective'
        for val in phi.values():
            assert val in range(m-1), \
            'codomain should be {0,..., (n^d - 1)/(n - 1) - 1}'
    for x in range(m):
        for line in L_i[x]:
            phi[(x, line)] = ParClasses[phi[(x, line)]]
    print 'finished phi at %f (+%f)' % (time() - t, time() - t1)
    t1 = time()

    #handle sigma
    if sigma == 'random':
        sigma = {}
        for x in range(m):
            for line in L_i[x]:
                [i, j] = line
                temp = phi[(j, line)][:]
                for hyp in phi[(i, line)]:
                    rand = randrange(0, len(temp))
                    sigma[(i, j, tuple(hyp))] = temp[rand]
                    sigma[(j, i, tuple(temp[rand]))] = hyp
                    del temp[rand]
    elif sigma == 'fixed':
        sigma = {}
        for x in range(m):
            for line in L_i[x]:
                [i, j] = line
                temp = phi[(j, line)][:]
                for hyp in phi[(i, line)]:
                    val = temp.pop()
                    sigma[(i, j, tuple(hyp))] = val
                    sigma[(j, i, tuple(val))] = hyp
    else:
        assert type(sigma) == dict, \
        'sigma must be a dictionary or \'fixed\': alternatively, \
        remove this argument and it will be generated randomly'
        correctKeys =      [(line[0], line[1], n) for line in L.edges()
                            for n in range(len(ParClasses)) if
                            ParClasses[n] in phi[(line[1], line)]]
        correctKeys.extend([(line[1], line[0], n) for line in L.edges()
                            for n in range(len(ParClasses)) if
                            ParClasses[n] in phi[(line[0], line)]])
        assert set(sigma.keys()) == set(correct(Keys)), \
        'the keys in sigma must be \
        {(i, j, n) | i, j in line in L, and n in phi[i, line]}'
        for key in sigma.keys():
            assert key == sigma[(key[1], key[0], sigma[key])], \
            'sigma_ij must be (sigma_ji)^(-1)'
        sigma = dict(((i, j, tuple(ParClasses[x])), ParClasses[y])
                             for ((i, j, x), y) in sigma.items())
    print 'finished sigma at %f (+%f)' % (time() - t, time() - t1)
    t1 = time()

    #build V
    edges = [] ###how many? *m^2*n^2
    for (i, j) in L.edges(labels=False):
        for hyp in phi[(i, (i, j))]:
            for x in hyp:
                newEdges = [((i, x), (j, y))
                            for y in sigma[(i, j, tuple(hyp))]]
                edges.extend(newEdges)
    print 'finished edges at %f (+%f)' % (time() - t, time() - t1)
    t1 = time()
    V = Graph(edges)
    print 'finished V at %f (+%f)' % (time() - t, time() - t1)
    t1 = time()

    #build D_i, F_i and A_i
    D_i = [0]*m
    for x in range(m):
        D_i[x] = sum([E[tuple(phi[x, line][0])] for line in L_i[x]])
    F_i = [1 - D_i[x] - ones/v for x in range(m)]
    #as the sum of (1/v)*J_\Omega_i, D_i, F_i is identity
    A_i = [0]*m
    for x in range(m):
        A_i[x] = ((v-k)/v)*ones - k*F_i[x]
        #we know A_i = k''*(1/v)*J_\Omega_i + r''*D_i + s''*F_i,
        #and (k'', s'', r'') = (v - k, 0, -k)
    print 'finished D, F and A at %f (+%f)' % (time() - t, time() - t1)
    t1 = time()

    #add the edges of the graph of B to V
    for i in range(m):
        V.add_edges([((i, x), (i, y)) for x in range(n^d)
                     for y in range(n^d) if not A_i[i][(x, y)]])
    print 'finished at %f (+%f)' % ((time() - t), time() - t1)
    return V
