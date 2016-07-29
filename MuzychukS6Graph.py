def MuzychukS6Graph(n, d, List=False, phi=None, sigma=None):
    '''
    Muzychuk strongly regular graph on n^d * ((n^d-1)/(n-1)+1) vertices, denoted S6 in [Mu07].
    n is a prime power.
    n is even or d is odd.
    phi is either:
        - None: phi_i are generated at random.
        - 'fixed': this will generate the same phi_i every time.
        - A dictionary describing the phi_i functions in Muzychuk's description: phi[(i, line)] should  be in {0,..., (n^d-1)/(n-1) - 1}. Note line is a tuple in ascending order. Also note each phi_i must be injective.
    sigma is similar, either:
        - None: sigma_ij are generated at random (from phi if phi is given).
        - 'fixed': this will generate the same sigma_ij every time the same phi_i are used.
        - A dictionary describing the sigma_ij: sigma[(i, j, n)] = m where i, j in some line in L, n in phi[(i, line)], m in phi[(j, line)]. Note a value must be given for all keys of this form. Also note sigma_ij must be (sigma_ji)^(-1).

    EXAMPLE:
        sage: MuzychukS6Graph(3, 3).is_strongly_regular(parameters=True)
        (378, 116, 34, 36)
    '''
    ### TO DO: optimise
    ###        add option to return phi, sigma?

    from sage.combinat.designs.block_design import AffineGeometryDesign
    from random import randrange
    from time import time

    assert is_even(n * (d-1)), 'n must be even or d must be odd'
    assert is_prime_power(n), 'n must be a prime power'
    assert phi or not sigma, 'sigma may only be given if phi is'
    t = time()

    #build L, L_i and the affine design
    m = int((n^d-1)/(n-1) + 1) #from m = p + 1, p = (n^d-1) / (n-1)
    L = graphs.CompleteGraph(m)
    L.delete_edges([(2*x, 2*x + 1) for x in range(m/2)])
    L_i = [0]*m
    for x in range(m):
        L_i[x] = [edge for edge in L.edges(labels=False) if x in edge]
    Design = AffineGeometryDesign(d, d-1, n)
    print 'finished preamble at %f' % (time() - t)

    #sort the hyperplanes into parallel classes
    ParClasses = [Design.blocks()]
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
    print 'finished ParClasses at %f' % (time() - t)

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
    print 'finished E at %f' % (time() - t)

    #If phi is given as an argument, check it is of the correct form
    if phi:
        phiRand = False
        if phi == 'fixed':
            phi = {}
            for x in range(m):
                val = 0
                for line in L_i[x]:
                    phi[(x, line)] = val
                    val+=1
        else:
            assert type(phi) == dict, 'phi must be a dictionary or \'fixed\': alternatively, remove this argument and it will be generated randomly'
            assert set(phi.keys()) == set([(x, line) for x in range(m) for line in L_i[x]]), 'each phi_i must have domain L_i'
            for x in range(m):
                assert m - 2 == len(set([val for (key, val) in phi.items() if key[0] == x])), 'each phi_i must be injective'
            for val in phi.values():
                assert val in range(m-1), 'codomain should be {0,..., (n^d - 1)/(n - 1) - 1}'
        #now put phi in the form we need
        for x in range(m):
            for line in L_i[x]:
                phi[(x, line)] = ParClasses[phi[(x, line)]]
    else:
        phiRand = True
    print 'finished phi checking at %f' % (time() - t)

    #handle sigma if it's given
    if sigma:
        sigmaRand = False
        if sigma == 'fixed':
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
            assert type(sigma) == dict, 'sigma must be a dictionary or \'fixed\': alternatively, remove this argument and it will be generated randomly'
            correctKeys =     [(line[0], line[1], n) for line in L.edges() for n in range(len(ParClasses)) if ParClasses[n] in phi[(line[1], line)]]
            correctKey.extend([(line[1], line[0], n) for line in L.edges() for n in range(len(ParClasses)) if ParClasses[n] in phi[(line[0], line)]])
            assert set(sigma.keys()) == set(correct(Keys)), 'the keys in sigma must be {(i, j, n) | i, j in line in L, and n in phi[i, line]}'
            for key in sigma.keys():
                assert key == sigma[(key[1], key[0], sigma[key])], 'sigma_ij must be (sigma_ji)^(-1)'
            sigma = dict(((i, j, ParClasses[x]), ParClasses[y]) for ((i, j, x), y) in sigma.items())
    else:
        sigmaRand = True
    print 'finished sigma checking at %f' % (time() - t)

    #set up a loop to build a list of graphs
    if List is False:
        N = 1
    else:
        N = List
    Glist = [0]*N
    for I in range(N):
        #build phi randomly if required
        if phiRand:
            phi = {}
            for x in range(m):
                temp = ParClasses[:]
                for line in L_i[x]:
                    rand = randrange(0, len(temp))
                    phi[(x, line)] = temp.pop(rand)
        print 'finished phi at %f' % (time() - t)

        #build sigma randomly if required
        if sigmaRand:
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
        print 'finished sigma at %f' % (time() - t)

        #build V
        edges = [] ###how many? *m^2*n^2
        for (i, j) in L.edges(labels=False):
            for hyp in phi[(i, (i, j))]:
                for x in hyp:
                    newEdges = [((i, x), (j, y)) for y in sigma[(i, j, tuple(hyp))]]
                    edges.extend(newEdges)
        V = Graph(edges)
        print 'finished V at %f' % (time() - t)

        #build D_i, F_i and A_i
        D_i = [0]*m
        for x in range(m):
            D_i[x] = sum([E[tuple(phi[x, line][0])] for line in L_i[x]])
        F_i = [1 - D_i[x] - ones/v for x in range(m)] #as the sum of (1/v)*J_\Omega_i, D_i, F_i is identity
        A_i = [0]*m
        for x in range(m):
            A_i[x] = ((v-k)/v)*ones - k*F_i[x] #we know A_i = k''*(1/v)*J_\Omega_i + r''*D_i + s''*F_i, and (k'', s'', r'') = (v - k, 0, -k)
        print 'finished D, F and A at %f' % (time() - t)

        #add the edges of the graph of B to V, and add the graph to the list
        for i in range(m):
            V.add_edges([((i, x), (i, y)) for x in range(n^d) for y in range(n^d) if not A_i[i][(x, y)]]) ### just make one list?
        Glist[I] = V
        print 'finished cycle at %f' % (time() - t)

    #return the answer
    if List is False:
        return Glist[0]
    else:
        return Glist