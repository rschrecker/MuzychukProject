def MuzychukS1_1Graph(n, d, phi='random', sigma='random'):
    from sage.combinat.designs.block_design import ProjectiveGeometryDesign
    from random import randrange
    from time import time
    t = time()

    #build L, L_i and the affine design
    ProjDesign = \
    ProjectiveGeometryDesign(d, 1, GF(n, 'a'), point_coordinates=False)
    L = [tuple(line) for line in ProjDesign.blocks()]
    m = ProjDesign.num_points()
    L_i = [0]*m
    for x in range(m):
        L_i[x] = [line for line in L if x in line]
    Design = \
    ProjectiveGeometryDesign(d, d-1, GF(n, 'a'), point_coordinates=False)
    projBlocks = Design.blocks()
    atInf = projBlocks[-1]
    Blocks = [[x for x in block if x not in atInf] for block in projBlocks[:-1]]
    print 'finished preamble at %f' % (time() - t)

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
    print 'finished ParClasses at %f' % (time() - t)

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
        assert type(phi) == dict, \
        'phi must be a dictionary or \'fixed\': alternatively, \
        remove this argument and it will be generated randomly'
        assert set(phi.keys()) == set([(x, line) for x in range(m)
            for line in L_i[x]]), 'each phi_i must have domain L_i'
        for x in range(m):
            assert m - 2 == len(set([val for (key, val) \
                in phi.items() if key[0] == x])), \
            'each phi_i must be injective'
        for val in phi.values():
            assert val in range(m-1), \
            'codomain should be {0,..., (n^d - 1)/(n - 1) - 1}'
    for x in range(m):
        for line in L_i[x]:
            phi[(x, line)] = ParClasses[phi[(x, line)]]
    print 'finished phi at %f' % (time() - t)

    #handle sigma
    ### VERY slow, order of 10s ###
    if sigma == 'random':
        sigma = {}
        for x in range(m):
            for line in L_i[x]:
                for j in line[1:]:
                    temp = phi[(j, line)][:]
                    for hyp in phi[(line[0], line)]:
                        rand = randrange(0, len(temp))
                        sigma[(line[0], j, tuple(hyp))] = temp[rand]
                        sigma[(j, line[0], tuple(temp[rand]))] = hyp
                        del temp[rand]
                    for i in line[1:]:
                        if i != j:
                            new_items = [((i, j, key[2]), \
                                sigma[(line[0], j, tuple(val))]) \
                                for (key, val) in sigma.items() \
                                if key[:2] == (i, line[0])]
                            sigma.update(new_items)
    elif sigma == 'fixed':
        sigma = {}
        for x in range(m):
            for line in L_i[x]:
                for j in line[1:]:
                    temp = phi[(j, line)][:]
                    for hyp in phi[(line[0], line)]:
                        val = temp.pop()
                        sigma[(line[0], j, tuple(hyp))] = val
                        sigma[(j, line[0], tuple(val))] = hyp
                        val +=1
                    for i in line[1:]:
                        if i != j:
                            new_items = [((i, j, key[2]), \
                                sigma[(line[0], j, tuple(val))]) \
                                for (key, val) in sigma.items() \
                                if key[:2] == (i, line[0])]
                            sigma.update(new_items)
    else:
        assert type(sigma) == dict, \
        'sigma must be a dictionary or \'fixed\': alternatively, \
        remove this argument and it will be generated randomly'
        correctKeys =      [(line[0], line[1], n) for line in L.edges() for n in range(len(ParClasses)) if ParClasses[n] in phi[(line[1], line)]]
        correctKeys.extend([(line[1], line[0], n) for line in L.edges() for n in range(len(ParClasses)) if ParClasses[n] in phi[(line[0], line)]])
        ### ^no ###
        assert set(sigma.keys()) == set(correct(Keys)), \
        'the keys in sigma must be \
        {(i, j, n) | i, j in line in L, and n in phi[i, line]}'
        for key in sigma.keys():
            assert key == sigma[(key[1], key[0], sigma[key])], \
            'sigma_ij must be (sigma_ji)^(-1)' ### no ###
        sigma = dict(((i, j, tuple(ParClasses[x]), ParClasses[y]) \
                              for ((i, j, x), y) in sigma.items()))
    print 'finished sigma at %f' % (time() - t)

    #build V
    edges = []
    for line in L:
        for index in range(n+1):
            i = line[index]
            for j in line[:index]:
                for hyp in phi[(i, line)]:
                    for x in hyp:
                        newEdges = [((i, x), (j, y)) for y \
                                    in sigma[(i, j, tuple(hyp))]]
                        edges.extend(newEdges)
    V = Graph(edges)
    print 'finished V at %f' % (time() - t)

    #add the edges of the graph O
    for i in range(m):
        newEdges = [((i, x), (i, y)) for x in range(n^d) for y in range(n^d)\
                    if x != y]
        #V.add_edges(newEdges)

    return V
