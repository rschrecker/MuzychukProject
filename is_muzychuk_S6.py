def is_muzychuk_S6(v, k, l, mu):
    '''
    Test whether there is a Muzychuk S6 graph which is (v, k, l, mu)-strongly
    regular.

    INPUT:

        v, k, l, mu: integers

    OUTPUT:

        A tuple t such that t[0](*t[1:]) builds the required graph if it exists,
        and None otherwise.

    EXAMPLES:

        sage: from sage.graphs.strongly_regular_db import is_muzychuk_S6
        sage: t = is_muzychuk_S6(378, 116, 34, 36)
        sage: G = t[0](*t[1:]); G
        Muzychuk S6 graph with parameters (n=3, d=3): Graph on 378 vertices
        sage: G.is_strongly_regular(parameters=True)
        (378, 116, 34, 36)
        sage: t = is_muzychuk_S6(5, 5, 5, 5); t
    '''
    n_list = [n for n in range(l-1) if Integer(n).is_prime_power()]
    for n in n_list:
        d = 2
        while n^d * ((n^d-1)/(n-1)+1) <= v:
            if v == n^d * ((n^d-1)/(n-1)+1) and k == n^(d-1)*(n^d-1)/(n-1) - 1\
            and l == mu - 2 and mu == n^(d-1) * (n^(d-1)-1) / (n-1):
                from sage.graphs.generators.classical_geometries\
                import MuzychuckS6Graph
                return (MuzychukS6Graph, n, d)
            d += 1
