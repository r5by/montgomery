# simulator of critical circuit design of [2]

from common import csa, compress


def compressor_2m4to2(TS, TC,
                        X_j_IN, Yi,
                        q, Mj_IN,
                        FBS, FBC,
                        w, m):
    '''
        (2m+4)-to-2 compressor -- simulate the HW implementation (maybe this one works?)
        ref [2], fig. 4
    '''

    os, oc = TS, TC

    for i in range(m):
        Yib = (Yi >> i) & 1
        t1 = X_j_IN * Yib << i
        os, oc = csa(os, oc, t1, w + m)

        qb = (q >> i) & 1
        t2 = qb * Mj_IN << i
        os, oc = csa(os, oc, t2, w + m)

    os, oc = csa(os, oc, FBS, w + m + 1)
    os, oc = csa(os, oc, FBC, w + m + 2)

    return os, oc



