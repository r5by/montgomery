import random
from mont.montgomery import Montgomery

"""
    Sample code to show-case usage of py-mont package on SM2
"""
if __name__ == '__main__':

    # Modulo for SM2
    p = 0xfffffffeffffffffffffffffffffffffffffffff00000000ffffffffffffffff
    M8 = Montgomery.factory(mod=p, mul_opt='real8').build(m=64, w=256)
    # M0 = Montgomery.factory(mod=p, mul_opt='real0').build(R=M8.R)

    # A fast real-8 mont domain (used for SW quick verification)
    M0 = Montgomery.factory(mod=p, mul_opt='real8').build(m=64, w=256)
    M0 = M0.config(mul_opt='real0')

    # x_, y_ = 8, 27  # x', y' as elements in prime field F_p
    x_ = random.randint(1, p - 1)
    y_ = random.randint(1, p - 1)

    # enter the Montgomery domain over F_{31}
    x, y = M8(x_), M8(y_)
    x0, y0 = M0(x_), M0(y_)

    R = M8.R
    assert R == M0.R

    # 1) multiplication
    c = x * y
    c0 = x0 * y0
    c_ = (x_ * y_) * R % p
    assert c == c_
    assert c0 == c

    # 2) inversion
    d = 1 / x
    d0 = 1 / x0
    assert d == d0

    # 3) division
    e = x / y
    e0 = x0 / y0
    assert e == e0

    # exit the Montomery domain
    _x, _y = int(x), int(y)
    _x0, _y0 = int(x0), int(y0)
    assert _x == x_
    assert _y == y_
    assert _x == _x0 and _y == _y0

    print(f'succeeded!')