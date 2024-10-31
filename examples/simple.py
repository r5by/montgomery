from mont.montgomery import Montgomery

if __name__ == '__main__':

    p = 61
    # M5 = Montgomery.factory(mod=p, mul_opt='real5').build(m=5)
    M8 = Montgomery.factory(mod=p, mul_opt='real8').build(m=64, w=256)
    # M0 = Montgomery.factory(mod=p, mul_opt='real0').build(R=M8.R)

    # A fast real-8 mont domain
    M0 = Montgomery.factory(mod=p, mul_opt='real8').build(m=64, w=256)
    M0 = M0.config(mul_opt='real0')

    x_, y_ = 8, 27  # x', y' as elements in prime field F_{31}

    # enter the Montgomery domain over F_{31}
    # x, y = M5(7), M5(21)
    # x0, y0 = M0(7), M0(21)
    x, y = M8(x_), M8(y_)
    x0, y0 = M0(x_), M0(y_)

    R = M8.R
    assert M0.ONE == R % p

    R_inv = M0._R_inv if M0._R_inv > 0 else p + M0._R_inv
    _one = (R * R_inv) % p
    assert _one == 1

    r = M8.r

    a = x + y
    a0 = x0 + y0
    a_ = (x_ + y_) * R % p
    assert a_ == a
    assert a == a0

    b = x - y
    b0 = x0 - y0
    b_ = (x_ - y_) * R % p
    assert  b == b_
    assert b == b0

    c = x * y
    c0 = x0 * y0
    c_ = (x_ * y_) * R % p
    assert c == c_
    assert c0 == c

    d = 1 / x
    d0 = 1 / x0
    assert d == d0
    assert M8.ONE == R % p
    assert d * x == M8.ONE


    e = x / y
    e0 = x0 / y0
    f = x ** int(e)
    g = y ** 65537
    assert e == e0

    # exit the Montomery domain
    _x, _y = int(x), int(y)
    _x0, _y0 = int(x0), int(y0)
    assert _x == x_
    assert _x == _x0 and _y == _y0
    assert _y == y_
    print(f'succeeded!')