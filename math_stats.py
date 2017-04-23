import math


def logAdd(a, b):
    """
    return log(exp(a) + exp(b))
    Source: https://github.com/drtconway/pykmer
    """
    x = max(a, b)
    y = min(a, b)
    w = y - x
    return x+math.log1p(math.exp(w))


def logSum(xs):
    """
    return log(sum([exp(x) for x in xs]))
    Source: https://github.com/drtconway/pykmer
    """
    assert len(xs) > 0
    y = xs[0]
    for x in xs[1:]:
        y = logAdd(y, x)
    return y