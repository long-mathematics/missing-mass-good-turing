"""Check the phase-transition manuscript and write reproducible numerical examples."""

import json
from pathlib import Path

import mpmath as mp
import sympy as sp


def check_series():
    z, s = sp.symbols("z s")
    q = sp.sqrt(5)
    rho = (1 + q) / 2
    gamma = sp.sqrt(2 + 4*q/5)
    beta = (6 + 2*q) / 15
    chi = gamma * (sp.Rational(31, 180) + 2*q/45)
    delta = (605 + 251*q) / 3375
    phi = z * (z + 1) * sp.exp(-z)
    coefficients = [
        sp.simplify(sp.exp(rho) * sp.diff(phi, z, j).subs(z, rho) / sp.factorial(j))
        for j in range(6)
    ]
    # Coefficients of phi(rho+h) - (1-s^2)phi(rho), after scaling by exp(rho).
    h = sp.Poly(gamma*s + beta*s**2 + chi*s**3 + delta*s**4, s)
    polynomial = sp.Poly(coefficients[0]*s**2, s)
    for j in range(2, 6):
        polynomial += coefficients[j] * h**j
    for j in range(2, 6):
        assert sp.simplify(polynomial.nth(j)) == 0
    historical_chi = (40*sp.sqrt(2) + 31*sp.sqrt(10))*sp.sqrt(5 + 2*q)/900
    assert sp.simplify(chi - historical_chi) == 0
    assert sp.Rational(5, 4)*sp.Rational(2, 3)**11*sp.Rational(11, 4)**4 < 1


def bisect(fun, lo, hi, steps):
    flo = fun(lo)
    assert flo * fun(hi) < 0
    for _ in range(steps):
        mid = (lo + hi) / 2
        value = fun(mid)
        if value == 0:
            return mid
        if (value > 0) == (flo > 0):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def check_pair_certificates():
    """Compare the midpoint integral against implicit differentiation of the marginal."""
    count = 0
    with mp.workdps(75):
        for theta in (1, 2, 3, 7, 11):
            b = (1 + mp.sqrt(1 + 4*theta**2)) / 2
            poles = [b, 1-b] + [mp.mpf(0)]*(theta-1)

            def g(z):
                return z**(theta-1) * (z*z-z-theta**2) * mp.exp(-z)

            def log_slope(z):
                return mp.fsum(1/(z-a) for a in poles) - 1

            c = bisect(log_slope, b + mp.mpf("1e-20"), b + 2*theta + 4, 275)
            for fraction in ("0.01", "0.25", "0.75", "0.99"):
                x = b + mp.mpf(fraction)*(c-b)
                hi = c+1
                while g(hi) > g(x):
                    hi += hi-c
                y = bisect(lambda z: g(z)-g(x), c, hi, 275)
                slope = mp.diff(g, x) / mp.diff(g, y)
                curvature = (mp.diff(g, x, 2) - mp.diff(g, y, 2)*slope**2) / mp.diff(g, y)
                m, h = (x+y)/2, (y-x)/2

                def integral(integrand):
                    return mp.quad(lambda u: mp.fsum(integrand(u, m-a+h*u) for a in poles),
                                   [-1, 0, 1])/2

                weight = integral(lambda u, q: q**-2)
                m1 = -integral(lambda u, q: u/q**2)/weight
                m2 = 2*integral(lambda u, q: (m1+u)**2/q**3)/weight
                assert abs(integral(lambda u, q: 1/q) - 1) < mp.mpf("1e-60")
                assert 0 < m1 < 1 and curvature > 0
                assert abs(curvature - 2*m2/(1-m1)**3) < mp.mpf("1e-58") * max(1, curvature)
                count += 1
    return count


def phase_at(n):
    # Direct gap parametrization from the quadratic marginal-equality equation.
    # The precision also resolves the tiny difference between the two thresholds.
    with mp.workdps(52 + int(mp.sqrt(n))):
        rho = (1 + mp.sqrt(5))/2
        gamma = mp.sqrt(2 + 4*mp.sqrt(5)/5)
        beta = (6 + 2*mp.sqrt(5))/15

        def phi(z):
            return z*(1+z)*mp.exp(-z)

        def light(h):
            return (mp.mpf("0.5") + 2*h/mp.expm1(2*h)
                    + mp.sqrt(mp.mpf("1.25") + (h/mp.sinh(h))**2))

        def total(h):
            return n*light(h) + 2*h

        def excess(h):
            x = light(h)
            return (n-1)*phi(x) + phi(x+2*h) - n*phi(total(h)/n)

        def slope(h):
            return mp.diff(total, h)

        lo, hi = mp.mpf("0.0001"), mp.mpf(1)
        while slope(hi) <= 0:
            hi *= 2
        steps = int(3.5 * mp.mp.dps)
        h_fold = bisect(slope, lo, hi, steps)
        hi = 2*h_fold
        while excess(hi) <= 0:
            hi *= 2
        h_cross = bisect(excess, h_fold, hi, steps)
        fold, cross = total(h_fold), total(h_cross)
        assert n*rho < fold < cross < 3*n
        x = light(h_cross)
        y = x + 2*h_cross
        residual = abs(mp.diff(phi, x) - mp.diff(phi, y))
        assert residual < mp.mpf("1e-40") and abs(excess(h_cross)) < mp.mpf("1e-40")
        assert slope(h_cross) > 0 and excess(h_fold) < 0

        def surrogate_equation(t):
            return n*phi(t/n) - (n-1)*phi(rho)

        upper = (n+10)*rho + 4*gamma*mp.sqrt(n)
        surrogate = bisect(surrogate_equation, n*rho, upper, steps)
        displacement = surrogate - cross
        assert displacement > 0
        predicted = (gamma*mp.exp(-beta)/mp.sqrt(5)
                     * n*mp.sqrt(n)*mp.exp(-gamma*mp.sqrt(n)))
        ratio = displacement / predicted
        if n == 10000:
            assert abs(ratio-1) < mp.mpf("0.05")
        return dict(n=n, fold=mp.nstr(fold, 12), transition=mp.nstr(cross, 12),
                    displacement=mp.nstr(displacement, 12), ratio_to_leading_term=mp.nstr(ratio, 12))


if __name__ == "__main__":
    check_series()
    results = {
        "symbolic_checks": ["ordinary Taylor coefficients through degree five",
                            "equality of the two radical expressions for chi",
                            "rational upper bound for the higher-order light example"],
        "midpoint_quadrature_pairs": check_pair_certificates(),
        "occupancy_orders_tested": [0, 1, 2, 6, 10],
        "ordinary_examples": [phase_at(n) for n in (3, 4, 10, 100, 1000, 10000)],
    }
    output = json.dumps(results, indent=2) + "\n"
    (Path(__file__).resolve().parents[1] / "checks" / "good-turing-phase.json").write_text(output)
    print(output, end="")
