"""Check the Schultze phase benchmark with exact SymPy algebra.

Run with Python and SymPy. All parameters in the coefficient regression are
rational complex numbers; the identity of the normalized columns is first
checked with symbolic parameters. No large-order numerical limit is used.
"""
import sympy as sp


def require_zero(expression):
    if sp.cancel(sp.expand(expression)) != 0:
        raise ValueError("Symbolic identity failed")


def verify_solution_and_normalizations():
    alpha, beta, gamma = sp.symbols("alpha beta gamma", nonzero=True)
    z, x, t = sp.symbols("z x t")

    A = alpha * (beta + 2)
    B = (beta + 1) * (gamma - alpha)
    rho = -A / alpha
    E = -rho - beta - gamma - 2

    # H = exp(-alpha*z)*(1-z)^(-gamma), on the branch analytic at zero.
    # Check L[H]/H via its rational logarithmic derivative.
    logH_z = -alpha - gamma / (z - 1)
    Hzz_over_H = sp.diff(logH_z, z) + logH_z**2
    first = alpha + (beta + 1)/z + (gamma + 1)/(z - 1)
    potential = (A*z+B)/(z*(z-1))
    require_zero(Hzz_over_H + first*logH_z + potential)
    require_zero(E + gamma)

    # At x=1-z>1, H(1-x)=exp(alpha*(x-1))*x^(-gamma).
    # Reflection flips the derivative sign, interchanges the simple poles,
    # and sends (A,B) to (-A,A+B).
    require_zero(-logH_z.subs(z, 1-x) - (alpha - gamma/x))
    require_zero(-first.subs(z, 1-x) -
                 (-alpha + (gamma+1)/x + (beta+1)/(x-1)))
    require_zero(potential.subs(z, 1-x) - (-A*x+A+B)/(x*(x-1)))

    # x^(-rho) H(1-x), x=1/(1-t), has analytic value 1 at t=0.
    # Its log derivative fixes the generating function without symbolic
    # branch simplifications: (1-t)^(rho+gamma)*exp(alpha*t/(1-t)).
    x_of_t = 1/(1-t)
    logY_t = (alpha-(rho+gamma)/x_of_t)*sp.diff(x_of_t, t)
    require_zero(logY_t - (-(rho+gamma)/(1-t) + alpha/(1-t)**2))

    # The normalized exponential in original z satisfies
    # exp(-alpha*z) z^E / [exp(alpha*x) x^E] -> exp(-alpha+i*pi*E)
    # on the upper negative z-ray. The exact solution on reflected x>1
    # has exponential coefficient exp(-alpha). Subtract the exponents
    # rather than forcing invalid power simplifications across a cut.
    reflected_log_coefficient = -alpha
    upper_chart_log_factor = -alpha + sp.I*sp.pi*E
    require_zero(reflected_log_coefficient - upper_chart_log_factor
                 - sp.I*sp.pi*gamma)
    return alpha, beta, gamma, rho


def verify_coefficients():
    alpha, beta, gamma, rho = verify_solution_and_normalizations()
    case = {alpha:(4+3*sp.I)/5, beta:-sp.Rational(1,5)+3*sp.I/20,
            gamma:sp.Rational(2,5)-sp.I/4}
    mu = rho.subs(case)
    g, b, c = (case[q] for q in (alpha, beta, gamma))
    A, B = g*(b+2), (b+1)*(c-g)
    a0, a1 = -c, -b
    b0 = -(A+B)+(1+c)*(1+b+g)/2
    b1 = B+(1+b)*(g-1-c)/2

    # The exact generating function gives L_k^(d)(-g).
    # Compare with Eq. (6.41) expressed in unit Taylor coefficients.
    d = -mu-c-1
    require_zero((b0+b1)/g+(d+1)/2)
    coeffs = {0:sp.Integer(1), -1:sp.Integer(0), -2:sp.Integer(0)}
    for k in range(8):
        phi1 = (3*k**2 + k*(-2*(a1+mu)+a0+g)
                - b1-(1-a1)*(mu+(1-a0)/2-g/2))
        phi2 = (-3*(k-1)**2+(k-1)*(4*mu-2*a0+a1-g)
                +(1-a0)*(1-a1+g)/2-mu**2+mu*(1+a0+g-a1)-b0)
        phi3 = (k-2)**2+(k-2)*(a0-2*mu)+mu**2-mu*a0
        coeffs[k+1] = sp.cancel((phi1*coeffs[k]+phi2*coeffs[k-1]
                                 +phi3*coeffs[k-2])/((k+1)*(k+1-a1)))
        require_zero(coeffs[k+1]-sp.assoc_laguerre(k+1, d, -g))


if __name__ == "__main__":
    verify_coefficients()
    print("CHE, reflection, lifted phase, generating function and eight coefficients: PASS")
