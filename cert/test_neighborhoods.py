"""Exact regression tests for normalization and complex interpolation bounds."""
import unittest
from fractions import Fraction as F

from .isolate import ComplexBox, RealBox
from .neighborhoods import D, NODES, S, coefficient_boxes, interpolation_errors


def point(real, imag):
    return ComplexBox(RealBox(real, real), RealBox(imag, imag))


def eval_poly(coeffs, t):
    result = F(0)
    for coefficient in reversed(coeffs):
        result = result*t+coefficient
    return result


class NeighborhoodsTests(unittest.TestCase):
    def test_full_degree_polynomial(self):
        real = [F(k+1, k+7) for k in range(9)]
        imag = [F((-1)**k, 2*k+1) for k in range(9)]
        samples = {t:point(eval_poly(real,t),eval_poly(imag,t)) for t in NODES}
        co = coefficient_boxes(samples)
        for k, box in enumerate(co):
            self.assertEqual(box, point(real[k],imag[k]))
        # B_1 is the t derivative: dividing by S is essential in physical frequency.
        self.assertEqual(co[1].divide_positive(S),point(real[1]/S,imag[1]/S))

    def test_pole_outside_outer_circle(self):
        # G(t)=1/(110-t) has an explicit closed-form interpolation error
        # (G-P)(t)=Q(t)/((110-t)Q(110)) at every complex point inside the circle.
        pole = 110
        M = 1/(pole-D/S)
        samples = {t:point(1/(pole-t), 0) for t in NODES}
        co = coefficient_boxes(samples)
        coeffs = [x.re.mid for x in co]
        deriv = [k*x for k,x in enumerate(coeffs) if k]
        for radius in (F(1),F(1,10),F(3,2)):
            error, derivative_error = interpolation_errors(M,radius)
            for j in range(-15,16):
                t = radius * j/15
                self.assertLessEqual(abs(1/(pole-t)-eval_poly(coeffs,t)), error)
                self.assertLessEqual(abs(1/(pole-t)**2-eval_poly(deriv,t)), derivative_error)

    def test_zero_domain_rejected(self):
        for radius in (0,100,101):
            with self.assertRaises(ValueError):
                interpolation_errors(F(1),radius)


if __name__ == '__main__':
    unittest.main()
