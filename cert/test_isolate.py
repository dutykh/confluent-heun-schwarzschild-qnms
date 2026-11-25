"""Meaningful exact checks for the root-counting layer (no physical evaluator)."""
import unittest
from fractions import Fraction as F
from cert.isolate import (RealBox,ComplexBox,sqrt_bounds,exp_upper,log_upper_ge_one,
                         make_certificate,schwarzschild_outer_bound,complex_box_from_iv)


def point(x,y=0,eps=F(0)):
    x,y=F(x),F(y)
    return ComplexBox(RealBox(x-eps,x+eps),RealBox(y-eps,y+eps))


class ExactIsolationTests(unittest.TestCase):
    def test_square_roots_are_enclosed_exactly(self):
        for x in [F(0),F(1),F(2),F(7,19),F(10**30,3),F(1,10**90)]:
            lo,hi=sqrt_bounds(x,256)
            self.assertLessEqual(lo*lo,x)
            self.assertLessEqual(x,hi*hi)
            self.assertLessEqual(hi-lo,F(1,2**256))

    def test_elementary_upper_bounds_sanity(self):
        self.assertEqual(exp_upper(F(0)),1)
        self.assertGreater(exp_upper(F(1)),F(8,3))
        self.assertLess(exp_upper(F(1)),F(11,4))
        self.assertGreater(log_upper_ge_one(F(2)),F(2,3))
        self.assertLess(log_upper_ge_one(F(2)),F(7,10))

    def test_quadratic_with_one_root_in_disk_passes(self):
        h=F(1,1000)
        fun=lambda z:z+2*z*z
        result=make_certificate(('0','0'),'0.01','0.1','0.12',h,
             point(fun(0),eps=F(1,10**30)),point(fun(h),eps=F(1,10**30)),
             point(fun(-h),eps=F(1,10**30)))
        self.assertTrue(result['proved'])
        # The two exact zeros are 0 and -1/2, so this checks the intended count.
        self.assertGreater(F(result['inequality']['strict_gap']['exact']),0)

    def test_offcenter_simple_root_passes(self):
        h=F(1,1000)
        fun=lambda z:z-F(1,1000)
        result=make_certificate((0,0),'0.01','0.1','0.101',h,
             point(fun(0)),point(fun(h)),point(fun(-h)))
        self.assertTrue(result['proved'])

    def test_two_roots_does_not_get_false_one_root_certificate(self):
        h=F(1,1000)
        fun=lambda z:z*z-F(1,10**6)
        result=make_certificate((0,0),'0.01','0.1','0.010001',h,
             point(fun(0)),point(fun(h)),point(fun(-h)))
        self.assertFalse(result['proved'])

    def test_wide_samples_fail_safely(self):
        result=make_certificate((0,0),'0.01','0.1','1','0.001',
             point(0,eps=F(1,10)),point(F(1,1000),eps=F(1,10)),
             point(-F(1,1000),eps=F(1,10)))
        self.assertFalse(result['proved'])

    def test_physical_outer_bound_and_phase_rejection(self):
        fundamental=('0.2209098781608393717509230','-0.2097914341737619175634781')
        result=schwarzschild_outer_bound(fundamental,'0.01',0,'1.5')
        self.assertEqual(result['outer_sup'],'62')
        with self.assertRaises(ValueError):
            schwarzschild_outer_bound(fundamental,'0.01',0,'1.5',pi_upper=3)
        with self.assertRaises(ValueError):
            schwarzschild_outer_bound(('0.1514838710703515','-1.202157180071607'),'0.01',0,'1.5')

    def test_iv_conversion_encloses_exact_decimal_input(self):
        try:
            from mpmath import iv
        except ImportError:
            self.skipTest('mpmath unavailable; checker itself has no mpmath dependency')
        old=iv.dps
        try:
            iv.dps=50
            box=complex_box_from_iv(iv.mpc('0.1','-0.2'))
            self.assertLessEqual(box.re.lo,F(1,10))
            self.assertGreaterEqual(box.re.hi,F(1,10))
            self.assertLessEqual(box.im.lo,-F(1,5))
            self.assertGreaterEqual(box.im.hi,-F(1,5))
        finally:
            iv.dps=old

if __name__=='__main__': unittest.main()
