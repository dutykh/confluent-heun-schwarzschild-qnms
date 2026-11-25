"""Outward-rounded MPFI elementary functions for mpmath.iv intervals.

All transfers use exact binary dyadic endpoint tuples. No binary-double or
rounded decimal conversion is used at the Python/C boundary. MPFI supplies
real transcendental enclosures; complex rectangles are formed with mpmath.iv
rational arithmetic. Requires system MPFI, MPFR and their normal GMP ABI.

API: cexp, clog, rexp, rlog, rsqrt, cabs, carg, sin, cos, pi.
Aliases: exp/log/sqrt dispatch on real versus complex input.
The logarithm is principal and rejects rectangles meeting its branch cut by
 default. An explicit allow_branch_cut=True returns a potentially wide hull,
 which must not be mistaken for a holomorphic branch on that rectangle.
"""
from __future__ import annotations
import ctypes as C
import ctypes.util
from contextlib import contextmanager
from fractions import Fraction
import mpmath as mp
from mpmath import libmp


class MPFR(C.Structure):
    _fields_ = [('prec', C.c_long), ('sign', C.c_int), ('exp', C.c_long),
                ('limbs', C.POINTER(C.c_ulong))]


class MPFI(C.Structure):
    _fields_ = [('left', MPFR), ('right', MPFR)]


_mpfr = C.CDLL(ctypes.util.find_library('mpfr') or 'libmpfr.so.6')
_mpfi = C.CDLL(ctypes.util.find_library('mpfi') or 'libmpfi.so.0')
PFR, PFI = C.POINTER(MPFR), C.POINTER(MPFI)


def _bind(lib, name, args, result=C.c_int):
    fn = getattr(lib, name)
    fn.argtypes, fn.restype = args, result
    return fn


_init = _bind(_mpfi, 'mpfi_init2', [PFI, C.c_long], None)
_clear = _bind(_mpfi, 'mpfi_clear', [PFI], None)
_setstr = _bind(_mpfr, 'mpfr_set_str', [PFR, C.c_char_p, C.c_int, C.c_int])
_mul2 = _bind(_mpfr, 'mpfr_mul_2si', [PFR, PFR, C.c_long, C.c_int])
_setzero = _bind(_mpfr, 'mpfr_set_zero', [PFR, C.c_int], None)
_setinf = _bind(_mpfr, 'mpfr_set_inf', [PFR, C.c_int], None)
_nan = _bind(_mpfr, 'mpfr_nan_p', [PFR])
_inf = _bind(_mpfr, 'mpfr_inf_p', [PFR])
_zero = _bind(_mpfr, 'mpfr_zero_p', [PFR])
_getstr = _bind(_mpfr, 'mpfr_get_str',
               [C.c_void_p, C.POINTER(C.c_long), C.c_int, C.c_size_t, PFR, C.c_int],
               C.c_void_p)
_freestr = _bind(_mpfr, 'mpfr_free_str', [C.c_void_p], None)
_mpfr_version = _bind(_mpfr, 'mpfr_get_version', [], C.c_char_p)
_mpfi_version = _bind(_mpfi, 'mpfi_get_version', [], C.c_char_p)
_unary = {name: _bind(_mpfi, 'mpfi_'+name, [PFI, PFI])
          for name in ('exp', 'log', 'sqrt', 'sin', 'cos', 'abs', 'sqr')}
_hypot = _bind(_mpfi, 'mpfi_hypot', [PFI, PFI, PFI])
_atan2 = _bind(_mpfi, 'mpfi_atan2', [PFI, PFI, PFI])
_pi = _bind(_mpfi, 'mpfi_const_pi', [PFI])
RNDN = 0


def raw_to_fraction(raw):
    """Exact finite mpmath mpf tuple -> fractions.Fraction."""
    sign, man, exponent, bc = raw
    if bc < 0:
        raise ValueError('nonfinite endpoint has no rational representation')
    value = Fraction((-1 if sign else 1)*int(man))
    return value*(2**exponent) if exponent >= 0 else value/Fraction(2**(-exponent))


def endpoints(x):
    """Return the two exact Fraction endpoints of a finite real IV interval."""
    x = _real(x)
    return tuple(raw_to_fraction(v) for v in x._mpi_)


def rational(p, q=1):
    """Outward interval enclosing an exact integer quotient p/q."""
    if not q:
        raise ZeroDivisionError
    return mp.iv.mpf(int(p))/mp.iv.mpf(int(q))


def _real(x):
    if hasattr(x, '_mpci_'):
        if x._mpci_[1] != (libmp.fzero, libmp.fzero):
            raise TypeError('expected a real interval')
        return mp.iv.make_mpf(x._mpci_[0])
    if isinstance(x, Fraction):
        return rational(x.numerator, x.denominator)
    return x if hasattr(x, '_mpi_') else mp.iv.mpf(x)


def _complex(z):
    return z if hasattr(z, '_mpci_') else mp.iv.mpc(_real(z), 0)


def _maxbits(*intervals):
    return max([int(mp.iv.prec), 2]+
               [int(raw[3]) for x in intervals for raw in x._mpi_ if raw[3] >= 0])


def _set_raw(dst, raw):
    sign, man, exponent, bc = raw
    ptr = C.byref(dst)
    if bc < 0:
        if raw == libmp.finf or raw == libmp.fninf:
            _setinf(ptr, -1 if sign else 1)
            return
        raise ValueError('NaN input interval')
    if not man:
        _setzero(ptr, 1)
        return
    # Destination precision is at least the significand bit count: both
    # operations are exact. set_str reports parse success; mul_2si returns
    # a ternary inexact indicator and must return zero here.
    binary = ('-' if sign else '')+format(int(man), 'b')
    if _setstr(ptr, binary.encode('ascii'), 2, RNDN) != 0:
        raise ArithmeticError('MPFR binary integer import failed')
    if _mul2(ptr, ptr, int(exponent), RNDN) != 0:
        raise ArithmeticError('MPFR dyadic scaling unexpectedly rounded')


def _get_raw(src):
    ptr = C.byref(src)
    if _nan(ptr):
        raise ArithmeticError('MPFI returned NaN')
    if _inf(ptr):
        return libmp.fninf if src.sign < 0 else libmp.finf
    if _zero(ptr):
        return libmp.fzero
    exponent = C.c_long()
    # In base 2, requesting exactly the MPFR precision gives the exact finite
    # binary significand. This is a serialization, not an approximation.
    address = _getstr(None, C.byref(exponent), 2, int(src.prec), ptr, RNDN)
    if not address:
        raise MemoryError('mpfr_get_str failed')
    try:
        digits = C.string_at(address).decode('ascii')
    finally:
        _freestr(address)
    negative = digits.startswith('-')
    digits = digits.lstrip('-')
    man = int(digits, 2)
    return libmp.from_man_exp(-man if negative else man, int(exponent.value)-len(digits))


@contextmanager
def _interval(precision, source=None):
    value = MPFI()
    _init(C.byref(value), int(precision))
    try:
        if source is not None:
            _set_raw(value.left, source._mpi_[0])
            _set_raw(value.right, source._mpi_[1])
        yield value
    finally:
        _clear(C.byref(value))


def _export(value):
    lo, hi = _get_raw(value.left), _get_raw(value.right)
    if libmp.mpf_cmp(lo, hi) > 0:
        raise ArithmeticError('MPFI returned reversed interval endpoints')
    return mp.iv.make_mpf((lo, hi))


def _run_unary(name, x):
    x = _real(x)
    precision = _maxbits(x)
    with _interval(precision, x) as a, _interval(precision) as result:
        _unary[name](C.byref(result), C.byref(a))
        return _export(result)


def _run_binary(fn, x, y):
    x, y = _real(x), _real(y)
    precision = _maxbits(x, y)
    with _interval(precision, x) as a, _interval(precision, y) as b, _interval(precision) as result:
        fn(C.byref(result), C.byref(a), C.byref(b))
        return _export(result)


def rexp(x):
    return _run_unary('exp', x)


def rlog(x):
    x = _real(x)
    if libmp.mpf_cmp(x._mpi_[0], libmp.fzero) <= 0:
        raise ValueError('rlog requires a strictly positive real interval')
    return _run_unary('log', x)


def rsqrt(x):
    x = _real(x)
    if libmp.mpf_cmp(x._mpi_[0], libmp.fzero) < 0:
        raise ValueError('rsqrt requires a nonnegative real interval')
    return _run_unary('sqrt', x)


def sin(x):
    return _run_unary('sin', x)


def cos(x):
    return _run_unary('cos', x)


def rabs(x):
    return _run_unary('abs', x)


def pi():
    with _interval(mp.iv.prec) as result:
        _pi(C.byref(result))
        return _export(result)


def cabs(z):
    """MPFI hypot encloses the modulus on the entire input rectangle."""
    z = _complex(z)
    return _run_binary(_hypot, z.real, z.imag)


def _contains_zero(x):
    lo, hi = x._mpi_
    return libmp.mpf_cmp(lo, libmp.fzero) <= 0 and libmp.mpf_cmp(hi, libmp.fzero) >= 0


def carg(z):
    """Principal argument hull; rejects rectangles containing the origin."""
    z = _complex(z)
    if _contains_zero(z.real) and _contains_zero(z.imag):
        raise ValueError('argument is undefined on a rectangle containing zero')
    return _run_binary(_atan2, z.imag, z.real)


def cexp(z):
    z = _complex(z)
    magnitude = rexp(z.real)
    return mp.iv.mpc(magnitude*cos(z.imag), magnitude*sin(z.imag))


def clog(z, *, allow_branch_cut=False):
    z = _complex(z)
    if not allow_branch_cut and _contains_zero(z.imag) and libmp.mpf_cmp(z.real._mpi_[0], libmp.fzero) <= 0:
        raise ValueError('principal complex log rectangle meets the branch cut')
    return mp.iv.mpc(rlog(cabs(z)), carg(z))


def _nonnegative_hull(x):
    """Intersect a known nonnegative expression's enclosure with [0,+inf)."""
    lo, hi = x._mpi_
    if libmp.mpf_cmp(hi, libmp.fzero) < 0:
        raise ArithmeticError('negative enclosure for a known nonnegative expression')
    return mp.iv.make_mpf((max(lo, libmp.fzero, key=lambda t: raw_to_fraction(t)), hi))


def csqrt(z):
    z = _complex(z)
    modulus = cabs(z)
    real = rsqrt(_nonnegative_hull((modulus+z.real)/2))
    imag_abs = rsqrt(_nonnegative_hull((modulus-z.real)/2))
    ilo, ihi = z.imag._mpi_
    if libmp.mpf_cmp(ilo, libmp.fzero) >= 0:
        imaginary = imag_abs
    elif (libmp.mpf_cmp(ihi, libmp.fzero) < 0 or
          (libmp.mpf_cmp(ihi, libmp.fzero) == 0 and
           libmp.mpf_cmp(z.real._mpi_[0], libmp.fzero) >= 0)):
        imaginary = -imag_abs
    else:
        imaginary = mp.iv.make_mpf((libmp.mpf_neg(imag_abs._mpi_[1]), imag_abs._mpi_[1]))
    return mp.iv.mpc(real, imaginary)


def exp(x):
    return cexp(x) if hasattr(x, '_mpci_') else rexp(x)


def log(x):
    return clog(x) if hasattr(x, '_mpci_') else rlog(x)


def sqrt(x):
    return csqrt(x) if hasattr(x, '_mpci_') else rsqrt(x)


# Clear descriptive aliases for callers not using the compact API.
real_exp, real_log, real_sqrt = rexp, rlog, rsqrt
abs_complex, phase = cabs, carg


def backend_info():
    return {'mpfi':_mpfi_version().decode(), 'mpfr':_mpfr_version().decode(),
            'interval_precision_bits':int(mp.iv.prec),
            'transfer':'exact binary dyadic endpoints',
            'mpfr_struct_bytes':C.sizeof(MPFR), 'mpfi_struct_bytes':C.sizeof(MPFI)}


def install_iv_patch():
    """Explicit opt-in patch for iv.exp/log/sqrt/sin/cos; returns originals.

    Use cabs rather than built-in abs for complex intervals, and pi() rather
    than a cached constant when changing precision. Rational iv operations
    and integer powers remain mpmath's directed-rounding implementations.
    """
    originals = {name:getattr(mp.iv,name) for name in ('exp','log','sqrt','sin','cos')}
    for name in originals:
        setattr(mp.iv,name,globals()[name])
    return originals


def self_test():
    # Exact round-trip including signs, widely separated exponents, and
    # endpoints whose bit counts exceed the current IV context precision.
    for bits in (53,113,256):
        mp.iv.prec=bits
        raws=[libmp.from_man_exp(0,0),libmp.from_man_exp(7,-11),
              libmp.from_man_exp(-123456789,-100),libmp.from_man_exp((1<<400)+17,-300),
              libmp.from_man_exp(-3,123)]
        for raw in raws:
            x=mp.iv.make_mpf((raw,raw))
            with _interval(_maxbits(x),x) as value:
                assert _get_raw(value.left)==raw and _get_raw(value.right)==raw
        assert rexp(0)._mpi_==(libmp.fone,libmp.fone)
        assert rlog(1)._mpi_==(libmp.fzero,libmp.fzero)
        assert rsqrt(4)._mpi_==(libmp.from_int(2),libmp.from_int(2))
        assert sin(0)._mpi_==(libmp.fzero,libmp.fzero)
        assert cos(0)._mpi_==(libmp.fone,libmp.fone)
        # Inclusion tests below are secondary diagnostics; exact transfer and
        # MPFI's directed-rounding contract provide the enclosure guarantee.
        with mp.workprec(bits+100):
            for x in ('0.1','1','-2.75','100'):
                y=rexp(mp.iv.mpf(x)); exact=mp.exp(mp.mpf(x))
                assert raw_to_fraction(y._mpi_[0]) <= raw_to_fraction(exact._mpf_) <= raw_to_fraction(y._mpi_[1])
            for real,imag in [('0.2','0.7'),('-1.2','-2.3'),('3.5','0.01')]:
                z=mp.iv.mpc(real,imag); result=cexp(z); exact=mp.exp(mp.mpc(real,imag))
                for enclosure,component in [(result.real,exact.real),(result.imag,exact.imag)]:
                    q=raw_to_fraction(component._mpf_)
                    assert raw_to_fraction(enclosure._mpi_[0]) <= q <= raw_to_fraction(enclosure._mpi_[1])
        broad=sin(mp.iv.mpf([-10,10]))
        assert endpoints(broad)==(Fraction(-1),Fraction(1))
        assert endpoints(cabs(mp.iv.mpc([3,3],[4,4])))==(Fraction(5),Fraction(5))
        # Branch protection and conservative phase hull at the negative cut.
        try: clog(mp.iv.mpc([-2,-1],[-1,1]))
        except ValueError: pass
        else: raise AssertionError('branch-cut check failed')
        angle=carg(mp.iv.mpc([-2,-1],[-1,1]))
        p=pi()
        assert libmp.mpf_cmp(angle._mpi_[0],libmp.mpf_neg(p._mpi_[0]))<=0
        assert libmp.mpf_cmp(angle._mpi_[1],p._mpi_[0])>=0
    return {'status':'PASS', **backend_info()}


if __name__=='__main__':
    import json
    print(json.dumps(self_test(),indent=2))
