"""Run the archived scalar certificate from the unpacked project root.

--check-saved rechecks exact-rational inequalities using saved amplitude boxes.
Without it, all horizon, propagation, quadrature and infinity enclosures are
regenerated before the Rouché checks. Use a fresh copy to retain original JSON.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from .campaign import MODES, run
from .check import check
from .compare import run as compare


def evaluate(task):
    run(*task)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-saved',action='store_true')
    parser.add_argument('--workers',type=int,default=1)
    args=parser.parse_args()
    if args.workers < 1: parser.error('workers must be positive')
    if not args.check_saved:
        samples=[(mode,pos,None) for mode in sorted(MODES) for pos in ('center','plus','minus')]
        samples.extend(('00','center',r) for r in (20,30,40))
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(evaluate,samples))
    for mode in sorted(MODES):
        check(mode, check_saved=args.check_saved)
    compare(check_saved=args.check_saved)


if __name__=='__main__': main()
