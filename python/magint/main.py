import argparse
from MagInt.run_dos import run_dos
from MagInt.run_shift_mu import run_shift_mu
from MagInt.run_magint import run_magint
from MagInt.run_shell import run_shell
import sys


def main():
    parser = argparse.ArgumentParser(description="Daje")
    parser.add_argument('command', choices=['dos', 'shift_mu', 'magint', 'shell'], help='The class to use')
    parser.add_argument('mu_shift', type=float, nargs='?', default=0.0, help='Shift in mu value')

    args = parser.parse_args()

    if args.command == 'dos':
        run_dos()
    elif args.command == 'shift_mu':
        run_shift_mu(args.mu_shift)
    elif args.command == 'magint':
        run_magint()
    elif args.command == 'shell':
        run_shell()


if __name__ == "__main__":
    main()
