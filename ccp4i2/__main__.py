from argparse import ArgumentParser

from . import __version__
from .bin import browser


def main():
    parser = ArgumentParser()
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.parse_args()
    browser.main()


if __name__ == '__main__':
    main()
