"""
Allow running as: python -m filepicker
"""

from .cli import main
import sys

if __name__ == '__main__':
    sys.exit(main())
