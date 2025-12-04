from random import choice
from string import ascii_letters, digits


_CHARS = ascii_letters + digits


def puid(length: int = 10):
    "Probably Unique Identifier"
    return "".join(choice(_CHARS) for _ in range(length))
