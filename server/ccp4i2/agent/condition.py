"""The ``when`` language of a judgement file: a condition on a job's results.

Deliberately tiny, and parsed here rather than handed to ``eval``, because
the files are data that anyone may edit:

    TFZ >= 8 and LLG > 60
    RFREE_START - RFREE >= 0.02 and RFREE - RWORK < 0.07
    outcome == "solved"
    not (RFREE > 0.35) or true

Names are result names; values are numbers, double-quoted strings, ``true``
and ``false``; ``+ - * /`` combine numbers. A comparison with a result that could not be read is unknown, and an
unknown condition does not hold, so a verdict list falls through to its
catch-all rather than claiming success on a missing number.
"""
import operator
import re

_TOKEN = re.compile(r"""
    \s*(?:
      (?P<num>\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)
    | "(?P<str>[^"]*)"
    | (?P<op><=|>=|==|!=|<|>)
    | (?P<arith>[-+*/])
    | (?P<paren>[()])
    | (?P<word>[A-Za-z_][A-Za-z0-9_]*)
    )""", re.X)

_COMPARE = {
    "<": operator.lt, "<=": operator.le, ">": operator.gt,
    ">=": operator.ge, "==": operator.eq, "!=": operator.ne,
}


_ARITH = {"+": operator.add, "-": operator.sub, "*": operator.mul, "/": operator.truediv}


class ConditionError(ValueError):
    pass


def _tokens(text):
    pos, out = 0, []
    text = text.rstrip()
    while pos < len(text):
        m = _TOKEN.match(text, pos)
        if not m or m.end() == pos:
            raise ConditionError(f"cannot read {text[pos:]!r} in {text!r}")
        kind = m.lastgroup
        value = m.group(kind)
        if kind == "num":
            value = float(value)
        out.append((kind, value))
        pos = m.end()
    return out


class _Missing:
    """A result that could not be read: every comparison with it is false."""

    def __repr__(self):
        return "missing"


MISSING = _Missing()


def parse(text):
    """Parse a condition into a tree; raise ConditionError if it is not one."""
    tokens = _tokens(str(text))
    pos = 0

    def peek():
        return tokens[pos] if pos < len(tokens) else (None, None)

    def take():
        nonlocal pos
        if pos >= len(tokens):
            raise ConditionError(f"{text!r} ends too soon")
        pos += 1
        return tokens[pos - 1]

    def expr():
        node = term()
        while peek() == ("word", "or"):
            take()
            node = ("or", node, term())
        return node

    def term():
        node = factor()
        while peek() == ("word", "and"):
            take()
            node = ("and", node, factor())
        return node

    def factor():
        if peek() == ("word", "not"):
            take()
            return ("not", factor())
        left = sum_()
        if peek()[0] == "op":
            op = take()[1]
            return ("cmp", op, left, sum_())
        return left

    def sum_():
        node = product()
        while peek() in (("arith", "+"), ("arith", "-")):
            op = take()[1]
            node = ("arith", op, node, product())
        return node

    def product():
        node = signed()
        while peek() in (("arith", "*"), ("arith", "/")):
            op = take()[1]
            node = ("arith", op, node, signed())
        return node

    def signed():
        if peek() == ("arith", "-"):
            take()
            return ("arith", "-", ("lit", 0.0), signed())
        return atom()

    def atom():
        kind, value = peek()
        if kind is None:
            raise ConditionError(f"{text!r} ends too soon")
        take()
        if (kind, value) == ("paren", "("):
            node = expr()
            if take() != ("paren", ")"):
                raise ConditionError(f"unbalanced parentheses in {text!r}")
            return node
        if kind in ("num", "str"):
            return ("lit", value)
        if kind == "word":
            if value in ("and", "or", "not"):
                raise ConditionError(f"{value!r} out of place in {text!r}")
            if value in ("true", "false"):
                return ("lit", value == "true")
            return ("name", value)
        raise ConditionError(f"{value!r} out of place in {text!r}")

    tree = expr()
    if pos != len(tokens):
        raise ConditionError(f"unexpected {tokens[pos][1]!r} in {text!r}")
    return tree


def names(tree):
    """The result names a condition refers to."""
    kind = tree[0]
    if kind == "name":
        return {tree[1]}
    if kind == "lit":
        return set()
    if kind == "not":
        return names(tree[1])
    if kind in ("cmp", "arith"):
        return names(tree[2]) | names(tree[3])
    return names(tree[1]) | names(tree[2])


def evaluate(tree, values):
    """Evaluate a parsed condition against ``values`` (name -> value).

    Three-valued: a comparison with a missing result is MISSING (unknown),
    ``not`` of unknown is unknown, and ``and``/``or`` follow Kleene's rules,
    so neither a condition nor its negation holds on a number not read.
    """
    kind = tree[0]
    if kind == "lit":
        return tree[1]
    if kind == "name":
        value = values.get(tree[1])
        return MISSING if value is None else value
    if kind == "not":
        inner = evaluate(tree[1], values)
        return MISSING if inner is MISSING else not bool(inner)
    if kind in ("and", "or"):
        left, right = evaluate(tree[1], values), evaluate(tree[2], values)
        decisive = kind == "or"  # the value that settles it: True for or, False for and
        sides = [x if x is MISSING else bool(x) for x in (left, right)]
        if decisive in sides:
            return decisive
        return MISSING if MISSING in sides else not decisive
    if kind == "arith":
        left, right = evaluate(tree[2], values), evaluate(tree[3], values)
        if left is MISSING or right is MISSING:
            return MISSING
        try:
            return _ARITH[tree[1]](left, right)
        except (TypeError, ZeroDivisionError):  # text, or a ratio to zero: unknown
            return MISSING
    left, right = evaluate(tree[2], values), evaluate(tree[3], values)
    if left is MISSING or right is MISSING:
        return MISSING
    try:
        return _COMPARE[tree[1]](left, right)
    except TypeError:  # a number against a string: not a match
        return False


def holds(text, values):
    """True only if the condition is known to hold."""
    return evaluate(parse(text), values) is True
