"""
The same invariants, over a container built the other way.

A CData tree comes into being by two independent routes: the def.xml parser,
and Phil2CData. Tiers 1 to 5 walk def.xml-built trees exclusively, so a change
to construction could hold for one and break the other with nothing to say so.

PHIL is the smaller route but the more searching test of the *contract*: it
builds containers without a def.xml, without qualifiers inherited from a task,
and with names mangled by scope (`refinement__resolution`). If the invariants
are properties of CData rather than of the parser, they hold here too.

Skips without libtbx, as the rest of the PHIL tier does.
"""
import pytest

parse = pytest.importorskip(
    "libtbx.phil", reason="needs libtbx (CCP4/cctbx)").parse

from ccp4i2.tests.unit.conformance import harness  # noqa: E402
from ccp4i2.utils.phil_to_cdata import Phil2CData  # noqa: E402


NESTED_PHIL = parse("""
    refinement {
      resolution = 2.0
        .type = float
      cycles = 10
        .type = int
      anisotropic = False
        .type = bool
      strategy = *individual_sites rigid_body
        .type = choice
      output {
        prefix = refined
          .type = str
        serial = 1
          .type = int
      }
    }
    scaling {
      resolution = 3.0
        .type = float
    }
""")


@pytest.fixture(name="container")
def container_fixture():
    return Phil2CData(NESTED_PHIL).convert()


def test_a_phil_container_is_walkable(container):
    walked = harness.Walk(task='phil')
    harness._visit(container, '', walked, 0, set())
    assert walked.paths, 'nothing was reached in a PHIL-built container'


def test_every_route_to_a_phil_parameter_agrees(container):
    """Tier 1's property, on the other construction route."""
    walked = harness.Walk(task='phil')
    harness._visit(container, '', walked, 0, set())
    assert not walked.divergences, (
        'access mechanisms disagree in a PHIL-built container:\n  '
        + '\n  '.join(str(d) for d in walked.divergences))


def test_phil_leaves_answer_the_questions_plugins_ask(container):
    """Tier 2's property: str(), isSet() and == must not raise."""
    walked = harness.Walk(task='phil')
    harness._visit(container, '', walked, 0, set())
    assert walked.leaves, 'a PHIL container with no leaves is not being measured'

    broken = []
    for path, obj in walked.leaves:
        for label, call in (('str', lambda o: str(o)),
                            ('isSet', lambda o: o.isSet()),
                            ('eq', lambda o: o == 'not the value')):
            try:
                call(obj)
            except AttributeError:
                if label == 'isSet':
                    continue          # not every leaf type carries it
                broken.append(f'{path}: {label}')
            except Exception as err:
                broken.append(f'{path}: {label}: {type(err).__name__}: {err}')
    assert not broken, f'{len(broken)} PHIL leaves misbehave:\n  ' + '\n  '.join(broken[:10])


def test_a_scope_becomes_a_container_and_a_definition_a_leaf(container):
    """The distinction tier 1 relies on, made by a different builder.

    A leaf that reported children is what broke phaser_MR.setKeywords, which
    treats a non-empty dataOrder() as "this is a sub-container".
    """
    walked = harness.Walk(task='phil')
    harness._visit(container, '', walked, 0, set())
    leaf_paths = {p for p, _o in walked.leaves}

    assert 'refinement' not in leaf_paths, 'a scope was taken for a leaf'
    assert any(p.startswith('refinement.') for p in leaf_paths), \
        'a scope yielded no leaves'


def test_two_scopes_may_hold_the_same_name(container):
    """`resolution` exists under both refinement and scaling.

    Name collisions across scopes are the PHIL-side equivalent of the
    inputData/outputData collisions that made i2run's tie-break matter.
    """
    walked = harness.Walk(task='phil')
    harness._visit(container, '', walked, 0, set())
    matching = [p for p, _o in walked.leaves if p.endswith('resolution')]
    assert len(matching) >= 2, f'expected the name in two scopes, got {matching}'
    assert len(set(matching)) == len(matching), 'the two collapsed into one'


def test_a_phil_container_survives_being_read_twice(container):
    """Reading must not perturb: __getattr__ can create."""
    first = harness.Walk(task='phil')
    harness._visit(container, '', first, 0, set())
    second = harness.Walk(task='phil')
    harness._visit(container, '', second, 0, set())
    assert first.paths == second.paths
