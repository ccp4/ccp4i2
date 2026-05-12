"""Shared API contract for CCP4i2 and consumers.

This package is the Python half of a bilingual API library; the TypeScript
half is published as ``@ccp4/ccp4i2-api``. Both halves implement the same
canonical bearer-token format, 401 response shape, and typed payloads
carried over the authenticated channel.

See the package README for current status and
``docs/CCP4I2_SERVICE_CONTRACT.md`` in the ccp4i2 monorepo for the
contract specification.
"""

__version__ = "0.3.0"
