"""What an agent needs to drive CCP4i2 well: each task's judgement, read
from the task's ``<task>.agent.yaml``, and evaluated against its jobs.

See docs/agentic-knowledge.md.
"""


def mcp_main():
    """The ``i2-mcp`` command: the MCP facade, if its extra is installed."""
    try:
        from .mcp_server import main
    except ImportError as err:
        raise SystemExit(f"i2-mcp needs the MCP package ({err}): pip install 'ccp4i2[agent]'")
    main()
