"""Natural-language query subsystem for the compounds app.

See apps/compounds/docs/NLP_QUERY_PROPOSAL.md for the design.

Slice 1 (current): QuerySpec dataclass + target-field resolver against the
Gene-hydrated matching pool. No LLM, no executor, no UI, no Protocol/Metric
resolution. See the proposal's §5, §6.1, §7.
"""
