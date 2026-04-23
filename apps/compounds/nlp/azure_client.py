"""AzureOpenAI client factory (§11.1).

Isolated so it can be monkey-patched cleanly by tests and so the heavyweight
``openai`` / ``azure.identity`` imports stay lazy — importing the rest of
``compounds.nlp`` from a test shell shouldn't require the SDK to be installed.

Production deployments authenticate via the user-assigned managed identity
already bound to the container app (`containerAppsIdentity`); no API key in
env. The ``AZURE_OPENAI_API_KEY`` fallback is dev-only — handy when a
developer's AAD account hasn't been granted the Cognitive Services OpenAI
User role on the dev resource, but never set in production.
"""

from __future__ import annotations

import os
from typing import Any

_OPENAI_SCOPE = "https://cognitiveservices.azure.com/.default"

# Environment-variable names the factory reads.
ENV_ENDPOINT = "AZURE_OPENAI_ENDPOINT"
ENV_API_VERSION = "AZURE_OPENAI_API_VERSION"
ENV_MODEL = "AZURE_OPENAI_MODEL"
ENV_API_KEY = "AZURE_OPENAI_API_KEY"

DEFAULT_API_VERSION = "2024-10-21"
DEFAULT_MODEL = "gpt-4o"


def get_model_name() -> str:
    return os.environ.get(ENV_MODEL, DEFAULT_MODEL)


def get_azure_openai_client() -> Any:
    """Return a configured ``AzureOpenAI`` client.

    Raises ``RuntimeError`` if ``AZURE_OPENAI_ENDPOINT`` is unset — there is
    no sensible fallback for the endpoint.
    """
    endpoint = os.environ.get(ENV_ENDPOINT)
    if not endpoint:
        raise RuntimeError(
            f"{ENV_ENDPOINT} is not set; cannot construct AzureOpenAI client."
        )
    api_version = os.environ.get(ENV_API_VERSION, DEFAULT_API_VERSION)

    # Deferred imports — keeps module import cheap and lets tests mock this
    # factory without requiring the SDK to be installed locally.
    from openai import AzureOpenAI

    api_key = os.environ.get(ENV_API_KEY)
    if api_key:
        return AzureOpenAI(
            api_version=api_version,
            azure_endpoint=endpoint,
            api_key=api_key,
        )

    from azure.identity import DefaultAzureCredential, get_bearer_token_provider

    credential = DefaultAzureCredential()
    token_provider = get_bearer_token_provider(credential, _OPENAI_SCOPE)
    return AzureOpenAI(
        api_version=api_version,
        azure_endpoint=endpoint,
        azure_ad_token_provider=token_provider,
    )
