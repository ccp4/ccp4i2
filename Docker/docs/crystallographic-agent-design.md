# Crystallographic AI Agent Design

## Overview

This document outlines the design for transforming CCP4i2 into a platform that supports an AI-powered crystallographic agent capable of autonomously working through structure solution challenges.

## Executive Summary

CCP4i2 already provides a robust foundation for AI-assisted crystallography:

- **176+ task plugins** covering the complete crystallographic pipeline
- **Structured performance metrics** that capture key quality indicators
- **Comprehensive REST API** for job management and execution
- **Real-time progress streaming** for monitoring running jobs
- **Workflow hints** via `WHATNEXT` attributes on each task

The primary work needed is building a **tool interface layer** that translates between the LLM's tool-calling format and the existing CCP4i2 API, plus **result interpretation** that converts raw XML outputs into agent-friendly summaries.

---

## Open Source Deployment Options

A key requirement for CCP4 is providing this functionality without requiring users to purchase additional licenses. The open-source AI ecosystem has matured significantly, making fully local, cost-free deployment viable.

### Recommended Stack

| Component | Role | License |
|-----------|------|---------|
| **Ollama** | Local LLM serving | MIT |
| **Llama 3.1/3.3** | Language model | Llama Community License (free for most uses) |
| **LangGraph** | Agent orchestration | MIT |
| **LangChain** | Tool integration | MIT |

### Local LLM Options

#### Ollama + Llama 3.1

[Ollama](https://ollama.com) provides the simplest path to running LLMs locally. It packages models with their dependencies and provides an OpenAI-compatible API.

```bash
# Install Ollama (macOS/Linux)
curl -fsSL https://ollama.com/install.sh | sh

# Pull Llama 3.1 (8B for most hardware, 70B for high-end)
ollama pull llama3.1:8b

# Or the smaller but capable 3.2
ollama pull llama3.2:3b
```

**Tool calling support**: [Ollama now natively supports tool/function calling](https://ollama.com/blog/tool-support) with Llama 3.1+, enabling the agent to invoke CCP4i2 tasks.

```python
import ollama

response = ollama.chat(
    model='llama3.1',
    messages=[{'role': 'user', 'content': 'Run refinement on job 123'}],
    tools=[{
        'type': 'function',
        'function': {
            'name': 'run_refinement',
            'description': 'Run REFMAC5 refinement',
            'parameters': {
                'type': 'object',
                'properties': {
                    'job_id': {'type': 'string'},
                    'num_cycles': {'type': 'integer'}
                }
            }
        }
    }]
)
```

#### Alternative: llama.cpp

For maximum control and performance, [llama.cpp](https://github.com/ggml-org/llama.cpp) provides direct model execution with [built-in function calling support](https://github.com/ggml-org/llama.cpp/blob/master/docs/function-calling.md) for Llama 3.1/3.2/3.3.

#### Alternative: vLLM

For production deployments with high throughput requirements, [vLLM](https://docs.vllm.ai/en/latest/features/tool_calling/) offers optimized inference with tool calling support.

### Model Selection Guide

| Model | Parameters | VRAM Required | Use Case |
|-------|------------|---------------|----------|
| Llama 3.2 3B | 3B | ~4GB | Laptops, quick responses |
| Llama 3.1 8B | 8B | ~8GB | Good balance of capability/speed |
| Llama 3.3 70B | 70B | ~48GB | Best open-source reasoning |
| Qwen 2.5 72B | 72B | ~48GB | Strong tool use |

For crystallographic reasoning, the **8B model is a reasonable minimum**; the 70B class models approach Claude-level capabilities for complex multi-step reasoning.

### Agent Framework: LangGraph

[LangGraph](https://www.langchain.com/langgraph) (MIT licensed) provides the orchestration layer for building stateful agents:

```python
from langgraph.graph import StateGraph
from langchain_community.llms import Ollama

# Use local Ollama model
llm = Ollama(model="llama3.1:8b")

# Define agent workflow
workflow = StateGraph(AgentState)
workflow.add_node("analyze", analyze_results)
workflow.add_node("decide", decide_next_step)
workflow.add_node("execute", execute_task)
workflow.add_edge("analyze", "decide")
workflow.add_conditional_edges("decide", route_to_task)

agent = workflow.compile()
```

Key LangGraph features for CCP4i2:
- **Stateful workflows**: Track structure solution progress across sessions
- **Human-in-the-loop**: Pause for user approval before expensive computations
- **Persistence**: Resume interrupted workflows
- **Tool integration**: 1000+ pre-built integrations, easy custom tools

### Architecture with Open Source Stack

```
┌─────────────────────────────────────────────────────────────────┐
│                     LangGraph Agent                              │
│               (Orchestration & State Management)                 │
└─────────────────────┬───────────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────────┐
│                   Ollama / llama.cpp                             │
│                (Local LLM - Llama 3.1/3.3)                       │
│                                                                  │
│  Capabilities:                                                   │
│  - Tool/function calling                                        │
│  - Multi-turn conversation                                       │
│  - Runs entirely offline                                        │
│  - No API costs                                                 │
└─────────────────────┬───────────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────────┐
│           CCP4i2 Tool Interface Layer                            │
│                                                                  │
│  - Tool schemas (OpenAI-compatible format)                      │
│  - Result interpreters                                          │
│  - Knowledge base                                                │
└─────────────────────┬───────────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────────┐
│           CCP4i2 Django REST API                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Integration Example

```python
# server/ccp4i2/agent/local_agent.py

from langchain_community.llms import Ollama
from langchain.agents import AgentExecutor, create_tool_calling_agent
from langchain_core.prompts import ChatPromptTemplate

from .tools import (
    run_refinement_tool,
    run_phaser_tool,
    get_job_metrics_tool,
    analyze_model_tool
)

def create_crystallographic_agent(model_name: str = "llama3.1:8b"):
    """
    Create a fully local crystallographic agent.

    No API keys or cloud services required.
    """
    # Local LLM via Ollama
    llm = Ollama(
        model=model_name,
        temperature=0.1,  # Lower temperature for more deterministic tool use
    )

    # Crystallographic tools
    tools = [
        run_refinement_tool,
        run_phaser_tool,
        get_job_metrics_tool,
        analyze_model_tool,
    ]

    # System prompt with crystallographic knowledge
    prompt = ChatPromptTemplate.from_messages([
        ("system", CRYSTALLOGRAPHIC_SYSTEM_PROMPT),
        ("human", "{input}"),
        ("placeholder", "{agent_scratchpad}")
    ])

    # Create agent
    agent = create_tool_calling_agent(llm, tools, prompt)

    return AgentExecutor(
        agent=agent,
        tools=tools,
        verbose=True,
        max_iterations=20,  # Allow multi-step workflows
        handle_parsing_errors=True
    )


# Usage
agent = create_crystallographic_agent()
result = agent.invoke({
    "input": "I have MR results in job 456. Please refine the model and "
             "tell me if the solution looks correct."
})
```

### Hardware Requirements

#### GPU Support

Ollama and llama.cpp support multiple GPU backends:

| Platform | Acceleration | Notes |
|----------|--------------|-------|
| **Apple Silicon (Metal)** | Native | Excellent support, uses unified memory |
| **NVIDIA (CUDA)** | Native | Best raw performance |
| **AMD (ROCm)** | Native | Linux primarily |
| **CPU only** | Yes | Slower but works everywhere |

#### Apple Silicon (Mac)

Apple Silicon Macs are particularly well-suited for local LLMs due to unified memory - the GPU can access all system RAM:

| Mac | Unified Memory | Recommended Model |
|-----|---------------|-------------------|
| M1/M2 (8GB) | 8GB | Llama 3.2 3B |
| M1/M2 Pro (16GB) | 16GB | Llama 3.1 8B |
| M1/M2/M3 Max (32-64GB) | 32-64GB | Llama 3.1 8B (fast), 70B (quantized) |
| M1/M2 Ultra (64-128GB) | 64-128GB | Llama 3.3 70B (full precision) |

#### General Requirements

| Configuration | Models Supported | Typical Use |
|---------------|-----------------|-------------|
| **Laptop** (16GB RAM, no GPU) | Llama 3.2 3B (CPU) | Basic assistance, simple queries |
| **Workstation** (32GB RAM, 8GB GPU) | Llama 3.1 8B | Full agent capability |
| **High-end** (64GB RAM, 24GB+ GPU) | Llama 3.3 70B (quantized) | Best reasoning |
| **Server** (128GB+ RAM, 48GB+ GPU) | Llama 3.3 70B (full) | Production deployment |

Most crystallographers with modern workstations (including Macs) can run the 8B model effectively.

### Comparison: Local vs Cloud

| Aspect | Local (Ollama + Llama) | Cloud (Claude API) |
|--------|------------------------|-------------------|
| **Cost** | Free (hardware only) | Per-token charges |
| **Privacy** | Data stays local | Data sent to cloud |
| **Latency** | Depends on hardware | Generally faster |
| **Reasoning** | Good (8B), Excellent (70B) | Excellent |
| **Setup** | Requires installation | API key only |
| **Offline** | Yes | No |

### Hybrid Approach

CCP4i2 could support both local and cloud models, letting users choose:

```python
# Configuration in settings
AGENT_CONFIG = {
    "backend": "ollama",  # or "anthropic", "openai"
    "model": "llama3.1:8b",  # or "claude-sonnet-4-20250514"
    "ollama_host": "http://localhost:11434",
}

def get_llm():
    """Get configured LLM backend."""
    if AGENT_CONFIG["backend"] == "ollama":
        return Ollama(
            model=AGENT_CONFIG["model"],
            base_url=AGENT_CONFIG["ollama_host"]
        )
    elif AGENT_CONFIG["backend"] == "anthropic":
        return ChatAnthropic(model=AGENT_CONFIG["model"])
```

### References

- [Ollama Tool Support](https://ollama.com/blog/tool-support)
- [LangGraph GitHub](https://github.com/langchain-ai/langgraph)
- [LangChain + LangGraph 1.0 Announcement](https://blog.langchain.com/langchain-langgraph-1dot0/)
- [Tool Calling in Llama 3 Guide](https://composio.dev/blog/tool-calling-in-llama-3-a-guide-to-build-agents)
- [Building Local AI Agents with LangGraph and Ollama](https://www.digitalocean.com/community/tutorials/local-ai-agents-with-langgraph-and-ollama)
- [llama.cpp Function Calling](https://github.com/ggml-org/llama.cpp/blob/master/docs/function-calling.md)
- [vLLM Tool Calling](https://docs.vllm.ai/en/latest/features/tool_calling/)
- [Open Source LLMs for 2025](https://blog.n8n.io/open-source-llm/)

---

## Global vs Local Optimization: The Refinement Duality

A fundamental insight for the agent architecture: crystallographic refinement requires both **global optimization** and **local problem-solving**, and these require fundamentally different tools.

### The Problem

Global refinement tools (Refmac, Servalcat) optimize all atoms simultaneously against diffraction data. They're excellent at:
- Minimizing R-factor and R-free
- Applying geometric restraints globally
- Optimizing B-factors and occupancies

But they **cannot fix local problems**:
- Wrong rotamers (the correct rotamer may be a local minimum away)
- Flipped peptides (180° rotation requires crossing an energy barrier)
- Misbuilt loops (topology is wrong, not just coordinates)
- Missing features (waters, ligands, alternate conformations)

Traditional refinement alternates between global and local:

```
┌─────────────────────────────────────────────────────────────────┐
│  Global Refinement (Refmac/Servalcat)                           │
│  - Optimizes all atoms simultaneously                            │
│  - Minimizes R-factor, R-free                                    │
│  - Applies geometric restraints globally                         │
│  - CANNOT fix local problems                                     │
└─────────────────────┬───────────────────────────────────────────┘
                      │ Produces maps + validation report
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│  Local Fixes (Interactive Coot session)                          │
│  - Identify outliers (Ramachandran, rotamers, density fit)       │
│  - Flip peptides, fix rotamers                                   │
│  - Rebuild loops, fit ligands                                    │
│  - Add/remove waters                                             │
└─────────────────────┬───────────────────────────────────────────┘
                      │ Improved model
                      ▼
                [Back to global refinement...]
```

An agent that only uses global tools will **plateau** - it cannot break through local minima.

### Coot Tools as CCP4i2 Jobs

The Crystallographic Object-Oriented Toolkit (Coot) provides the local problem-solving tools traditionally used in interactive graphics sessions. Critically, **Coot exposes all its automatic tools through an API**, which means they can be wrapped as standard CCP4i2 jobs. This dramatically simplifies the architecture - the agent uses the same job system for both global and local operations:

```
┌─────────────────────────────────────────────────────────────────┐
│                    Strategic Agent                               │
│                                                                  │
│  "R-free stalled at 0.28. Validation shows 5 Ramachandran       │
│  outliers. I should fix those with Coot before running          │
│  more Refmac cycles."                                            │
└─────────────────────────────────┬───────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────┐
│                    CCP4i2 Job System                             │
│                    (Unified Interface)                           │
├─────────────────────────────┬───────────────────────────────────┤
│  Global Optimization        │  Local Problem-Solving            │
│  (Batch jobs)               │  (Coot wrapper jobs)              │
│                             │                                   │
│  - refmac                   │  - coot_flip_peptide              │
│  - servalcat                │  - coot_auto_fit_rotamer          │
│  - modelcraft               │  - coot_real_space_refine         │
│  - phaser                   │  - coot_fit_loop                  │
│                             │  - coot_add_waters                │
│  Metrics:                   │  - coot_delete_sidechain          │
│  - R-factor, R-free         │                                   │
│  - Overall CC               │  Metrics:                         │
│  - Geometry RMS             │  - Per-residue RSCC               │
│                             │  - Ramachandran status            │
│                             │  - Rotamer status                 │
└─────────────────────────────┴───────────────────────────────────┘
                                  │
                                  ▼
                    Same API: POST /api/jobs/
                              POST /api/jobs/{id}/run/
                              GET /api/jobs/{id}/report_xml/
```

This unified approach means:
- **Same tool schema format** for global and local tools
- **Same execution pattern** (create job → set parameters → run → get results)
- **Same result interpretation** infrastructure
- **Job chaining** works naturally (Coot job takes output from Refmac job)
- **History tracking** through the existing job database

### Agent Workflow with Global/Local Duality

```python
def intelligent_refinement_loop(project_id, initial_model_job):
    """
    Refinement strategy combining global optimization and local fixes.

    This mimics expert crystallographer behavior:
    1. Run global refinement until it stalls
    2. Identify local problems via validation
    3. Fix local problems with Coot tools
    4. Return to global refinement
    """
    current_job = initial_model_job

    for iteration in range(max_iterations):
        # === GLOBAL PHASE ===
        refine_result = run_refinement(current_job)

        # Check if global refinement is still making progress
        if not significant_improvement(refine_result):

            # === DIAGNOSTIC PHASE ===
            validation = analyze_model(refine_result)
            local_problems = identify_local_problems(validation)

            if not local_problems:
                # Model is as good as we can make it
                return refine_result

            # === LOCAL FIX PHASE ===
            for problem in prioritize_problems(local_problems):

                if problem.type == "ramachandran_outlier":
                    # Peptide flip often fixes Ramachandran outliers
                    coot_flip_peptide(problem.chain, problem.residue)

                elif problem.type == "rotamer_outlier":
                    # Search for better rotamer fit to density
                    coot_auto_fit_rotamer(problem.chain, problem.residue)

                elif problem.type == "poor_density_fit":
                    # Local real-space refinement
                    coot_real_space_refine_zone(
                        problem.chain,
                        problem.residue - 2,
                        problem.residue + 2
                    )

                elif problem.type == "missing_sidechain_density":
                    # Truncate to Ala - don't model what isn't there
                    coot_delete_sidechain(problem.chain, problem.residue)

                elif problem.type == "unmodeled_water":
                    # Add water at difference density peak
                    coot_add_water(problem.position)

            # Model changed - return to global refinement
            current_job = save_model_and_create_refinement_job()

        else:
            # Still improving globally - continue
            current_job = refine_result

    return current_job
```

### Coot Tool Definitions

```python
COOT_LOCAL_TOOLS = [
    {
        "name": "coot_flip_peptide",
        "description": """Flip a peptide bond by 180 degrees.

        USE WHEN:
        - Ramachandran outlier that could be fixed by flipping
        - After MR when chain direction may be ambiguous
        - Residue is in a disallowed region of Ramachandran plot

        EFFECT: Dramatically changes phi/psi of this and neighboring
        residues. Often converts outlier to favored region.

        CAUTION: Will affect density fit - follow with real-space refine.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "model_job_id": {"type": "string", "description": "Job with model"},
                "chain": {"type": "string", "description": "Chain ID (e.g., 'A')"},
                "residue_number": {"type": "integer"}
            },
            "required": ["model_job_id", "chain", "residue_number"]
        }
    },

    {
        "name": "coot_auto_fit_rotamer",
        "description": """Automatically fit the best rotamer for a residue.

        USE WHEN:
        - Rotamer outlier flagged by validation
        - Sidechain doesn't fit density well
        - After backbone changes that may invalidate sidechain position

        EFFECT: Searches rotamer library, picks best fit to electron
        density. May significantly improve local geometry.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "model_job_id": {"type": "string"},
                "chain": {"type": "string"},
                "residue_number": {"type": "integer"}
            },
            "required": ["model_job_id", "chain", "residue_number"]
        }
    },

    {
        "name": "coot_real_space_refine_zone",
        "description": """Real-space refinement of a local zone.

        USE WHEN:
        - Region has correct topology but poor density fit
        - After rotamer/peptide changes to optimize local geometry
        - Backbone needs adjustment to fit density

        EFFECT: Optimizes atom positions against electron density
        while maintaining reasonable geometry. Affects all atoms
        in the zone.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "model_job_id": {"type": "string"},
                "chain": {"type": "string"},
                "start_residue": {"type": "integer"},
                "end_residue": {"type": "integer"}
            },
            "required": ["model_job_id", "chain", "start_residue", "end_residue"]
        }
    },

    {
        "name": "coot_fit_loop",
        "description": """Build or rebuild a loop region.

        USE WHEN:
        - Gap in the model where loop should be
        - Existing loop doesn't fit density
        - Loop traced with wrong connectivity

        EFFECT: Uses loop database or ab initio building to
        generate conformations, selects best fit to density.

        NOTE: May fail if density is too weak - loops are often
        flexible and poorly ordered.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "model_job_id": {"type": "string"},
                "chain": {"type": "string"},
                "start_residue": {"type": "integer", "description": "Anchor before loop"},
                "end_residue": {"type": "integer", "description": "Anchor after loop"},
                "sequence": {"type": "string", "description": "Loop sequence if known"}
            },
            "required": ["model_job_id", "chain", "start_residue", "end_residue"]
        }
    },

    {
        "name": "coot_add_water",
        "description": """Add water molecule(s) at density peaks.

        USE WHEN:
        - Positive Fo-Fc density peaks > 3σ
        - Peak is roughly spherical
        - Peak is at hydrogen-bonding distance from protein/other waters

        EFFECT: Adds HOH residues. Each water typically improves
        R-factors by 0.001-0.01. Ordered waters are real features.

        MODE: 'auto' finds all suitable peaks; or specify position.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "model_job_id": {"type": "string"},
                "mode": {
                    "type": "string",
                    "enum": ["auto", "position"],
                    "default": "auto"
                },
                "position": {
                    "type": "array",
                    "items": {"type": "number"},
                    "description": "[x, y, z] if mode is 'position'"
                },
                "sigma_cutoff": {"type": "number", "default": 3.0}
            },
            "required": ["model_job_id"]
        }
    },

    {
        "name": "coot_delete_sidechain",
        "description": """Truncate residue to alanine (remove sidechain).

        USE WHEN:
        - No density for sidechain atoms
        - Sidechain is clearly disordered
        - High B-factors on sidechain with no density support

        EFFECT: Removes atoms beyond CB. Better to have no atoms
        than wrong atoms - reduces overfitting and may improve R-free.

        NOTE: This is good practice for exposed, flexible residues
        where the sidechain genuinely is disordered.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "model_job_id": {"type": "string"},
                "chain": {"type": "string"},
                "residue_number": {"type": "integer"}
            },
            "required": ["model_job_id", "chain", "residue_number"]
        }
    },

    {
        "name": "coot_regularize_zone",
        "description": """Idealize geometry in a zone without density fitting.

        USE WHEN:
        - After manual edits that distorted geometry
        - To clean up geometry before refinement
        - Zone has good density fit but poor geometry

        EFFECT: Adjusts atom positions to ideal bond lengths and
        angles. Does NOT consider electron density.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "model_job_id": {"type": "string"},
                "chain": {"type": "string"},
                "start_residue": {"type": "integer"},
                "end_residue": {"type": "integer"}
            },
            "required": ["model_job_id", "chain", "start_residue", "end_residue"]
        }
    }
]
```

### Local Problem Prioritization

The agent needs heuristics for which problems to fix first:

```python
LOCAL_PROBLEM_PRIORITY = {
    "peptide_flip_candidate": {
        "priority": 1,  # Highest - often big impact
        "description": "Ramachandran outlier likely fixable by 180° flip",
        "detection": "Outlier where flipped phi/psi would be in favored region",
        "tool": "coot_flip_peptide",
        "expected_impact": "Removes outlier, may improve R by 0.002-0.01",
        "risk": "Low - easy to undo if wrong"
    },

    "rotamer_outlier": {
        "priority": 2,
        "description": "Sidechain in non-preferred rotamer with density elsewhere",
        "detection": "Rotamer outlier + density blob in alternative position",
        "tool": "coot_auto_fit_rotamer",
        "expected_impact": "Fixes geometry, improves local density fit",
        "risk": "Low"
    },

    "poor_backbone_density_fit": {
        "priority": 3,
        "description": "Backbone atoms don't match density",
        "detection": "Main chain RSR > 0.4 or RSCC < 0.7",
        "tool": "coot_real_space_refine_zone",
        "expected_impact": "Better density fit, may help R-factors",
        "risk": "Medium - can distort geometry if density is weak"
    },

    "unmodeled_water_peaks": {
        "priority": 4,
        "description": "Clear positive difference density for waters",
        "detection": "Spherical Fo-Fc peaks > 3σ at H-bond distances",
        "tool": "coot_add_water",
        "expected_impact": "~0.001-0.01 R improvement per water",
        "risk": "Low - easy to remove if wrong"
    },

    "disordered_sidechain": {
        "priority": 5,
        "description": "Sidechain modeled but no density",
        "detection": "Sidechain RSCC < 0.5 and B-factor > 60",
        "tool": "coot_delete_sidechain",
        "expected_impact": "May improve R-free by reducing overfitting",
        "risk": "Low - can rebuild later if density appears"
    },

    "broken_loop": {
        "priority": 6,  # Lower - often difficult
        "description": "Loop region with poor or no density",
        "detection": "Continuous region with RSR > 0.5 or gaps",
        "tool": "coot_fit_loop",
        "expected_impact": "Highly variable - depends on density quality",
        "risk": "High - may need human judgment"
    }
}
```

### Integration Architecture

Since Coot tools are wrapped as CCP4i2 jobs, the architecture is beautifully unified:

```
┌─────────────────────────────────────────────────────────────────┐
│                    LangGraph Orchestrator                        │
│                                                                  │
│  States: ANALYZE → DECIDE → EXECUTE → INTERPRET                 │
└─────────────────────────────────┬───────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────┐
│                 CCP4i2 Tool Interface Layer                      │
│                                                                  │
│  - Auto-generated tool schemas from plugin metadata              │
│  - Result interpreters for all task types                        │
│  - Crystallographic knowledge base                               │
└─────────────────────────────────┬───────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────┐
│                    CCP4i2 REST API                               │
│                                                                  │
│  POST /api/jobs/                    Create any job type          │
│  POST /api/jobs/{id}/set_context_job/   Chain from previous job  │
│  POST /api/jobs/{id}/run/           Execute (global or local)    │
│  GET /api/jobs/{id}/report_xml/     Get results                  │
│  GET /api/jobs/{id}/what_next/      Get suggestions              │
└─────────────────────────────────┬───────────────────────────────┘
                                  │
                    ┌─────────────┴─────────────┐
                    │                           │
                    ▼                           ▼
        ┌───────────────────────┐   ┌───────────────────────┐
        │  Global Tasks         │   │  Local Tasks          │
        │  (existing wrappers)  │   │  (Coot wrappers)      │
        │                       │   │                       │
        │  refmac               │   │  coot_flip_peptide    │
        │  servalcat            │   │  coot_auto_fit_rotamer│
        │  modelcraft           │   │  coot_real_space_refine│
        │  phaser               │   │  coot_fit_loop        │
        │  aimless              │   │  coot_add_waters      │
        └───────────────────────┘   └───────────────────────┘
```

**Key benefit**: The agent doesn't need to know the difference between global and local tools - it just creates jobs, runs them, and interprets results. The `WHATNEXT` attribute on each task can suggest appropriate follow-ups (e.g., `coot_flip_peptide.WHATNEXT = ['refmac']`).

### Coot Session Architecture

A critical insight: Coot maintains **in-memory state** (model coordinates, maps, etc.) and your colleague's implementation uses **RPC/JSON** to communicate with it. This is far more efficient than file-based job wrapping because:

1. **No file I/O per operation** - Model stays in memory
2. **No map recalculation** - Maps persist between operations
3. **Fast iteration** - Multiple fixes can be applied in milliseconds
4. **Undo support** - Coot's internal undo stack is available

The agent should interact with a **persistent Coot session** for local fixes, only creating CCP4i2 jobs at meaningful checkpoints (e.g., "after fixing all Ramachandran outliers, save and refine").

### Coot-Centric Architecture (Recommended)

A radical but elegant simplification: **Coot is the agent's primary interface**, and Coot calls CCP4i2's REST API when global operations are needed.

This matches the crystallographer's mental model: "I'm working in Coot, and I send the model off for refinement when ready."

```
┌─────────────────────────────────────────────────────────────────┐
│                    LLM Agent (Ollama/Claude)                     │
│                                                                  │
│  Strategic reasoning:                                            │
│  - "R-free stalled, fix outliers then refine"                   │
│  - "Model incomplete, run ModelCraft"                           │
│  - "Structure solved, prepare for deposition"                   │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      │ JSON-RPC
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Coot RPC Server                               │
│                    (Primary Agent Interface)                     │
│                                                                  │
│  In-memory state:                                                │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐              │
│  │   Model     │  │    Maps     │  │ Validation  │              │
│  │ coordinates │  │  2Fo-Fc     │  │  outliers   │              │
│  │             │  │  Fo-Fc      │  │  geometry   │              │
│  └─────────────┘  └─────────────┘  └─────────────┘              │
│                                                                  │
│  Local operations (fast, in-memory):                            │
│  - flip_peptide()         - auto_fit_rotamer()                  │
│  - real_space_refine()    - add_waters()                        │
│  - delete_sidechain()     - fit_loop()                          │
│  - undo() / redo()        - get_validation()                    │
│                                                                  │
│  Global operations (delegates to CCP4i2):                       │
│  - refine()  ──────────────────┐                                │
│  - build()   ──────────────────┼──► CCP4i2 REST API             │
│  - phase()   ──────────────────┘                                │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      │ REST API calls
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                    CCP4i2 Django Server                          │
│                                                                  │
│  POST /api/jobs/                    Create job                   │
│  POST /api/jobs/{id}/run/           Execute                      │
│  GET  /api/jobs/{id}/status/        Poll completion              │
│  GET  /api/jobs/{id}/report_xml/    Get results                  │
│                                                                  │
│  Heavy batch operations:                                         │
│  - refmac / servalcat (global refinement)                       │
│  - modelcraft / buccaneer (model building)                      │
│  - phaser (molecular replacement)                               │
│  - aimless (data processing)                                    │
│                                                                  │
│  Provides:                                                       │
│  - Job history & audit trail                                    │
│  - Output file management                                        │
│  - Performance metrics                                           │
└─────────────────────────────────────────────────────────────────┘
```

### Why This Architecture Works

1. **Single source of truth**: The model lives in Coot's memory. No synchronization issues.

2. **Natural fit**: Crystallographers think "I'm in Coot, I'll refine this". The agent works the same way.

3. **Efficient**: Local operations are instant (RPC). Global operations are naturally batch (REST API).

4. **Clean separation**:
   - Coot: Interactive, in-memory, fast feedback
   - CCP4i2: Batch processing, file management, audit trail

5. **Existing infrastructure**: Your colleague's Coot RPC work + CCP4i2's REST API = complete system.

### Agent Tool Interface (Coot-Centric)

The agent talks only to Coot. Coot handles everything:

```python
COOT_AGENT_TOOLS = [
    # === Session Management ===
    {
        "name": "load_structure",
        "description": "Load a model and maps into the Coot session",
        "params": ["pdb_path", "mtz_path"]
    },

    # === Queries (instant) ===
    {
        "name": "get_validation",
        "description": "Get current validation: Ramachandran, rotamers, density fit",
        "params": []
    },
    {
        "name": "get_metrics",
        "description": "Get current R-factors and geometry statistics",
        "params": []
    },

    # === Local Fixes (fast, in-memory) ===
    {
        "name": "flip_peptide",
        "description": "Flip peptide bond 180°",
        "params": ["chain", "residue"]
    },
    {
        "name": "auto_fit_rotamer",
        "description": "Fit best rotamer to density",
        "params": ["chain", "residue"]
    },
    {
        "name": "real_space_refine_zone",
        "description": "Refine a zone against density",
        "params": ["chain", "start", "end"]
    },
    {
        "name": "add_waters",
        "description": "Add waters at difference peaks",
        "params": ["sigma_cutoff"]
    },
    {
        "name": "undo",
        "description": "Undo last operation",
        "params": []
    },

    # === Global Operations (delegates to CCP4i2) ===
    {
        "name": "refine",
        "description": """Run global refinement (Refmac/Servalcat).

        Coot saves current model, calls CCP4i2 API, waits for completion,
        loads refined model and new maps back into session.

        This is a blocking call - may take 30-120 seconds.""",
        "params": ["program", "cycles"]
    },
    {
        "name": "auto_build",
        "description": """Run automated model building (ModelCraft).

        Coot saves current state, calls CCP4i2 API, loads result.
        May take several minutes.""",
        "params": ["sequence_file"]
    },

    # === Checkpoint ===
    {
        "name": "save_checkpoint",
        "description": "Save current model to a named checkpoint",
        "params": ["name"]
    }
]
```

### Coot's CCP4i2 Client

Coot needs a thin client to call CCP4i2's API:

```python
# Inside Coot's Python environment

class CCP4i2Client:
    """Client for Coot to call CCP4i2 API for global operations."""

    def __init__(self, base_url: str = "http://localhost:8000/api"):
        self.base_url = base_url
        self.project_id = None

    def refine(self, model_path: str, mtz_path: str,
               program: str = "refmac", cycles: int = 10) -> dict:
        """
        Run global refinement via CCP4i2.

        1. Create job
        2. Upload/reference input files
        3. Run and wait for completion
        4. Return output file paths
        """
        # Create refinement job
        job = requests.post(f"{self.base_url}/jobs/", json={
            "project": self.project_id,
            "task_name": program,
        }).json()

        job_id = job["uuid"]

        # Set input files
        requests.post(f"{self.base_url}/jobs/{job_id}/set_parameter/", json={
            "path": "inputData.XYZIN",
            "value": model_path
        })
        requests.post(f"{self.base_url}/jobs/{job_id}/set_parameter/", json={
            "path": "inputData.F_SIGF",
            "value": mtz_path
        })
        requests.post(f"{self.base_url}/jobs/{job_id}/set_parameter/", json={
            "path": "controlParameters.NCYCLES",
            "value": cycles
        })

        # Run job
        requests.post(f"{self.base_url}/jobs/{job_id}/run/")

        # Poll for completion
        while True:
            status = requests.get(f"{self.base_url}/jobs/{job_id}/").json()
            if status["status"] in ["Finished", "Failed"]:
                break
            time.sleep(2)

        if status["status"] == "Failed":
            raise RuntimeError(f"Refinement failed: {status.get('error')}")

        # Get output paths
        container = requests.get(
            f"{self.base_url}/jobs/{job_id}/container/"
        ).json()

        return {
            "model": container["outputData"]["XYZOUT"]["fullPath"],
            "map_coeffs": container["outputData"]["FPHIOUT"]["fullPath"],
            "diff_map": container["outputData"]["DIFFPHIOUT"]["fullPath"],
            "r_factor": container["performance"]["RFactor"],
            "r_free": container["performance"]["RFree"],
        }
```

### Coot RPC Methods for Global Operations

The Coot RPC server exposes methods that orchestrate the full cycle:

```python
# In Coot's RPC server

class CootAgentRPC:
    """RPC methods for agent interaction."""

    def __init__(self):
        self.ccp4i2 = CCP4i2Client()
        self.model_mol = None
        self.map_mol = None

    def refine(self, program: str = "refmac", cycles: int = 10) -> dict:
        """
        Run global refinement and reload results.

        Agent calls this single method - Coot handles all the orchestration.
        """
        # Save current model to temp file
        temp_pdb = self._save_temp_model()

        # Get current map coefficients path (or extract from loaded map)
        temp_mtz = self._get_map_coefficients_path()

        # Call CCP4i2 for refinement
        result = self.ccp4i2.refine(
            model_path=temp_pdb,
            mtz_path=temp_mtz,
            program=program,
            cycles=cycles
        )

        # Load refined model (replaces current)
        coot.close_molecule(self.model_mol)
        self.model_mol = coot.read_pdb(result["model"])

        # Load new maps
        coot.close_molecule(self.map_mol)
        self.map_mol = coot.auto_read_mtz(result["map_coeffs"])

        # Return summary for agent
        return {
            "success": True,
            "r_factor": result["r_factor"],
            "r_free": result["r_free"],
            "message": f"Refined with {program}: R={result['r_factor']:.3f}, R-free={result['r_free']:.3f}"
        }

    def get_validation(self) -> dict:
        """Get current model validation status."""
        return {
            "ramachandran_outliers": self._get_rama_outliers(),
            "rotamer_outliers": self._get_rotamer_outliers(),
            "density_fit_problems": self._get_density_problems(),
            "clashes": self._get_clashes(),
        }

    def flip_peptide(self, chain: str, residue: int) -> dict:
        """Flip peptide and return before/after status."""
        before = self._get_residue_rama(chain, residue)
        coot.pepflip(self.model_mol, chain, residue, "")
        after = self._get_residue_rama(chain, residue)

        return {
            "operation": "flip_peptide",
            "residue": f"{chain}/{residue}",
            "before": before,
            "after": after,
            "improved": after != "OUTLIER"
        }

    # ... other local operations ...
```

### Example Agent Session

```
Agent: I'll analyze the current model and improve it.

[Calls: get_validation()]

Current validation:
- 3 Ramachandran outliers: A/45, A/123, B/89
- 5 rotamer outliers
- 12 residues with poor density fit
- R-factor: 0.28, R-free: 0.32

I'll fix the Ramachandran outliers first, then refine.

[Calls: flip_peptide(chain="A", residue=45)]
→ A/45: OUTLIER → FAVORED ✓

[Calls: flip_peptide(chain="A", residue=123)]
→ A/123: OUTLIER → ALLOWED ✓

[Calls: flip_peptide(chain="B", residue=89)]
→ B/89: OUTLIER → OUTLIER (flip didn't help)

[Calls: undo()]
[Calls: real_space_refine_zone(chain="B", start=87, end=91)]
→ Zone B/87-91: RSCC 0.65 → 0.78, Ramachandran now ALLOWED ✓

Now I'll fix rotamer outliers...
[Calls: auto_fit_rotamer for 5 residues...]

All local fixes applied. Running global refinement:

[Calls: refine(program="refmac", cycles=10)]

Refinement complete:
- R-factor: 0.28 → 0.24
- R-free: 0.32 → 0.27

Excellent improvement! The local fixes allowed global refinement to
find a better minimum. Validation now shows 0 Ramachandran outliers.
```

### Architecture Benefits

| Aspect | Coot-Centric | CCP4i2-Centric |
|--------|--------------|----------------|
| **Primary interface** | Coot RPC | CCP4i2 REST API |
| **Model state** | In Coot memory | In files |
| **Local operations** | Native, instant | Wrapped, slow |
| **Global operations** | Via CCP4i2 API | Native |
| **Undo** | Native | Complex |
| **Mental model** | "Working in Coot" | "Running jobs" |
| **Implementation** | Uses existing RPC work | New wrappers needed |

### Resilience: Handling Coot Crashes

Stateful services are fragile. Coot can crash from memory leaks, corrupted models, or edge cases in local operations. The architecture must handle this gracefully.

#### The Risk

If Coot "falls over", you lose:
- Current model state (if not checkpointed recently)
- Undo history
- Calculated maps
- In-flight operation results

#### Mitigation Strategy: Frequent Checkpoints to CCP4i2

The key insight: **CCP4i2 is already a persistent, stateless job store**. Use it as the recovery mechanism.

```python
class ResilientCootSession:
    """
    Coot session wrapper with automatic checkpointing and crash recovery.
    """

    def __init__(self, coot_rpc_url: str, ccp4i2_client: CCP4i2Client):
        self.coot_url = coot_rpc_url
        self.ccp4i2 = ccp4i2_client

        # Checkpoint tracking
        self.last_checkpoint_job_id: str = None
        self.operations_since_checkpoint: int = 0
        self.checkpoint_threshold: int = 5  # Checkpoint every N operations

    def _ensure_coot_alive(self) -> bool:
        """Check if Coot RPC is responding."""
        try:
            response = requests.get(f"{self.coot_url}/health", timeout=2)
            return response.ok
        except:
            return False

    def _restart_coot_from_checkpoint(self):
        """Restart Coot and reload from last checkpoint."""
        # Start new Coot instance (implementation depends on deployment)
        self._start_coot_process()

        if self.last_checkpoint_job_id:
            # Get checkpoint files from CCP4i2
            job = self.ccp4i2.get_job(self.last_checkpoint_job_id)
            model_path = job["outputData"]["XYZOUT"]["fullPath"]
            mtz_path = job["outputData"]["FPHIOUT"]["fullPath"]

            # Reload into fresh Coot
            self._rpc_call("read_pdb", {"filename": model_path})
            self._rpc_call("auto_read_mtz", {"filename": mtz_path})

            return {
                "recovered": True,
                "from_checkpoint": self.last_checkpoint_job_id,
                "operations_lost": self.operations_since_checkpoint
            }

        return {"recovered": False, "error": "No checkpoint available"}

    def _auto_checkpoint_if_needed(self):
        """Checkpoint to CCP4i2 if threshold reached."""
        if self.operations_since_checkpoint >= self.checkpoint_threshold:
            self._create_checkpoint("auto")

    def _create_checkpoint(self, name: str) -> str:
        """Save current state to CCP4i2 job."""
        # Save model to temp file
        temp_pdb = f"/tmp/coot_checkpoint_{name}.pdb"
        self._rpc_call("write_pdb_file", {"imol": 0, "filename": temp_pdb})

        # Create CCP4i2 job to store checkpoint
        job = self.ccp4i2.create_job(
            task_name="coot_checkpoint",
            metadata={
                "checkpoint_name": name,
                "operations_count": self.operations_since_checkpoint,
                "timestamp": datetime.now().isoformat()
            }
        )

        # Copy model to job directory
        job_model_path = self.ccp4i2.store_file(job["uuid"], temp_pdb, "XYZOUT")

        self.last_checkpoint_job_id = job["uuid"]
        self.operations_since_checkpoint = 0

        return job["uuid"]

    def execute_operation(self, method: str, params: dict) -> dict:
        """
        Execute a Coot operation with crash recovery.
        """
        # Check Coot is alive
        if not self._ensure_coot_alive():
            recovery = self._restart_coot_from_checkpoint()
            if not recovery["recovered"]:
                raise RuntimeError("Coot crashed and no checkpoint available")

            # Inform caller about recovery
            return {
                "recovered_from_crash": True,
                "operations_lost": recovery["operations_lost"],
                "message": "Coot restarted from checkpoint. Please retry operation."
            }

        # Execute the operation
        try:
            result = self._rpc_call(method, params)
            self.operations_since_checkpoint += 1
            self._auto_checkpoint_if_needed()
            return result

        except requests.exceptions.ConnectionError:
            # Coot died during operation
            recovery = self._restart_coot_from_checkpoint()
            return {
                "recovered_from_crash": True,
                "operation_failed": method,
                "message": f"Coot crashed during {method}. Recovered from checkpoint."
            }
```

#### Checkpoint Strategy

| Event | Action |
|-------|--------|
| **After N local operations** | Auto-checkpoint (configurable, default 5) |
| **Before global refinement** | Mandatory checkpoint |
| **After global refinement** | Checkpoint (new model from CCP4i2) |
| **User requests checkpoint** | Explicit checkpoint |
| **Coot crash detected** | Restart + reload from last checkpoint |

#### Agent Awareness

The agent should be aware of checkpointing and handle recovery gracefully:

```python
{
    "name": "execute_coot_operation",
    "description": """Execute a Coot operation with automatic crash recovery.

    If Coot crashes, it will be restarted from the last checkpoint.
    Operations since the checkpoint will be lost and need to be retried.

    The response may include:
    - recovered_from_crash: true if Coot was restarted
    - operations_lost: number of operations lost since checkpoint

    IMPORTANT: If recovery occurs, re-query validation state before
    continuing, as some fixes may have been lost."""
}
```

Example agent handling:

```
[Calls: flip_peptide(chain="A", residue=45)]

Response: {
    "recovered_from_crash": true,
    "operations_lost": 3,
    "message": "Coot crashed during flip_peptide. Recovered from checkpoint."
}

Agent: Coot crashed and was restarted. I lost 3 operations since the last
checkpoint. Let me re-check the current validation state and redo the fixes.

[Calls: get_validation()]

Current validation shows 3 Ramachandran outliers (same as before my fixes).
I'll redo the peptide flips that were lost...
```

#### Stateless Alternative: Per-Operation File I/O

For maximum resilience (at the cost of performance), operations can be made stateless:

```python
class StatelessCootOperations:
    """
    Each operation reads model, modifies, writes result.
    No persistent Coot session - maximum reliability.
    """

    def flip_peptide(self, model_path: str, mtz_path: str,
                     chain: str, residue: int) -> dict:
        """
        Stateless peptide flip.

        1. Start fresh Coot process
        2. Load model and maps
        3. Perform flip
        4. Save result
        5. Exit Coot

        Slow but crash-proof.
        """
        output_path = f"{model_path}.flipped.pdb"

        script = f'''
import coot
mol = coot.read_pdb("{model_path}")
coot.auto_read_mtz("{mtz_path}")

# Get before state
before = get_ramachandran_status(mol, "{chain}", {residue})

# Flip
coot.pepflip(mol, "{chain}", {residue}, "")

# Get after state
after = get_ramachandran_status(mol, "{chain}", {residue})

# Save
coot.write_pdb_file(mol, "{output_path}")

# Write result JSON
import json
with open("{output_path}.result.json", "w") as f:
    json.dump({{"before": before, "after": after}}, f)

coot.coot_no_state_real_exit(0)
'''

        # Run Coot with script
        subprocess.run(["coot", "--no-graphics", "-s", "-"],
                       input=script, capture_output=True)

        # Read result
        with open(f"{output_path}.result.json") as f:
            return json.load(f)
```

#### Hybrid: Session with Aggressive Checkpointing

The recommended approach balances performance with resilience:

| Aspect | Session Mode | Stateless Mode |
|--------|--------------|----------------|
| **Performance** | Fast (in-memory) | Slow (file I/O each op) |
| **Crash impact** | Lose ops since checkpoint | None (each op independent) |
| **Complexity** | Checkpoint management | Simple |
| **Recommended for** | Interactive refinement | Critical/unattended ops |

```python
COOT_MODE = {
    "session": {
        "use_when": "Interactive work, user present",
        "checkpoint_every": 5,
        "pros": "Fast, undo support",
        "cons": "Can lose recent work on crash"
    },
    "stateless": {
        "use_when": "Unattended batch processing, overnight runs",
        "checkpoint_every": 1,  # Each operation is its own checkpoint
        "pros": "Crash-proof",
        "cons": "Slow, no undo"
    }
}
```

### Fallback: File-Based Wrapper

For environments where Coot RPC isn't available, a simplified file-based wrapper can be used as a fallback. See appendix for the `coot_autofix` CCP4i2 wrapper implementation.

### Local Fix Performance Metrics

A critical requirement: each Coot wrapper must report **before/after quality indicators** so the agent can assess whether the fix helped. This requires new performance classes:

```python
# core/cdata_stubs/CCP4PerformanceData.py (additions)

@cdata_class(
    attributes={
        "residue_id": attribute(AttributeType.STRING),  # "A/45"
        "fix_type": attribute(AttributeType.STRING),    # "peptide_flip", "rotamer_fit", etc.

        # Before metrics
        "ramachandran_before": attribute(AttributeType.STRING),  # "OUTLIER", "ALLOWED", "FAVORED"
        "rotamer_before": attribute(AttributeType.STRING),       # "OUTLIER", "ALLOWED", "FAVORED"
        "rscc_before": attribute(AttributeType.FLOAT),           # Real-space CC (0-1)
        "rsr_before": attribute(AttributeType.FLOAT),            # Real-space R (lower=better)

        # After metrics
        "ramachandran_after": attribute(AttributeType.STRING),
        "rotamer_after": attribute(AttributeType.STRING),
        "rscc_after": attribute(AttributeType.FLOAT),
        "rsr_after": attribute(AttributeType.FLOAT),

        # Summary
        "success": attribute(AttributeType.BOOL),
        "improvement_summary": attribute(AttributeType.STRING),
    }
)
class CLocalFixPerformance(CPerformanceIndicatorStub):
    """Performance metrics for local model fixes (Coot operations)."""
    pass


@cdata_class(
    attributes={
        "waters_added": attribute(AttributeType.INT),
        "waters_rejected": attribute(AttributeType.INT),
        "average_peak_sigma": attribute(AttributeType.FLOAT),
        "average_distance_to_protein": attribute(AttributeType.FLOAT),
    }
)
class CWaterAdditionPerformance(CPerformanceIndicatorStub):
    """Performance metrics for automated water addition."""
    pass


@cdata_class(
    attributes={
        "zone_start": attribute(AttributeType.STRING),  # "A/45"
        "zone_end": attribute(AttributeType.STRING),    # "A/52"
        "residues_affected": attribute(AttributeType.INT),
        "rscc_before": attribute(AttributeType.FLOAT),
        "rscc_after": attribute(AttributeType.FLOAT),
        "geometry_score_before": attribute(AttributeType.FLOAT),
        "geometry_score_after": attribute(AttributeType.FLOAT),
    }
)
class CZoneRefinePerformance(CPerformanceIndicatorStub):
    """Performance metrics for real-space zone refinement."""
    pass
```

Each Coot wrapper then populates these in `processOutputFiles()`:

```python
# wrappers/coot_flip_peptide/script/coot_flip_peptide.py

class coot_flip_peptide(CPluginScript):
    PERFORMANCECLASS = 'CLocalFixPerformance'

    def processOutputFiles(self):
        # Calculate before/after metrics
        chain = self.container.controlParameters.CHAIN
        resno = self.container.controlParameters.RESIDUE

        # These would use gemmi or Coot's own validation
        before_rama = self.stored_before_ramachandran  # Saved in processInputFiles
        after_rama = calculate_ramachandran(
            self.container.outputData.XYZOUT.fullPath,
            chain, resno
        )

        before_rscc = self.stored_before_rscc
        after_rscc = calculate_rscc(
            self.container.outputData.XYZOUT.fullPath,
            self.container.inputData.FPHIOUT.fullPath,  # Map coefficients
            chain, resno
        )

        # Populate performance indicators
        perf = self.container.outputData.performance
        perf.residue_id = f"{chain}/{resno}"
        perf.fix_type = "peptide_flip"
        perf.ramachandran_before = before_rama
        perf.ramachandran_after = after_rama
        perf.rscc_before = before_rscc
        perf.rscc_after = after_rscc

        # Determine success
        improved_rama = (before_rama == "OUTLIER" and after_rama != "OUTLIER")
        improved_density = (after_rscc > before_rscc)

        perf.success = improved_rama or improved_density

        if perf.success:
            perf.improvement_summary = f"Ramachandran: {before_rama}→{after_rama}, RSCC: {before_rscc:.2f}→{after_rscc:.2f}"
        else:
            perf.improvement_summary = f"No improvement (Rama: {after_rama}, RSCC: {after_rscc:.2f})"

        return CPluginScript.SUCCEEDED
```

### Agent Interpretation of Local Fix Results

The result interpreter for local fixes provides actionable feedback:

```python
# server/ccp4i2/agent/interpreters/local_fix.py

def interpret_local_fix_result(job_id: str) -> dict:
    """
    Interpret results from a Coot local fix job.

    Returns structured assessment the agent can reason about.
    """
    job = get_job(job_id)
    perf = job.performance  # CLocalFixPerformance

    result = {
        "job_id": str(job_id),
        "task": job.task_name,
        "residue": perf.residue_id,
        "success": perf.success,

        "changes": {
            "ramachandran": {
                "before": perf.ramachandran_before,
                "after": perf.ramachandran_after,
                "improved": is_ramachandran_improvement(
                    perf.ramachandran_before,
                    perf.ramachandran_after
                )
            },
            "density_fit": {
                "rscc_before": perf.rscc_before,
                "rscc_after": perf.rscc_after,
                "improved": perf.rscc_after > perf.rscc_before + 0.02
            }
        },

        "summary": perf.improvement_summary,
    }

    # Generate recommendations based on outcome
    if not perf.success:
        result["recommendations"] = generate_fix_failure_recommendations(job, perf)
    else:
        result["recommendations"] = ["proceed_to_refinement"]

    return result


def generate_fix_failure_recommendations(job, perf) -> list:
    """
    When a local fix didn't help, suggest alternatives.
    """
    recommendations = []

    if job.task_name == "coot_flip_peptide":
        if perf.rscc_after < 0.7:
            recommendations.append({
                "action": "coot_real_space_refine_zone",
                "reason": "Density fit still poor after flip - try RSR",
                "params": {"chain": perf.residue_id.split("/")[0],
                          "start": int(perf.residue_id.split("/")[1]) - 2,
                          "end": int(perf.residue_id.split("/")[1]) + 2}
            })
        recommendations.append({
            "action": "undo_and_skip",
            "reason": "Flip didn't improve - may be genuine outlier"
        })

    elif job.task_name == "coot_auto_fit_rotamer":
        if perf.rscc_after < 0.5:
            recommendations.append({
                "action": "coot_delete_sidechain",
                "reason": "No good rotamer fits density - sidechain may be disordered"
            })

    return recommendations
```

### Example Agent Reasoning with Metrics

```
Agent: Attempting to fix Ramachandran outlier at A/45...

[Calls: coot_flip_peptide(chain="A", residue=45)]

Result:
  - Ramachandran: OUTLIER → FAVORED ✓
  - RSCC: 0.72 → 0.68 (slight decrease)
  - Success: True (Ramachandran improved)

The peptide flip fixed the Ramachandran outlier. The slight RSCC decrease
is expected - the atoms moved. I'll include this residue in a local
real-space refinement zone to optimize density fit.

[Calls: coot_real_space_refine_zone(chain="A", start=43, end=47)]

Result:
  - Zone RSCC: 0.68 → 0.81 ✓
  - Geometry score: 0.92 → 0.95 ✓

Excellent - the zone now has good density fit and geometry.

---

[Later...]

Agent: Attempting to fix rotamer outlier at B/123 (Leu)...

[Calls: coot_auto_fit_rotamer(chain="B", residue=123)]

Result:
  - Rotamer: OUTLIER → OUTLIER (unchanged)
  - RSCC: 0.45 → 0.42 (no improvement)
  - Success: False

The auto-fit couldn't find a good rotamer. With RSCC < 0.5, this sidechain
likely has no ordered density. I'll truncate it to alanine.

[Calls: coot_delete_sidechain(chain="B", residue=123)]

Result:
  - Atoms removed: 4 (Leu sidechain beyond CB)
  - This will reduce overfitting and should help R-free
```

### Required Metrics Infrastructure

To support this, each Coot wrapper needs access to:

1. **Map coefficients** - To calculate RSCC/RSR (passed as input from previous refinement)
2. **Validation libraries** - Ramachandran/rotamer analysis (gemmi or Coot's own)
3. **Pre-fix state** - Stored in `processInputFiles()` before modification

```python
def processInputFiles(self):
    """Store pre-fix quality metrics before Coot modifies the model."""

    chain = self.container.controlParameters.CHAIN
    resno = self.container.controlParameters.RESIDUE

    # Calculate and store "before" metrics
    self.stored_before_ramachandran = calculate_ramachandran(
        self.container.inputData.XYZIN.fullPath,
        chain, resno
    )
    self.stored_before_rscc = calculate_rscc(
        self.container.inputData.XYZIN.fullPath,
        self.container.inputData.FPHIOUT.fullPath,
        chain, resno
    )

    return CPluginScript.SUCCEEDED
```

This creates a complete feedback loop: **attempt fix → measure outcome → decide next action**.

### Key Insights for Implementation

1. **Unified Job System**: Since Coot tools are CCP4i2 jobs, model handoff is automatic via `set_context_job`. The agent just chains jobs: `refmac → coot_flip_peptide → refmac`. No special integration needed.

2. **Validation as Trigger**: Model validation (MolProbity, Coot's own checks) becomes the key decision point - it tells the agent when to switch from global to local mode.

3. **WHATNEXT Encoding**: Each Coot wrapper's `WHATNEXT` attribute guides the agent. For example:
   - `coot_flip_peptide.WHATNEXT = ['refmac', 'coot_real_space_refine']`
   - `coot_add_waters.WHATNEXT = ['refmac']`
   - `refmac.WHATNEXT = ['coot_flip_peptide', 'coot_auto_fit_rotamer', 'modelcraft']`

4. **Incremental Fixes**: Fix one category of problem, then re-refine. Don't try to fix everything at once - some problems may resolve themselves after refinement.

5. **Confidence Levels**: Some Coot operations (flip peptide, add water) are low-risk and can be automated. Others (loop building, major rebuilds) may need human confirmation - this can be encoded in wrapper metadata.

6. **State Tracking**: The job database already tracks history. The agent can query previous jobs to avoid repeating failed fixes.

---

## Current Infrastructure Analysis

### Plugin System

The plugin registry (`core/task_manager/plugin_registry.py`) provides lazy-loading access to all crystallographic tasks:

```python
from ccp4i2.core.task_manager.plugin_registry import get_registry

registry = get_registry()
plugins = registry.list_plugins()        # Returns ~176 plugin names
metadata = registry.get_plugin_metadata('refmac')  # Task metadata
plugin_class = registry.get_plugin_class('refmac') # Lazy-loaded class
```

Each plugin defines rich metadata:

| Attribute | Purpose | Example |
|-----------|---------|---------|
| `TASKNAME` | Unique identifier | `refmac` |
| `TASKTITLE` | Human-readable name | `Refinement (Refmac5)` |
| `TASKMODULE` | Category | `wrappers` |
| `WHATNEXT` | Suggested follow-up tasks | `['prosmart_refmac', 'modelcraft']` |
| `PERFORMANCECLASS` | Metrics class | `CRefinementPerformance` |
| `ERROR_CODES` | Task-specific errors | `{201: {'description': '...'}}` |

### Performance Metrics

Structured classes capture key crystallographic quality indicators:

#### CRefinementPerformance
```python
RFactor: float      # Working R-factor
RFree: float        # Free R-factor
RMSBond: float      # Bond length RMSD (target ~0.02 Å)
RMSAngle: float     # Bond angle RMSD (target ~2.0°)
weightUsed: float   # Geometry weight used
```

#### CDataReductionPerformance
```python
spaceGroup: str     # Space group assignment
highResLimit: float # High resolution limit
rMeas: float        # Redundancy-independent R-factor
```

#### CModelBuildPerformance
```python
RFactor: float      # R-factor after building
completeness: float # Model completeness (%)
```

#### CExpPhasPerformance
```python
FOM: float          # Figure of merit
CFOM: float         # Combined FOM
Hand1Score: float   # Score for hand 1
Hand2Score: float   # Score for hand 2
CC: float           # Map correlation
```

#### CPhaseErrorPerformance
```python
phaseError: float           # Mean phase error
weightedPhaseError: float   # Weighted phase error
reflectionCorrelation: float # Map-model correlation
```

### REST API Endpoints

The `JobViewSet` provides all necessary operations:

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/api/jobs/` | POST | Create new job |
| `/api/jobs/{id}/run/` | POST | Execute job |
| `/api/jobs/{id}/container/` | GET | Get full job state as JSON |
| `/api/jobs/{id}/report_xml/` | GET | Get results report |
| `/api/jobs/{id}/validation/` | GET | Validate parameters |
| `/api/jobs/{id}/what_next/` | GET | Get suggested next steps |
| `/api/jobs/{id}/set_parameter/` | POST | Set job parameter |
| `/api/jobs/{id}/set_context_job/` | POST | Chain jobs together |
| `/api/jobs/{id}/diagnostic_xml/` | GET | Get error diagnostics |

### Task Definition Files

Each task has a `.def.xml` file defining inputs, outputs, and parameters with validation rules:

```xml
<container id="inputData">
  <content id="XYZIN">
    <className>CPdbDataFile</className>
    <qualifiers>
      <mustExist>True</mustExist>
      <toolTip>File containing model coordinates</toolTip>
    </qualifiers>
  </content>
</container>

<container id="controlParameters">
  <content id="NCYCLES">
    <className>CInt</className>
    <qualifiers>
      <default>10</default>
      <min>1</min>
      <max>100</max>
    </qualifiers>
  </content>
</container>
```

---

## Proposed Agent Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     Claude Agent                                 │
│                                                                  │
│  Capabilities:                                                   │
│  - Understands crystallographic concepts                         │
│  - Plans multi-step workflows                                    │
│  - Interprets quality metrics                                    │
│  - Diagnoses problems and suggests fixes                         │
└─────────────────────┬───────────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────────┐
│              Agent Tool Layer                                    │
│                                                                  │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐  │
│  │  Tool Schemas   │  │  Interpreters   │  │ Knowledge Base  │  │
│  │                 │  │                 │  │                 │  │
│  │ - run_refine    │  │ - refinement    │  │ - decision      │  │
│  │ - run_phaser    │  │ - phasing       │  │   rules         │  │
│  │ - run_build     │  │ - validation    │  │ - thresholds    │  │
│  │ - get_metrics   │  │ - model stats   │  │ - strategies    │  │
│  └─────────────────┘  └─────────────────┘  └─────────────────┘  │
└─────────────────────┬───────────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────────┐
│           CCP4i2 Infrastructure                                  │
│                                                                  │
│  Django REST API → Plugin System → CCP4 Programs                │
└─────────────────────────────────────────────────────────────────┘
```

---

## Components to Build

### 1. Tool Schema Generator

**Purpose**: Auto-generate Claude tool definitions from plugin metadata.

**Input**: `plugin_lookup.json` + `.def.xml` files

**Output**: Claude-compatible tool schemas

```python
# Example generated tool schema
{
    "name": "run_refinement",
    "description": """Run REFMAC5 crystallographic refinement.

    Use when:
    - Model needs refinement against diffraction data
    - R-factors are high and need improvement
    - After molecular replacement or model building

    Interprets: R-factor, R-free, geometry statistics""",

    "input_schema": {
        "type": "object",
        "properties": {
            "project_id": {
                "type": "string",
                "description": "Project UUID"
            },
            "model_job_id": {
                "type": "string",
                "description": "Job ID containing input model (XYZOUT)"
            },
            "data_job_id": {
                "type": "string",
                "description": "Job ID containing reflection data (F_SIGF)"
            },
            "num_cycles": {
                "type": "integer",
                "default": 10,
                "minimum": 1,
                "maximum": 100,
                "description": "Number of refinement cycles"
            },
            "use_tls": {
                "type": "boolean",
                "default": false,
                "description": "Use TLS refinement for flexible regions"
            }
        },
        "required": ["project_id", "model_job_id", "data_job_id"]
    }
}
```

**Implementation location**: `server/ccp4i2/agent/tool_generator.py`

### 2. Result Interpreters

**Purpose**: Transform raw job outputs into agent-friendly summaries with diagnostic insights.

```python
# server/ccp4i2/agent/interpreters/refinement.py

def interpret_refinement_result(job_id: str) -> dict:
    """
    Convert REFMAC output to agent-interpretable metrics.

    Returns structured analysis that Claude can reason about.
    """
    job = Job.objects.get(uuid=job_id)
    report = get_job_report_xml(job)

    # Extract metrics from XML
    initial_r = extract_initial_r_factor(report)
    final_r = extract_final_r_factor(report)
    initial_rfree = extract_initial_r_free(report)
    final_rfree = extract_final_r_free(report)
    rms_bond = extract_rms_bond(report)
    rms_angle = extract_rms_angle(report)

    # Calculate derived metrics
    r_free_gap = final_rfree - final_r
    r_improvement = initial_r - final_r

    # Generate diagnosis
    warnings = []
    diagnosis_parts = []

    if r_free_gap > 0.07:
        warnings.append("Large R-free gap suggests possible overfitting")
        diagnosis_parts.append("Consider reducing refinement parameters")

    if rms_bond > 0.025:
        warnings.append("Bond length deviations are high")
        diagnosis_parts.append("May need to increase geometry weight")

    if r_improvement < 0.001:
        warnings.append("R-factors did not improve significantly")
        diagnosis_parts.append("Model may need rebuilding rather than more refinement")

    return {
        "success": True,
        "job_id": str(job_id),
        "task": "refmac",

        "summary": f"Refinement {'improved' if r_improvement > 0.01 else 'converged'}: "
                   f"R={final_r:.3f} R-free={final_rfree:.3f}",

        "metrics": {
            "r_factor": {
                "initial": initial_r,
                "final": final_r,
                "improved": r_improvement > 0.005,
                "change": -r_improvement
            },
            "r_free": {
                "initial": initial_rfree,
                "final": final_rfree,
                "improved": (initial_rfree - final_rfree) > 0.005
            },
            "r_free_gap": {
                "value": r_free_gap,
                "status": "good" if r_free_gap < 0.05 else
                         "acceptable" if r_free_gap < 0.07 else "concerning"
            },
            "geometry": {
                "rms_bond": rms_bond,
                "rms_angle": rms_angle,
                "bond_status": "good" if rms_bond < 0.02 else "high",
                "angle_status": "good" if rms_angle < 2.5 else "high"
            }
        },

        "warnings": warnings,
        "diagnosis": " ".join(diagnosis_parts) if diagnosis_parts else
                    "Refinement progressed normally",

        "suggested_actions": _suggest_refinement_actions(
            r_free_gap, r_improvement, rms_bond, rms_angle
        ),

        "output_files": {
            "model": str(job.directory / "XYZOUT.pdb"),
            "map_coefficients": str(job.directory / "FPHIOUT.mtz"),
            "difference_map": str(job.directory / "DIFFPHIOUT.mtz")
        }
    }


def _suggest_refinement_actions(r_free_gap, r_improvement, rms_bond, rms_angle):
    """Generate contextual suggestions for next steps."""
    actions = []

    if r_improvement < 0.001:
        actions.append({
            "action": "run_model_building",
            "reason": "Refinement has converged; try automated model building",
            "priority": "high"
        })

    if r_free_gap > 0.06:
        actions.append({
            "action": "run_validation",
            "reason": "Check for model errors causing overfitting",
            "priority": "high"
        })

    actions.append({
        "action": "add_waters",
        "reason": "Adding ordered waters may improve R-factors",
        "priority": "medium"
    })

    return actions
```

**Interpreters needed**:
- `refinement.py` - REFMAC, Servalcat results
- `phasing.py` - Phaser MR, experimental phasing
- `model_building.py` - ModelCraft, Buccaneer results
- `validation.py` - MolProbity, geometry checks
- `data_processing.py` - Aimless, scaling results

### 3. Crystallographic Knowledge Base

**Purpose**: Encode domain expertise as structured rules for decision-making.

```python
# server/ccp4i2/agent/knowledge/refinement_rules.py

REFINEMENT_DECISION_RULES = {
    "high_r_free_gap": {
        "condition": lambda m: m["r_free_gap"]["value"] > 0.07,
        "severity": "warning",
        "diagnosis": "Large R-free gap (>{gap:.2f}) suggests overfitting",
        "possible_causes": [
            "Too many refinement parameters",
            "Incorrect test set selection",
            "Model errors being fitted"
        ],
        "recommended_actions": [
            "reduce_parameters",
            "check_test_set",
            "run_validation"
        ]
    },

    "stalled_refinement": {
        "condition": lambda m: abs(m["r_factor"]["change"]) < 0.002,
        "severity": "info",
        "diagnosis": "Refinement has converged or stalled",
        "possible_causes": [
            "Model is well refined",
            "Model has errors preventing improvement",
            "Resolution limit reached"
        ],
        "recommended_actions": [
            "run_model_building",
            "inspect_density",
            "try_different_weight"
        ]
    },

    "poor_geometry": {
        "condition": lambda m: m["geometry"]["rms_bond"] > 0.025,
        "severity": "warning",
        "diagnosis": "Bond lengths deviate significantly from ideal",
        "possible_causes": [
            "Geometry weight too low",
            "Missing or incorrect restraints",
            "Genuine structural distortion"
        ],
        "recommended_actions": [
            "increase_geometry_weight",
            "check_ligand_restraints"
        ]
    },

    "excellent_refinement": {
        "condition": lambda m: (
            m["r_factor"]["final"] < 0.20 and
            m["r_free_gap"]["value"] < 0.04 and
            m["geometry"]["rms_bond"] < 0.015
        ),
        "severity": "success",
        "diagnosis": "Refinement statistics are excellent",
        "recommended_actions": [
            "final_validation",
            "prepare_deposition"
        ]
    }
}


# Phasing strategy selection
PHASING_STRATEGY = {
    "molecular_replacement": {
        "requirements": [
            "Search model with >25% sequence identity",
            "Reasonable resolution (<4Å preferred)"
        ],
        "indicators_of_success": {
            "llg": ">60 strongly suggests correct solution",
            "tfz": ">8 indicates likely correct placement",
            "r_factor_after_refinement": "<0.45 promising"
        }
    },

    "sad_phasing": {
        "requirements": [
            "Anomalous signal present",
            "High redundancy data",
            "Wavelength near absorption edge"
        ],
        "indicators_of_success": {
            "anomalous_signal": "Measurable to at least 4Å",
            "sites_found": "Heavy atom substructure located",
            "fom": ">0.3 for initial phases"
        }
    }
}
```

### 4. Structure Analysis Tools

**Purpose**: Tools for inspecting models and maps that an agent can call.

```python
# server/ccp4i2/agent/tools/analysis.py

def analyze_model_geometry(model_path: str) -> dict:
    """
    Analyze model geometry using gemmi.

    Returns Ramachandran statistics, rotamer outliers, clashes.
    """
    import gemmi

    structure = gemmi.read_structure(model_path)

    # Ramachandran analysis
    rama_favored = 0
    rama_allowed = 0
    rama_outliers = []

    for model in structure:
        for chain in model:
            for residue in chain:
                # Calculate phi/psi and classify
                phi, psi = calculate_phi_psi(residue)
                classification = classify_ramachandran(phi, psi, residue.name)

                if classification == "favored":
                    rama_favored += 1
                elif classification == "allowed":
                    rama_allowed += 1
                else:
                    rama_outliers.append({
                        "chain": chain.name,
                        "residue": residue.seqid.num,
                        "name": residue.name,
                        "phi": phi,
                        "psi": psi
                    })

    total = rama_favored + rama_allowed + len(rama_outliers)

    return {
        "ramachandran": {
            "favored_percent": 100 * rama_favored / total,
            "allowed_percent": 100 * rama_allowed / total,
            "outlier_percent": 100 * len(rama_outliers) / total,
            "outliers": rama_outliers[:10],  # Top 10 outliers
            "status": "good" if len(rama_outliers) / total < 0.02 else "needs_attention"
        },
        "summary": f"Ramachandran: {100*rama_favored/total:.1f}% favored, "
                   f"{len(rama_outliers)} outliers"
    }


def find_unmodeled_density(map_path: str, model_path: str) -> dict:
    """
    Find significant density peaks not explained by the model.

    Useful for identifying missing ligands, waters, or alternate conformations.
    """
    import gemmi

    # Load map and model
    mtz = gemmi.read_mtz_file(map_path)
    structure = gemmi.read_structure(model_path)

    # Calculate difference map peaks
    peaks = find_difference_peaks(mtz, structure, sigma_cutoff=3.5)

    return {
        "positive_peaks": [
            {
                "position": [p.x, p.y, p.z],
                "sigma": p.height,
                "nearby_residue": find_nearest_residue(p, structure),
                "interpretation": guess_peak_identity(p)
            }
            for p in peaks if p.height > 0
        ][:20],
        "negative_peaks": [
            {
                "position": [p.x, p.y, p.z],
                "sigma": abs(p.height),
                "nearby_residue": find_nearest_residue(p, structure)
            }
            for p in peaks if p.height < 0
        ][:10],
        "summary": f"Found {len([p for p in peaks if p.height > 0])} positive peaks "
                   f"and {len([p for p in peaks if p.height < 0])} negative peaks"
    }
```

### 5. Session State Manager

**Purpose**: Track the structure solution progress across multiple interactions.

```python
# server/ccp4i2/agent/session.py

from dataclasses import dataclass, field
from typing import List, Optional
from datetime import datetime


@dataclass
class AttemptRecord:
    """Record of a single task attempt."""
    timestamp: datetime
    task_name: str
    job_id: str
    parameters: dict
    result_summary: str
    metrics: dict
    success: bool
    notes: str = ""


@dataclass
class StructureSolutionSession:
    """
    Tracks the complete structure solution journey.

    Provides context for the agent to make informed decisions.
    """
    session_id: str
    project_id: str
    created_at: datetime = field(default_factory=datetime.now)

    # Current best results
    best_model_job_id: Optional[str] = None
    best_r_factor: Optional[float] = None
    best_r_free: Optional[float] = None

    # History
    attempts: List[AttemptRecord] = field(default_factory=list)

    # Problem tracking
    known_issues: List[dict] = field(default_factory=list)
    resolved_issues: List[dict] = field(default_factory=list)

    def record_attempt(self, task_name: str, job_id: str,
                       parameters: dict, result: dict):
        """Record a task attempt and update best results if improved."""

        record = AttemptRecord(
            timestamp=datetime.now(),
            task_name=task_name,
            job_id=job_id,
            parameters=parameters,
            result_summary=result.get("summary", ""),
            metrics=result.get("metrics", {}),
            success=result.get("success", False)
        )
        self.attempts.append(record)

        # Update best if this is a refinement with better R-free
        if task_name in ["refmac", "servalcat"] and result.get("success"):
            r_free = result.get("metrics", {}).get("r_free", {}).get("final")
            if r_free and (self.best_r_free is None or r_free < self.best_r_free):
                self.best_model_job_id = job_id
                self.best_r_free = r_free
                self.best_r_factor = result["metrics"]["r_factor"]["final"]

    def get_context_summary(self) -> str:
        """
        Generate a summary of the session for the agent.

        Used to provide context at the start of each interaction.
        """
        lines = [
            f"## Structure Solution Session",
            f"Project: {self.project_id}",
            f"Started: {self.created_at.strftime('%Y-%m-%d %H:%M')}",
            f"Attempts: {len(self.attempts)}",
            ""
        ]

        if self.best_r_free:
            lines.extend([
                "### Current Best Model",
                f"- R-factor: {self.best_r_factor:.4f}",
                f"- R-free: {self.best_r_free:.4f}",
                f"- From job: {self.best_model_job_id}",
                ""
            ])

        if self.known_issues:
            lines.append("### Known Issues")
            for issue in self.known_issues:
                lines.append(f"- {issue['description']}")
            lines.append("")

        # Recent history
        if self.attempts:
            lines.append("### Recent Attempts")
            for attempt in self.attempts[-5:]:
                status = "✓" if attempt.success else "✗"
                lines.append(
                    f"- {status} {attempt.task_name}: {attempt.result_summary}"
                )

        return "\n".join(lines)

    def what_hasnt_been_tried(self) -> List[str]:
        """Suggest tasks that haven't been attempted yet."""
        tried_tasks = {a.task_name for a in self.attempts}

        suggestions = []
        if "modelcraft" not in tried_tasks:
            suggestions.append("Automated model building with ModelCraft")
        if "molprobity" not in tried_tasks:
            suggestions.append("Geometry validation with MolProbity")
        if not any("water" in t for t in tried_tasks):
            suggestions.append("Add ordered water molecules")

        return suggestions
```

---

## Tool Definitions for Claude

### Core Crystallographic Tools

```python
CRYSTALLOGRAPHIC_TOOLS = [
    {
        "name": "create_and_run_refinement",
        "description": """Run crystallographic refinement with REFMAC5.

        Refines atomic coordinates and B-factors against diffraction data.

        USE WHEN:
        - After molecular replacement to improve initial model
        - After model building to optimize geometry and fit
        - R-factors need improvement

        INTERPRETS:
        - R-factor and R-free (lower is better, gap should be <0.05)
        - Geometry: RMS bonds (~0.02Å) and angles (~2.0°)
        - Convergence: whether refinement is still improving
        """,
        "input_schema": {
            "type": "object",
            "properties": {
                "project_id": {"type": "string"},
                "context_job_id": {
                    "type": "string",
                    "description": "Job to take input files from"
                },
                "num_cycles": {
                    "type": "integer",
                    "default": 10,
                    "description": "Number of refinement cycles"
                },
                "weight": {
                    "type": "number",
                    "description": "Geometry weight (auto if not specified)"
                }
            },
            "required": ["project_id", "context_job_id"]
        }
    },

    {
        "name": "create_and_run_model_building",
        "description": """Run automated model building with ModelCraft.

        Automatically builds protein model into electron density.

        USE WHEN:
        - After phasing to build initial model
        - After refinement stalls to improve/extend model
        - Model completeness is low

        REQUIRES:
        - Reflection data with phases or map coefficients
        - Sequence file (FASTA format)
        """,
        "input_schema": {
            "type": "object",
            "properties": {
                "project_id": {"type": "string"},
                "context_job_id": {"type": "string"},
                "sequence_file": {
                    "type": "string",
                    "description": "Path or file ID for sequence"
                }
            },
            "required": ["project_id", "context_job_id"]
        }
    },

    {
        "name": "create_and_run_molecular_replacement",
        "description": """Run molecular replacement with Phaser.

        Places a search model in the unit cell by matching Patterson functions.

        USE WHEN:
        - Starting structure solution with a homologous model
        - Search model has >25% sequence identity

        SUCCESS INDICATORS:
        - LLG (Log-Likelihood Gain) > 60 is strong
        - TFZ (Translation Function Z-score) > 8 is good
        - R-factor after rigid body < 0.50
        """,
        "input_schema": {
            "type": "object",
            "properties": {
                "project_id": {"type": "string"},
                "data_job_id": {"type": "string"},
                "search_model_path": {"type": "string"},
                "sequence_identity": {
                    "type": "number",
                    "description": "Expected sequence identity (0-1)"
                },
                "copies": {
                    "type": "integer",
                    "default": 1,
                    "description": "Number of copies in asymmetric unit"
                }
            },
            "required": ["project_id", "data_job_id", "search_model_path"]
        }
    },

    {
        "name": "get_job_metrics",
        "description": """Get interpreted metrics from a completed job.

        Returns structured analysis with diagnosis and suggestions.
        """,
        "input_schema": {
            "type": "object",
            "properties": {
                "job_id": {"type": "string"}
            },
            "required": ["job_id"]
        }
    },

    {
        "name": "analyze_model",
        "description": """Analyze a model for geometry issues and problem regions.

        Returns Ramachandran statistics, clashes, rotamer outliers.
        """,
        "input_schema": {
            "type": "object",
            "properties": {
                "job_id": {
                    "type": "string",
                    "description": "Job containing the model to analyze"
                }
            },
            "required": ["job_id"]
        }
    },

    {
        "name": "find_density_features",
        "description": """Find unmodeled density peaks and negative density.

        Useful for finding missing ligands, waters, or model errors.
        """,
        "input_schema": {
            "type": "object",
            "properties": {
                "job_id": {"type": "string"},
                "sigma_cutoff": {
                    "type": "number",
                    "default": 3.5,
                    "description": "Sigma cutoff for peak detection"
                }
            },
            "required": ["job_id"]
        }
    },

    {
        "name": "get_session_context",
        "description": """Get summary of current structure solution progress.

        Returns history of attempts, current best model, known issues.
        """,
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"}
            },
            "required": ["session_id"]
        }
    },

    {
        "name": "get_suggested_next_steps",
        "description": """Get AI-informed suggestions for what to try next.

        Analyzes current state and recommends actions with reasoning.
        """,
        "input_schema": {
            "type": "object",
            "properties": {
                "job_id": {
                    "type": "string",
                    "description": "Most recent job to base suggestions on"
                }
            },
            "required": ["job_id"]
        }
    }
]
```

---

## Example Agent Interaction

```
User: I have diffraction data and a search model. Help me solve this structure.

Agent: I'll help you solve this structure using molecular replacement. Let me first
understand what data you have.

[Calls: get_project_files(project_id)]

I can see you have:
- Reflection data: data.mtz (2.1Å resolution)
- Search model: search_model.pdb (45% sequence identity)

This looks like a good candidate for molecular replacement. The 45% sequence
identity should give a strong signal. Let me run Phaser.

[Calls: create_and_run_molecular_replacement(
    project_id="...",
    data_job_id="...",
    search_model_path="search_model.pdb",
    sequence_identity=0.45
)]

Molecular replacement completed successfully:
- LLG: 847 (excellent - well above 60 threshold)
- TFZ: 24.3 (very strong - above 8 threshold)
- 1 copy placed in the asymmetric unit

The solution looks correct. Now I'll refine the placed model to improve the fit.

[Calls: create_and_run_refinement(
    project_id="...",
    context_job_id="mr_job_123",
    num_cycles=20
)]

Initial refinement results:
- R-factor: 0.312 → 0.247 (improved by 0.065)
- R-free: 0.358 → 0.298 (improved by 0.060)
- R-free gap: 0.051 (acceptable)
- Geometry: bonds 0.014Å, angles 1.7° (good)

The refinement improved significantly, but we can do better. The R-factors are
still relatively high. I recommend automated model building to improve and
complete the model.

[Calls: create_and_run_model_building(
    project_id="...",
    context_job_id="refine_job_456"
)]

Model building completed:
- Completeness: 94% of residues built
- R-factor improved to 0.218
- R-free improved to 0.261

The model is nearly complete. Let me run a final round of refinement and then
validation to identify any remaining issues.

[Continues with refinement and validation...]
```

---

## Implementation Roadmap

### Phase 1: Core Tool Interface (2-3 weeks)
1. Tool schema generator from plugin metadata
2. Basic job creation and execution tools
3. Result interpreter for refinement tasks

### Phase 2: Expanded Tools (2-3 weeks)
1. Interpreters for phasing, model building, validation
2. Model analysis tools (geometry, density)
3. Session state management

### Phase 3: Knowledge Integration (2-3 weeks)
1. Decision rules for common scenarios
2. Diagnostic reasoning
3. Strategy selection logic

### Phase 4: Testing & Refinement (ongoing)
1. Test with real structure solution cases
2. Refine prompts and knowledge base
3. Add handling for edge cases

---

## File Structure

```
server/ccp4i2/agent/
├── __init__.py
├── tools/
│   ├── __init__.py
│   ├── generator.py          # Tool schema generator
│   ├── job_tools.py          # Job creation/execution
│   └── analysis_tools.py     # Model/map analysis
├── interpreters/
│   ├── __init__.py
│   ├── refinement.py
│   ├── phasing.py
│   ├── model_building.py
│   └── validation.py
├── knowledge/
│   ├── __init__.py
│   ├── refinement_rules.py
│   ├── phasing_strategy.py
│   └── thresholds.py
├── session.py                 # Session state manager
└── api.py                     # Agent API endpoints
```

---

## References

- [CCP4i2 JobViewSet](../server/ccp4i2/api/JobViewSet.py) - REST API implementation
- [Plugin Registry](../server/ccp4i2/core/task_manager/plugin_registry.py) - Task discovery
- [Performance Classes](../server/ccp4i2/core/cdata_stubs/CCP4PerformanceData.py) - Metrics definitions
- [REFMAC Wrapper](../server/ccp4i2/wrappers/refmac_i2/script/refmac_i2.py) - Example task implementation
