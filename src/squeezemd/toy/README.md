# toy/ — experimental scripts

This directory holds **experimental, unmaintained** scripts (trajectory streaming,
clustering prototypes, dissociation analysis). They are **not** part of the
installed `squeezemd` workflow: none are wired into the Snakefile or exposed as
console entry points, and some have known broken imports.

They are kept for reference only. They are intentionally excluded from linting
(`ruff` `extend-exclude` in `pyproject.toml`) and are not covered by tests. Treat
anything here as a scratchpad, not as a supported API.
