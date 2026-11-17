---
name: issue-graph
description: Read-only context for issues, pull requests, and backlogs. Use to trace related work, check PR counts, prioritize follow-ups, filter saved captures, open exact dashboard views, replay triage queries, or compare snapshots.
---

# issue-graph

Before running operational commands, load and read the guidance bundled with the CLI:

```bash
issue-graph skills get core
```

Follow the returned core guidance. When detailed workflows or references are needed,
load and read them through the CLI:

```bash
issue-graph skills get core --full
```

Discover available guidance with `issue-graph skills list`. Inspect command syntax
with `issue-graph skills --help`, `issue-graph skills list --help`, or
`issue-graph skills get --help`.

If the executable, skills command, core, or referenced assets are unavailable, stop
and report the CLI/skill mismatch and the observed error. Do not fabricate operational
guidance, fall back to remembered instructions, or automatically install or upgrade
anything. Ask for an explicitly authorized setup correction before proceeding.
