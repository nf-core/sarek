# Session Management Guide

## Session Brief Template

Before starting a session, create `TODO-rewrite-session-{N}.md` with:

```markdown
# Session {N} — {Date}

## Goal
{One sentence: what this session will accomplish}

## Context
- Plan phase: {which phase from TODO-rewrite-plan.md}
- Previous session: {link to TODO-rewrite-session-{N-1}.md or "first session"}
- Blockers from last session: {any known issues}

## Files to Touch
- {list of files that will be created or modified}

## Acceptance Criteria
- {specific, testable criteria for this session being "done"}
```

## Session Log Template

At the end of a session, update the session file with:

```markdown
## What Was Done
- {bullet list of accomplishments}

## What Works
- {bullet list of tests/validations that pass}

## What Doesn't Work Yet
- {bullet list of remaining issues, with specifics}

## Discoveries
- {anything learned about the upstream tool, edge cases, etc.}

## Next Session
- {what the next session should focus on}
- Estimated effort: {low/medium/high}
- Recommended model: {sonnet/opus}

## Git Commits
- {commit hash}: {message}
```

## Estimating Remaining Work

After each session, update the plan with revised estimates:

| Phase | Status | Sessions Est. | Sessions Used | Model |
|-------|--------|--------------|---------------|-------|
| 1     | Done   | 2            | 2             | sonnet |
| 2     | In progress | 3       | 1             | sonnet |
| 3     | Not started | 2      | 0             | opus   |

## When to Use Opus vs Sonnet

**Use Sonnet (default) for:**
- Writing code from a clear specification
- Fixing bugs where the cause is known
- Writing tests, validation scripts, documentation
- Mechanical refactoring

**Escalate to Opus for:**
- Reverse-engineering upstream tool behavior from source code
- Diagnosing discrepancies with no clear cause
- Architectural decisions (module structure, data flow)
- Situations where you've been stuck for > 30 minutes
