---
name: rewrite-finish
description: Final cleanup, commit, push, and project summary. Use as the last step after integration is tested and working. Handles git history cleanup, documentation review, temporary file removal, branch push, PR creation, and writing a project summary with lessons learned.
model: sonnet
user-invocable: true
argument-hint: "[--branch <name>]"
---

# Finish and Deliver

Final cleanup, documentation review, commit, push, and write project summary.
$ARGUMENTS may specify a branch name (default: `feature/rewrite-<tool>`).

## Pre-Finish Checklist

Verify ALL of these before proceeding:

```bash
# Code compiles without warnings
cargo build --release 2>&1 | grep -c "warning"

# Formatting is clean
cargo fmt --check

# Clippy passes
cargo clippy -- -D warnings

# All tests pass
cargo test --release

# Validation against reference outputs passes
# (run whatever validation scripts were created in /rewrite-validate)
```

If any fail, fix them first. Do not proceed with incomplete work.

## Git History Cleanup

1. Review the commit history:
   ```bash
   git log --oneline -20
   ```

2. Ensure commits are logical and well-described. If there are many small fixup
   commits from iterative development, consider squash-merging when creating the PR
   (GitHub's "Squash and merge" option), or use non-interactive squash:
   ```bash
   # Squash all commits since branching from main (non-interactive)
   git reset --soft $(git merge-base HEAD main) && git commit -m "feat: Add rewritten <tool>"
   ```
   Only do this on your own unpushed feature branch.

   > **WARNING**: `git reset --soft` is dangerous if you miscount commits or target the wrong merge-base. Before running, verify with `git log --oneline` that you're only squashing your own commits on an unpushed branch. If in doubt, skip squashing — a clean branch with multiple well-messaged commits is better than a botched squash.

3. **Never force-push to shared branches.** Only clean up your own feature branch.

4. Verify no sensitive data in commits:
   ```bash
   git log --all --diff-filter=A -- "*.env" "*.key" "*.pem" "*credentials*"
   ```

## Documentation Final Check

### README.md
- [ ] Describes what the tool does and which upstream tools it replaces
- [ ] Build instructions are complete and tested
- [ ] Usage examples with real command lines
- [ ] Benchmark results with specific numbers, hardware, and test conditions
- [ ] Numbers in README match numbers in `benchmark/README.md`
- [ ] Known limitations documented
- [ ] License specified

### AGENTS.md
- [ ] Build / lint / test commands documented
- [ ] Project structure documented with module descriptions
- [ ] Code style conventions documented
- [ ] Key dependencies listed with purpose
- [ ] CI pipeline described
- [ ] Notes for future AI-assisted development

### Code Comments
- [ ] No stale TODO/FIXME comments (convert to GitHub issues or remove)
- [ ] No commented-out code blocks
- [ ] No debug print statements (`println!`, `eprintln!`, `dbg!`)
- [ ] All `#[allow(...)]` annotations have justification comments

### Benchmark Documentation
- [ ] `benchmark/README.md` exists with methodology and results
- [ ] Test conditions documented (hardware, OS, input data size)
- [ ] Speedup numbers are specific (e.g., "3.7x faster" not "much faster")
- [ ] Memory usage compared
- [ ] Thread scaling documented if applicable

## Remove Temporary Files

```bash
# Remove any TODO session files that are no longer needed
ls TODO-rewrite-session-*.md 2>/dev/null

# Remove any temporary diagnostic scripts
ls scripts/debug_*.py scripts/tmp_*.sh 2>/dev/null

# Check for large files that shouldn't be committed
find . -size +10M -not -path "./.git/*" -not -path "./target/*"
```

Keep useful TODO files that document lessons learned and bug tracking — these are valuable project history. Only remove scratch/debug files.

## Branch and Push

```bash
# Create feature branch (if not already on one)
BRANCH="${1:-feature/rewrite-<tool>}"
git checkout -b "$BRANCH" 2>/dev/null || git checkout "$BRANCH"

# Stage all changes
git add -A

# Commit with descriptive message
git commit -m "Add rewritten <tool> implementation

Replaces <upstream-tool-1>, <upstream-tool-2> with a single Rust binary.
Key improvements:
- <X>x faster on typical inputs
- Single-pass BAM processing
- <other improvements>

Output is validated against upstream tool references."

# Push
git push -u origin "$BRANCH"
```

## Create Pull Request

```bash
gh pr create \
  --title "feat: Add rewritten <tool> (replaces <upstream-tools>)" \
  --body "$(cat <<'EOF'
## Summary

Adds a Rust reimplementation of <upstream-tools> as a single binary (`<binary>`).

## Performance

| Metric | Upstream | Rewritten | Improvement |
|--------|----------|-----------|-------------|
| Wall time | Xm Ys | Xm Ys | X.Xx faster |
| Peak memory | X GB | X GB | X.Xx less |

## Output Matching

All outputs validated against upstream tool references:
- <file1>: exact match
- <file2>: exact match (within float tolerance)
- <etc>

## Integration

- New process `REWRITTEN_TOOL` in `modules/local/<tool>/`
- Toggle: `--use_rewritten_<tool> true|false` (default: true)
- Original processes preserved, unchanged
- Uses Wave module binaries (no container publication needed)

## Validation

- [ ] All outputs match upstream within defined tolerances
- [ ] Link to validation report: `TODO-rewrite-validate.md`
- [ ] Known acceptable differences documented with evidence

## Benchmark Results

- [ ] Link to benchmark report: `benchmark/README.md`

## Testing

- [x] `cargo test --release` passes
- [x] Pipeline `-profile test` passes with rewritten tool
- [x] Output matches original tool output
- [x] MultiQC parses output correctly
EOF
)"
```

## Write Project Summary

Write `TODO-rewrite-summary.md`:

```markdown
# Rewrite Summary: <tool>

## What Was Done
- Replaced: <list of upstream tools>
- Language: <Python/R/shell> → Rust
- Architecture: <single-pass BAM, bundled analyses, etc.>

## Results
- Speedup: <X.Xx> on <test data description>
- Memory: <comparison>
- Output accuracy: <exact match / within tolerance>

## Cost
- Estimated sessions: <N>
- Model usage: <N Opus sessions, M Sonnet sessions>
- Calendar time: <X days>

## Approach
- <Brief description of implementation strategy>
- <Key decisions made and why>

## Lessons Learned
- <What worked well>
- <What was wasteful>
- <Tips for next time>

Structure the lessons learned section with these categories:
- **Upstream source analysis**: How was the original tool's source obtained and understood? What was undocumented?
- **Algorithmic challenges**: Which algorithms required statement-by-statement porting? What numerical precision issues arose?
- **Validation issues**: What output mismatches were found? What caused them? How long did matching take?
- **Performance insights**: What was the actual bottleneck? What optimizations had the most impact?
- **Integration surprises**: What Nextflow/pipeline issues were unexpected?
- **Estimation accuracy**: How did actual effort compare to estimates? What was underestimated?
- **Plot generation**: If applicable, what charting challenges arose?

## Future Opportunities
- <More tools that could be rewritten>
- <Further optimizations possible>
- <Adjacent improvements identified>
- <Tools in other pipelines with similar patterns>
```

## Final Verification

Run one last end-to-end check:

```bash
# Clean build from scratch
cargo clean && cargo build --release

# Full test suite
cargo test --release

# Pipeline test (if integrated)
cd <pipeline-path>
nextflow run . -profile test -resume --outdir final_check
```

Confirm everything works. Then the rewrite is complete.

## Post-Merge Considerations

- Set up regression testing: ensure the pipeline's CI runs with the rewritten tool enabled so future pipeline changes don't break compatibility
- Consider creating a GitHub release with the compiled binary for the rewritten tool
- Update the pipeline's changelog/version notes to document the new tool option
- If the rewrite covers a subset of the upstream tool's functionality, document what is NOT supported
