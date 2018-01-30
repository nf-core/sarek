# Release checklist
This checklist is for our own reference

1. Check that everything is up to date and ready to go
2. Increase version numbers.
3. Update version numbers in code: `main.nf`, `buildContainers.nf`, `buildReferences.nf`
4. If any changes on any containers, match the tag to current version `docker.config`, `singularity.config`, `singularity-path.config`.
5. Build, and get the containers.
  - `./scripts/do_all.sh --push`
  - `./scripts/do_all.sh --pull`
6. Test against sample data.
  - Check for any command line errors
  - Check version numbers are printed correctly
  - `./scripts/test.sh -p docker`
  - `./scripts/test.sh -p singularity`
  - `./scripts/test.sh -p singularityPath`
7. Commit and push version updates
8. Make a [release](https://github.com/SciLifeLab/CAW/releases) on GitHub - list PRs as changelog.
9. Tweet that new version is released
10. Commit and push. Continue making more awesome :metal:
11. Have fika :cake:
