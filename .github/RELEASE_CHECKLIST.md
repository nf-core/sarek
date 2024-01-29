# Release checklist

> This checklist is for our own reference, to help us prepare a new release

1. Check that everything is ready to go

-   Desired [PRs](https://github.com/nf-core/sarek/pulls) are merged
-   [GHA CI](https://github.com/nf-core/sarek/actions/workflows/ci.yml) are passing on `dev`
-   [GHA linting](https://github.com/nf-core/sarek/actions/workflows/linting.yml) are passing on `dev`

2. Increase version number following [semantic versioning](http://semver.org/spec/v2.0.0.html)
3. Choose an appropriate codename for the release (if major or minor)

-   i.e. Peaks in [Sarek National Park](https://en.wikipedia.org/wiki/Sarek_National_Park#Topography)

4. Sync `dev` and checkout a new branch for the release
5. Bump version:

-   `nf-core bump-version . <VERSION>`
-   edit `CHANGELOG`

6. Make a PR to `master`
7. Wait for reviews
8. Merge said PR
9. Make a [release](https://github.com/nf-core/sarek/releases) on GitHub
10. Update [bio.tools](https://bio.tools/nf-core-sarek) with the new release details
11. RT the nf-core automated tweet about the new released version
12. Make a new branch from `dev`
13. Checkout the `CHANGELOG.md` from `master`

-   `git checkout upstream/master -- CHANGELOG.md`

14. Add a new `Unreleased` section in `CHANGELOG.md` for the `dev` version
15. Checkout figures from `master`

-   `git checkout upstream/master -- docs/images/sarek_indices_subway.svg`
-   `git checkout upstream/master -- docs/images/sarek_indices_subway.png`
-   `git checkout upstream/master -- docs/images/sarek_subway.svg`
-   `git checkout upstream/master -- docs/images/sarek_subway.png`
-   `git checkout upstream/master -- docs/images/sarek_workflow.svg`
-   `git checkout upstream/master -- docs/images/sarek_workflow.png`

16. Make a PR to `dev`
17. Wait for review
18. Merge said PR
19. Commit and push. Continue making more awesome :metal:
20. Have fika :cake:
