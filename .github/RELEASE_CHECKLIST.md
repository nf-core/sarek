# Release checklist

> This checklist is for our own reference, to help us prepare a new release

1.  Check that everything is ready to go

    -   Desired [PRs](https://github.com/SciLifeLab/Sarek/pulls) are merged
    -   [Travis tests](https://travis-ci.org/SciLifeLab/Sarek/branches) are passing on `dev`

2.  Increase version number following [semantic versioning](http://semver.org/spec/v2.0.0.html)
3.  Choose an appropriate codename for the release (if major or minor)
    -   i.e. Peaks in [Sarek National Park](https://en.wikipedia.org/wiki/Sarek_National_Park#Topography)
4.  Build docker containers.

    -   `./scripts/do_all.sh --tag <VERSION>`

5.  Test against sample data.

    -   `./scripts/test.sh -p docker --tag <VERSION>`
    -   Check for any command line errors

6.  Use script to update version in files:

    -   `./scripts/do_release.sh -r "<VERSION>" -c "<CODENAME>"`

7.  Push latest updates
8.  Make a PR against `dev`
9.  Merge said PR
10.  Make a PR from `dev` to `master`
11.  Merge said PR
12. Make a [release](https://github.com/SciLifeLab/Sarek/releases) on GitHub
13. Update [bio.tools](https://bio.tools/Sarek) with the new release details
14. Tweet that a new version is released
15. Add a new `Unreleased` section in `CHANGELOG.md` for the `dev` version
16. Download all new containers to `/sw/data/uppnex/ToolBox/sarek` on `rackham`
17. Download all new containers to `/btb/containers` on `munin`
18. Update symbolic links to `latest`
19. Commit and push. Continue making more awesome :metal:
20. Have fika :cake:
