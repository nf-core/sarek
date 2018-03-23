# Release checklist
This checklist is for our own reference

1. Check that everything is up to date and ready to go
  - Travis test is passing
  - Manual testing on Bianca is passing
2. Increase version numbers.
3. Update version numbers in code: `configuration/base.config`
4. Build, and get the containers.
  - `./scripts/do_all.sh --push --tag <VERSION>`
  - `./scripts/do_all.sh --pull --tag <VERSION>`
5. Test against sample data.
  - Check for any command line errors
  - Check version numbers are printed correctly
  - `./scripts/test.sh -p docker --tag <VERSION>`
  - `./scripts/test.sh -p singularity --tag <VERSION>`
  - `./scripts/test.sh -p singularityPath --tag <VERSION>`
6. Commit and push version updates
7. Make a [release](https://github.com/SciLifeLab/Sarek/releases) on GitHub
8. Tweet that new version is released
9. Commit and push. Continue making more awesome :metal:
10. Have fika :cake:
