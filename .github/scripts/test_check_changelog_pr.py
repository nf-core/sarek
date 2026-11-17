import pathlib
import subprocess
import sys
import tempfile
import unittest


CHECK_CHANGELOG_PR = pathlib.Path(__file__).with_name("check_changelog_pr.py")


class CheckChangelogPrTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.repo = pathlib.Path(self.directory.name)
        self.git("init", "--quiet")
        self.git("config", "user.name", "Changelog test")
        self.git("config", "user.email", "changelog@example.invalid")
        self.changelog = self.repo / "CHANGELOG.md"
        self.changelog.write_text("# Changelog\n\n- [#XXX](https://github.com/nf-core/sarek/pull/XXX) - Existing entry\n")
        self.commit_changelog()
        self.base = self.git("rev-parse", "HEAD").stdout.strip()

    def git(self, *args):
        return subprocess.run(["git", *args], cwd=self.repo, check=True, capture_output=True, text=True)

    def commit_changelog(self):
        self.git("add", "CHANGELOG.md")
        self.git("-c", "core.hooksPath=/dev/null", "commit", "--quiet", "-m", "Write changelog fixture")

    def check_changelog(self, added="", pr_number="2468"):
        self.changelog.write_text(self.changelog.read_text() + added)
        self.commit_changelog()
        return subprocess.run(
            [sys.executable, str(CHECK_CHANGELOG_PR), "--base", self.base, "--pr-number", pr_number],
            cwd=self.repo,
            capture_output=True,
            text=True,
        )

    def test_current_pr_passes_with_unchanged_legacy_placeholder(self):
        result = self.check_changelog("- [#2468](https://github.com/nf-core/sarek/pull/2468) - Add check\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("PR #2468", result.stdout)

    def test_placeholder_regression_fails(self):
        result = self.check_changelog("- [#XXX](https://github.com/nf-core/sarek/pull/XXX) - Add check\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn("replace the placeholder", result.stderr)
        self.assertIn("[#2468](https://github.com/nf-core/sarek/pull/2468)", result.stderr)

    def test_placeholder_fails_even_with_current_pr_present(self):
        result = self.check_changelog(
            "- [#2468](https://github.com/nf-core/sarek/pull/2468) - Add check\n"
            "- [#XXX](https://github.com/nf-core/sarek/pull/XXX) - Another change\n"
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("replace the placeholder", result.stderr)

    def test_wrong_or_mismatched_references_fail(self):
        for reference in (
            "[#1357](https://github.com/nf-core/sarek/pull/1357)",
            "[#2468](https://github.com/nf-core/sarek/pull/1357)",
            "[#1357](https://github.com/nf-core/sarek/pull/2468)",
        ):
            with self.subTest(reference=reference):
                self.git("restore", "--source", self.base, "--", "CHANGELOG.md")
                result = self.check_changelog(f"- {reference} - Add check\n")
                self.assertEqual(result.returncode, 1)
                self.assertIn("add a CHANGELOG.md entry", result.stderr)

    def test_label_and_url_placeholders_fail(self):
        for reference in (
            "[#XXX](https://github.com/nf-core/sarek/pull/2468)",
            "[#2468](https://github.com/nf-core/sarek/pull/XXX)",
            "[#xxx](https://github.com/nf-core/sarek/pull/xxx)",
        ):
            with self.subTest(reference=reference):
                self.git("restore", "--source", self.base, "--", "CHANGELOG.md")
                result = self.check_changelog(f"- {reference} - Add check\n")
                self.assertEqual(result.returncode, 1)
                self.assertIn("replace the placeholder", result.stderr)

    def test_missing_changelog_entry_fails(self):
        self.git("-c", "core.hooksPath=/dev/null", "commit", "--quiet", "--allow-empty", "-m", "No changelog")
        result = subprocess.run(
            [sys.executable, str(CHECK_CHANGELOG_PR), "--base", self.base, "--pr-number", "2468"],
            cwd=self.repo,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("add a CHANGELOG.md entry", result.stderr)

    def test_reference_outside_bullet_does_not_count(self):
        result = self.check_changelog("## [#2468](https://github.com/nf-core/sarek/pull/2468)\n")
        self.assertEqual(result.returncode, 1)

    def test_parent_reference_is_allowed_alongside_current_pr(self):
        result = self.check_changelog(
            "- [#1357](https://github.com/nf-core/sarek/pull/1357) - Parent change\n"
            "- [#2468](https://github.com/nf-core/sarek/pull/2468) - Add check\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_base_only_change_does_not_count(self):
        self.git("checkout", "--quiet", "-b", "base-advanced")
        self.changelog.write_text(self.changelog.read_text() + "- [#2468](https://github.com/nf-core/sarek/pull/2468) - Base change\n")
        self.commit_changelog()
        self.base = self.git("rev-parse", "HEAD").stdout.strip()
        self.git("checkout", "--quiet", "--detach", "HEAD~1")
        result = self.check_changelog("- [#1357](https://github.com/nf-core/sarek/pull/1357) - Head change\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn("add a CHANGELOG.md entry", result.stderr)


if __name__ == "__main__":
    unittest.main()
