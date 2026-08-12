import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("plan_weighted_shards.py")
SPEC = importlib.util.spec_from_file_location("plan_weighted_shards", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PlanWeightedShardsTest(unittest.TestCase):
    def test_pairs_discovery_rows_with_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            discovery = Path(directory) / "discovery.csv"
            log = Path(directory) / "dry-run.log"
            with discovery.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["filename", "test"])
                writer.writeheader()
                writer.writerow({"filename": "/workspace/tests/example.nf.test", "test": "first test"})
                writer.writerow({"filename": "/workspace/tests/example.nf.test", "test": "second test"})
            log.write_text("Test [1234abcd] 'first test' PASSED\nTest [9876fedc] 'second test' PASSED\n")

            tests = MODULE.read_discovery(discovery, log)

            self.assertEqual(tests[0]["selector"], "tests/example.nf.test@1234abcd")
            self.assertEqual(tests[1]["key"], "tests/example.nf.test::second test")

    def test_assigns_longest_tests_to_different_shards(self):
        tests = [
            {"key": "slow", "selector": "tests/a.nf.test@11111111"},
            {"key": "medium", "selector": "tests/b.nf.test@22222222"},
            {"key": "fast", "selector": "tests/c.nf.test@33333333"},
        ]

        shards = MODULE.assign_shards(tests, {"slow": 10, "medium": 6, "fast": 4}, 2)

        MODULE.validate_assignments(tests, shards)
        self.assertEqual(shards[0]["estimated_seconds"], 10)
        self.assertEqual(shards[1]["estimated_seconds"], 10)

    def test_rejects_mismatched_discovery_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            discovery = Path(directory) / "discovery.csv"
            log = Path(directory) / "dry-run.log"
            discovery.write_text("filename,test\n/workspace/tests/example.nf.test,example\n")
            log.write_text("")

            with self.assertRaisesRegex(ValueError, "1 tests.*0 hashes"):
                MODULE.read_discovery(discovery, log)


if __name__ == "__main__":
    unittest.main()
