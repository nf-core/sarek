import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("summarise_runtimes.py")
SPEC = importlib.util.spec_from_file_location("summarise_runtimes", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SummariseRuntimesTest(unittest.TestCase):
    def test_summarises_passing_observations(self):
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "results.csv"
            with report.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["filename", "testsuite", "type", "test", "result", "start", "end", "time"],
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "filename": "/workspace/tests/example.nf.test",
                            "testsuite": "Test pipeline",
                            "type": "PipelineTestSuite",
                            "test": "example",
                            "result": "PASSED",
                            "start": "0",
                            "end": "1000",
                            "time": "1.0",
                        },
                        {
                            "filename": "/workspace/tests/example.nf.test",
                            "testsuite": "Test pipeline",
                            "type": "PipelineTestSuite",
                            "test": "example",
                            "result": "PASSED",
                            "start": "0",
                            "end": "3000",
                            "time": "3.0",
                        },
                        {
                            "filename": "/workspace/tests/example.nf.test",
                            "testsuite": "Test pipeline",
                            "type": "PipelineTestSuite",
                            "test": "ignored failure",
                            "result": "FAILED",
                            "start": "0",
                            "end": "9000",
                            "time": "9.0",
                        },
                    ]
                )

            rows = MODULE.summarise(MODULE.read_reports(directory))

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["test"], "tests/example.nf.test::example")
            self.assertEqual(rows[0]["observations"], 2)
            self.assertEqual(rows[0]["median_seconds"], 2.0)
            self.assertEqual(rows[0]["p75_seconds"], 3.0)


if __name__ == "__main__":
    unittest.main()
