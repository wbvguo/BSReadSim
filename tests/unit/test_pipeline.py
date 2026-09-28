"""Thread planning for a run, and its batches processed in worker processes, in order."""

from dataclasses import dataclass
import os
import time
import unittest

from bsreadsim.errors import BSReadSimError
from bsreadsim.pipeline import ThreadPlan, _formatted_batches, thread_plan


class ThreadPlanTests(unittest.TestCase):
    def test_threads_are_split_without_exceeding_the_budget(self):
        for threads in range(1, 257):
            for output_format in ("fastq", "fastq.gz", "bam", "sam"):
                plan = thread_plan(threads, output_format)
                self.assertGreaterEqual(min(plan.workers, plan.core), 1)
                self.assertLessEqual(plan.bam, 64 if output_format == "bam" else 0)
                if threads > 2:  # compression threads count a third: they idle most of the time
                    self.assertLessEqual(plan.workers + plan.core + plan.bam // 3, threads)
        self.assertEqual(thread_plan(1, "fastq.gz"), ThreadPlan(1, 1, 0))
        self.assertEqual(thread_plan(1, "bam"), ThreadPlan(1, 1, 0))

    def test_measured_plans(self):
        # Fastest plans measured on the example genome (2M WGBS reads)
        self.assertEqual(thread_plan(4, "bam"), ThreadPlan(workers=3, core=1, bam=1))
        self.assertEqual(thread_plan(8, "bam"), ThreadPlan(workers=5, core=2, bam=3))
        self.assertEqual(thread_plan(4, "fastq.gz"), ThreadPlan(workers=3, core=1, bam=0))
        self.assertEqual(thread_plan(4, "fastq"), ThreadPlan(workers=2, core=2, bam=0))


@dataclass(frozen=True)
class Repeat:
    times: int

    def __call__(self, payload):
        time.sleep((8 - payload[0] % 8) * 0.002)  # later payloads tend to finish first
        return payload * self.times


def fail_on_three(payload):
    if payload[0] == 3:
        raise ValueError("no threes")
    return payload


def die_on_three(payload):
    if payload[0] == 3:
        os._exit(7)
    return payload


class WorkerProcessTests(unittest.TestCase):
    def test_results_keep_the_batch_order(self):
        payloads = [bytes([index]) * (1 + index * 50_000) for index in range(40)]
        self.assertEqual(list(_formatted_batches(payloads, Repeat(2), 3)),
                         [payload * 2 for payload in payloads])

    def test_an_exception_in_a_worker_is_raised(self):
        with self.assertRaisesRegex(ValueError, "no threes") as caught:
            list(_formatted_batches((bytes([index]) for index in range(10)), fail_on_three, 2))
        self.assertIn("fail_on_three", str(caught.exception.__cause__))

    def test_a_worker_that_dies_is_reported(self):
        with self.assertRaisesRegex(BSReadSimError, "worker process stopped"):
            list(_formatted_batches((bytes([index]) for index in range(10)), die_on_three, 2))

    def test_an_error_reading_the_stream_is_raised(self):
        def payloads():
            yield b"\x00"
            raise BSReadSimError("htsim stream was truncated")

        with self.assertRaisesRegex(BSReadSimError, "truncated"):
            list(_formatted_batches(payloads(), Repeat(1), 2))

    def test_leaving_early_stops_the_workers(self):
        batches = _formatted_batches(
            (bytes([index % 256]) * 100_000 for index in range(1000)), Repeat(1), 2)
        next(batches)
        batches.close()  # returns once the reader and the workers have stopped


if __name__ == "__main__":
    unittest.main()
