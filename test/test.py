#! /usr/bin/python3

import shlex
import subprocess
import sys
import os
import resource
import time
import tempfile
import argparse
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Union, List, Dict, Tuple, TextIO
from dataclasses import dataclass


BIN_DIR = "bin/test"

_MAXRSS_SCALE = 1 if sys.platform == "darwin" else 1024

@dataclass
class TestResult:
    ret: int
    stdout: str
    stderr: str
    runtime: float    # wall seconds
    cpu_time: float   # user + system seconds
    max_rss: int      # peak resident set, bytes

Results = Dict[Tuple[str, Union[int, None]], TestResult]


@dataclass
class TestCase:
    name: str
    input: Union[List[int], None]

tests = [
    TestCase("Taylor-Green-MG",  [8, 16, 64]),
    TestCase("Taylor-Green-FFT", [8, 16, 64]),
    TestCase("Multigrid",        [32, 64, 128, 512, 1024]),
    TestCase("Polar-Couette",    [8, 16, 32]),
    TestCase("Channel-MG",       [16, 32, 64]),
    TestCase("Channel-FFT",      [16, 32, 64]),
]


def run_cmd(cmd, echo=False):
    if echo:
        print(f"[CMD] {' '.join(map(shlex.quote, cmd))}")
    ret = subprocess.run(cmd, capture_output=True)
    return ret.returncode, ret.stdout.decode("utf-8"), ret.stderr.decode("utf-8")


def build_tests(jobs: int = 1, parallel: bool = False) -> bool:
    targets = [f"{BIN_DIR}/{t.name}" for t in tests]

    print(f"[INFO] Build {len(targets)} target(s) with -j{jobs}")
    parallel_flag = "PARALLEL=1" if parallel else "PARALLEL=0"
    ret, stdout, stderr = run_cmd(["make", parallel_flag, f"-j{jobs}", *targets])
    if ret != 0:
        print("Build failed", file=sys.stderr)
        print(stderr, file=sys.stderr)
        return False
    return True


def run_name(name: str, inp: Union[int, None]):
    return name if inp is None else f"{name}-{inp}"


_print_lock = threading.Lock()
def log(msg: str) -> None:
    with _print_lock:
        sys.stdout.write(msg + "\n")
        sys.stdout.flush()


def run_test(cmd: List[str], name: str) -> TestResult:
    log(f"[INFO] Run {name}...")

    with tempfile.TemporaryFile() as out, tempfile.TemporaryFile() as err:
        start = time.perf_counter()
        proc = subprocess.Popen(cmd, stdout=out, stderr=err, env={"OMP_NUM_THREADS": "4"})
        _pid, status, ru = os.wait4(proc.pid, 0)
        runtime = time.perf_counter() - start

        proc.returncode = os.waitstatus_to_exitcode(status)

        out.seek(0)
        err.seek(0)
        stdout = out.read().decode("utf-8")
        stderr = err.read().decode("utf-8")

    return TestResult(
        ret=proc.returncode,
        stdout=stdout,
        stderr=stderr,
        runtime=runtime,
        cpu_time=ru.ru_utime + ru.ru_stime,
        max_rss=ru.ru_maxrss * _MAXRSS_SCALE,
    )


def run_tests(jobs: int = 1) -> Results:
    work: List[Tuple[Tuple[str, Union[int, None]], List[str]]] = []
    for test in tests:
        exe = f"{BIN_DIR}/{test.name}"
        if test.input is None:
            work.append(((test.name, None), [exe]))
        else:
            for n in test.input:
                work.append(((test.name, n), [exe, f"{n}"]))

    if jobs == 1:
        return {key: run_test(cmd, run_name(*key)) for key, cmd in work}

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(run_test, cmd, run_name(*key)) for key, cmd in work]
        return {key: fut.result() for (key, _), fut in zip(work, futures)}


_UNITS = ("B", "kB", "MB", "GB", "TB", "PB", "EB")
def format_bytes(n: int, precision: int = 1, base: int = 1024) -> str:
    if n < 0:
        return "-" + format_bytes(-n, precision, base)

    size = float(n)
    unit = _UNITS[0]
    for i, unit in enumerate(_UNITS):
        if i == len(_UNITS) - 1 or round(size, precision) < base:
            break
        size /= base

    if unit == "B":
        return f"{int(size)}{unit}"

    text = f"{size:.{precision}f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return f"{text}{unit}"


GREEN = "\033[32m"
RED   = "\033[31m"
RESET = "\033[0m"
STATUS_WIDTH = 4
MIN_DOTS = 6         # minimum dots next to the longest name
METRICS_INDENT = 6   # spaces before the metrics line, after the leading space

def render_results(results: Results, file: TextIO = sys.stdout) -> None:
    color = file.isatty()

    rows: List[Tuple[str, bool, str]] = []
    for (name, inp), res in results.items():
        metrics = (
            f"wall={res.runtime:.2f}s, "
            f"cpu={res.cpu_time:.2f}s, "
            f"mem={format_bytes(res.max_rss)}"
        )
        rows.append((run_name(name, inp), res.ret == 0, metrics))

    if not rows:
        return

    longest_name = max(len(n) for n, _, _ in rows)
    longest_metrics = max(len(m) for _, _, m in rows)

    inner = max(
        longest_name + 1 + MIN_DOTS + 1 + STATUS_WIDTH,
        METRICS_INDENT + longest_metrics + STATUS_WIDTH,
    )

    border = "=" * inner
    print(f" {border}", file=file)

    for name, ok, metrics in rows:
        dots = "." * (inner - len(name) - 2 - STATUS_WIDTH)
        status = "PASS" if ok else "FAIL"
        if color:
            status = f"{GREEN if ok else RED}{status}{RESET}"
        print(f" {name} {dots} {status}", file=file)
        print(f" {'':{METRICS_INDENT}}{metrics}", file=file)

    print(f" {border}", file=file)


def dump_failed(results: Results):
    for (name, inp), res in results.items():
        if res.ret == 0:
            continue
        print("-"*100)
        print(f"{run_name(name, inp)} stderr:")
        print(res.stderr)
        print("-"*100)


def parse_args(argv: Union[List[str], None] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build and run the test suite.")
    p.add_argument(
        "-j", "--jobs",
        type=int, nargs="?", const=os.cpu_count() or 1, default=1, metavar="N",
        help="run N build/test jobs in parallel "
             "(default: 1; bare -j uses all %d cores)" % (os.cpu_count() or 1),
    )
    p.add_argument(
        "-p", "--parallel",
        action="store_true",
        help="Build the test case with PARALLEL=1",
    )
    args = p.parse_args(argv)
    if args.jobs < 1:
        p.error("-j must be at least 1")
    return args


def main():
    args = parse_args()
    if not build_tests(args.jobs, args.parallel):
        sys.exit(1)
    results = run_tests(args.jobs)
    render_results(results)
    dump_failed(results)


if __name__ == "__main__":
    main()
