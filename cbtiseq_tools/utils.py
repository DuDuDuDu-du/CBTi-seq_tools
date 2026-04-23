import subprocess
import os


def run_cmd(cmd, step_name):
    print(f"[RUNNING] {step_name}")
    print(cmd)

    result = subprocess.run(
        cmd,
        shell=True,
        text=True,
        capture_output=True
    )

    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"{step_name} failed")

    print(f"[DONE] {step_name}\n")


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
