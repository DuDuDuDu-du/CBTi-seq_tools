import argparse
from cbtiseq_tools.pipeline import run_pipeline

def main():
    parser = argparse.ArgumentParser(
        prog="CBTi-seq_tools",
        description="CBTi-seq pipeline"
    )

    parser.add_argument("--workdir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)

    parser.add_argument("--bcA", default=None)
    parser.add_argument("--umiA", default=None)
    parser.add_argument("--bcB", default=None)
    parser.add_argument("--umiB", default=None)

    parser.add_argument("--whitelist", default=None)

    args = parser.parse_args()

    run_pipeline(args)

if __name__ == "__main__":
    main()
