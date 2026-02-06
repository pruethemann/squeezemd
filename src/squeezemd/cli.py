"""Console entrypoint placeholder for squeezemd.

This module currently provides a minimal CLI stub. It is kept for
packaging and to reserve the top-level entrypoint while the full
command interface is implemented elsewhere (e.g., Snakemake wrapper).
"""
import argparse
import sys


def main():
    """Parse arguments and print a placeholder message.

    Notes
    -----
    This function intentionally does not implement behavior yet. It
    exists to keep the CLI entrypoint stable during development.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('_', nargs='*')
    args = parser.parse_args()

    print("Arguments: " + str(args._))
    print("Replace this message by putting your code into "
          "squeezemd.cli.main")
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
