"""Export a same-run verified derivation; original validity is never emitted."""

import argparse

from llm_verifier.files import clear_output, distinct_paths, read_json, write_json

from score_derivation import projection


def export(result_path, output, *, intake_path=None, preview=False):
    inputs = [result_path] + ([intake_path] if intake_path else [])
    distinct_paths(inputs, [output])
    clear_output(output)
    expected = read_json(intake_path)["sha256"] if intake_path else None
    result = projection(read_json(result_path), expected_intake=expected, preview=preview)
    write_json(output, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--derivation", required=True)
    parser.add_argument("--intake")
    parser.add_argument("--out", required=True)
    parser.add_argument("--preview", action="store_true",
                        help="nonpublishing preview without a top-level scores field")
    args = parser.parse_args()
    export(args.derivation, args.out, intake_path=args.intake, preview=args.preview)


if __name__ == "__main__":
    main()
