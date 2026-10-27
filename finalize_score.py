"""Thin local benchmark adapter: shared finalization followed by metric export."""

import argparse

from llm_verifier import pipeline
from llm_verifier.files import clear_output, distinct_paths, read_json

from export_metrics import export


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("intake", "initial", "initial-review", "score", "dossier"):
        parser.add_argument("--" + field, required=True)
    for field in ("followup", "final-review"):
        parser.add_argument("--" + field)
    args = parser.parse_args()
    inputs = list(filter(None, (args.intake, args.initial, args.initial_review,
                               args.followup, args.final_review)))
    distinct_paths(inputs, [args.score, args.dossier])
    clear_output(args.score)
    clear_output(args.dossier)
    try:
        pipeline.finalize(
            read_json(args.intake), read_json(args.initial), read_json(args.initial_review),
            args.score, args.dossier,
            read_json(args.followup) if args.followup else None,
            read_json(args.final_review) if args.final_review else None,
        )
        export(args.dossier, args.score, intake_path=args.intake)
    except Exception:
        clear_output(args.score)
        raise


if __name__ == "__main__":
    main()
