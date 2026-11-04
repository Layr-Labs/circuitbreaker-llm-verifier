"""Select explicit ordinary or retained-result execution before any review job."""

import json
import os


MARKER = "yukon-retained-result-v1"


def _unique(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate evaluation key")
        result[key] = value
    return result


def select_route(raw):
    if raw in (None, ""):
        return "ordinary"
    if type(raw) is not str or len(raw.encode("utf-8")) > 65536:
        raise ValueError("invalid evaluation input")
    value = json.loads(raw, object_pairs_hook=_unique,
                       parse_constant=lambda _: (_ for _ in ()).throw(ValueError("nonfinite JSON")))
    if type(value) is not dict or not {"policies", "fields"} <= set(value) \
            or set(value) - {"policies", "fields", "admission", "retained_result"}:
        raise ValueError("invalid evaluation fields")
    fields = value["fields"]
    if type(fields) is not dict or set(fields) != {"scores", "validity"} \
            or any(type(flag) is not bool for flag in fields.values()) \
            or not any(fields.values()):
        raise ValueError("invalid requested fields")
    policies = value["policies"]
    if type(policies) is not dict or set(policies) != {"scoring", "validity"}:
        raise ValueError("invalid policy selection")
    scoring = policies["scoring"]
    marker = scoring.get("retained_result") if type(scoring) is dict else None
    if marker is not None and marker != MARKER:
        raise ValueError("unsupported retained-result capability")
    if fields["validity"]:
        if "retained_result" in value:
            raise ValueError("fresh validity cannot claim retained derivation")
        # Baselines and submissions without selected validity use this explicit
        # full-review branch. Descriptor absence never selects a paid fallback.
        return "ordinary"
    if marker != MARKER:
        raise ValueError("score-only execution requires the retained-result policy")
    retained = value.get("retained_result")
    if type(retained) is not dict or retained.get("version") != MARKER:
        raise ValueError("score-only execution is missing its retained-result descriptor")
    # The generic prepare phase validates the entire closed descriptor and its
    # artifact/plan bindings before constructing any provider.
    return "derive"


def main():
    from pathlib import Path
    from tools.derivation_profile import read_profile
    read_profile(Path(__file__).resolve().parent)
    route = select_route(os.environ.get("YUKON_EVALUATION", ""))
    output = os.environ.get("GITHUB_OUTPUT")
    if output:
        with open(output, "a", encoding="utf-8") as handle:
            handle.write("route=" + route + "\n")
    print(route)


if __name__ == "__main__":
    main()
