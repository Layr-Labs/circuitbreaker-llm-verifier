"""Representation-boundary checks; synthetic inputs assert no theorem validity."""

from copy import deepcopy
import unittest

from llm_verifier.errors import VerifierError

from expressions import CONVENTION, VERSION, VERSION_V2, polynomial, ratio
from score_derivation import validate_transcription


def original():
    return {
        "version": VERSION,
        "convention": CONVENTION,
        "family": "Synthetic fixed family; no mathematical acceptance.",
        "theorem": "Synthetic source text for a separately reviewed transcription.",
        "guarantee": polynomial(4, 1),
        "projection": {"kind": "finite", "coefficient": ratio(4), "exponent": ratio(1)},
    }


def term(coefficient, numerator, denominator=1, log=0, loglog=0):
    return {
        "coefficient": ratio(coefficient),
        "exponent": ratio(numerator, denominator),
        "log_exponent": ratio(log),
        "loglog_exponent": ratio(loglog),
    }


def target(bound, remainder):
    return {
        **deepcopy(bound),
        "version": VERSION_V2,
        "guarantee": {
            "kind": "eventual_lower_bound",
            "leading": deepcopy(bound["guarantee"]),
            "remainder": remainder,
        },
    }


class TranscriptionBoundaryTests(unittest.TestCase):
    def test_explicit_big_o_and_multiple_terms_preserve_original_bytes(self):
        source = original()
        frozen = deepcopy(source)
        remainders = [
            {"kind": "explicit", "terms": [term(182, 4, 5)]},
            {"kind": "explicit", "terms": [term(8, 1, 100), term(6, 3, 50000), term(23, 0)]},
            {"kind": "big_o", "exponent": ratio(0), "log_exponent": ratio(1),
             "loglog_exponent": ratio(0)},
            {"kind": "explicit", "terms": []},
        ]
        for remainder in remainders:
            with self.subTest(remainder=remainder):
                proposed = target(source, remainder)
                result = validate_transcription(source, proposed)
                self.assertEqual(result, proposed)
                self.assertIsNot(result, proposed)
                self.assertEqual(source, frozen)
                result["family"] = "mutated result"
                self.assertEqual(proposed["family"], frozen["family"])

    def test_metadata_and_projected_leading_changes_reject(self):
        source = original()
        proposed = target(source, {"kind": "explicit", "terms": [term(7, 1, 2)]})
        for name in ("family", "theorem"):
            with self.subTest(field=name):
                changed = deepcopy(proposed)
                changed[name] += " altered"
                with self.assertRaises(VerifierError):
                    validate_transcription(source, changed)
        changed = deepcopy(proposed)
        changed["guarantee"]["leading"] = polynomial(5, 1)
        changed["projection"]["coefficient"] = ratio(5)
        with self.assertRaises(VerifierError):
            validate_transcription(source, changed)
        changed = deepcopy(proposed)
        changed["guarantee"]["leading"] = polynomial(4, 2)
        changed["projection"]["exponent"] = ratio(2)
        with self.assertRaises(VerifierError):
            validate_transcription(source, changed)

    def test_same_projection_does_not_allow_a_different_log_leading_term(self):
        source = original()
        proposed = target(source, {"kind": "explicit", "terms": [term(7, 1, 2)]})
        proposed["guarantee"]["leading"]["log_exponent"] = ratio(1)
        with self.assertRaises(VerifierError):
            validate_transcription(source, proposed)

    def test_existing_eventual_guarantee_cannot_be_tightened(self):
        source = target(original(), {"kind": "explicit", "terms": [term(7, 1, 2)]})
        proposed = deepcopy(source)
        proposed["guarantee"]["remainder"]["terms"][0]["coefficient"] = ratio(6)
        with self.assertRaises(VerifierError):
            validate_transcription(source, proposed)
        self.assertEqual(validate_transcription(source, source), source)

    def test_unsupported_or_unchanged_declaration_retains_its_exact_version(self):
        source = original()
        self.assertEqual(validate_transcription(source, source), source)
        changed = {**source, "version": VERSION_V2}
        with self.assertRaises(VerifierError):
            validate_transcription(source, changed)

    def test_wrong_version_and_non_subleading_loss_reject(self):
        source = original()
        changed = target(source, {"kind": "explicit", "terms": [term(7, 1, 2)]})
        changed["version"] = VERSION
        with self.assertRaises(VerifierError):
            validate_transcription(source, changed)
        changed = target(source, {"kind": "explicit", "terms": [term(7, 1)]})
        with self.assertRaises(VerifierError):
            validate_transcription(source, changed)


class ScoreProjectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import json
        from pathlib import Path
        import tempfile
        from tests.test_refinements import synthetic_origin

        cls.temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        _, bundle, cls.public, cls.target = synthetic_origin(Path(cls.temporary.name))
        read = lambda path: json.loads(path.read_text())
        intake = read(bundle / "llmv-intake/intake.json")
        accepted = read(bundle / "verifier-accepted-results/accepted-results.json")
        cls.material = {"intake": intake, "accepted_results": accepted,
                        "manifest": intake["data"]["manifest"]}
        cls.original_score = read(bundle / "circuitbreaker-score/.yukon/score.json")

    def result(self, status="derived"):
        from score_derivation import SCORING_POLICY
        return {
            "sha256": "d" * 64,
            "simulation": True,
            "data": {
                "intake": {"sha256": "e" * 64, "data": {"evaluation": {
                    "policies": {"scoring": deepcopy(SCORING_POLICY)},
                    "retained_result": {"version": "yukon-retained-result-v1",
                                        "validity": {"valid": True, "result_id": "synthetic-old-validity"}},
                }}},
                "decision": {"status": status, "reason": "Synthetic projection check only.",
                             "payload": {"target_bound": self.target} if status == "derived" else None,
                             "source_refs": [{"source": "claim:0"}] if status == "derived" else []},
            },
        }

    def test_retains_original_rich_metric_and_binds_new_frontier_to_scoring(self):
        from score_derivation import _project_verified
        before = deepcopy(self.material)
        result = self.result()
        output = _project_verified(result, self.material)
        self.assertEqual(output["metrics"]["circuitbreaker"], self.original_score["metrics"]["circuitbreaker"])
        self.assertEqual(output["scores"], self.original_score["scores"])
        self.assertNotIn("valid", output)
        self.assertEqual(output["derived_from"], result["data"]["intake"]["data"]["evaluation"]["retained_result"])
        frontier = output["metrics"]["circuitbreaker_frontier"]
        self.assertEqual(frontier["version"], "circuitbreaker-frontier-v2")
        self.assertEqual(frontier["resources"]["circuit_size"]["bound"], self.target["guarantee"])
        self.assertEqual(frontier["resources"]["circuit_size"]["statement"], self.material["manifest"]["claims"][0]["statement"])
        self.assertEqual(frontier["dossier_sha256"], result["sha256"])
        self.assertEqual(frontier["intake_sha256"], result["data"]["intake"]["sha256"])
        self.assertEqual(self.material, before)

    def test_unsupported_retains_original_v1_without_guessing_loss(self):
        from score_derivation import _project_verified
        output = _project_verified(self.result("unsupported"), self.material)
        frontier = output["metrics"]["circuitbreaker_frontier"]
        self.assertEqual(frontier["version"], "circuitbreaker-frontier-v1")
        self.assertEqual(frontier["resources"]["circuit_size"]["bound"],
                         self.original_score["metrics"]["circuitbreaker"]["bound"]["guarantee"])
        self.assertEqual(output["metrics"]["circuitbreaker_derivation"]["status"], "unsupported")

    def test_main_claim_cannot_be_silently_dropped_as_not_applicable(self):
        from score_derivation import _project_verified
        with self.assertRaises(VerifierError):
            _project_verified(self.result("not_applicable"), self.material)

    def test_rejected_source_has_null_scores_and_no_competitive_metric(self):
        from score_derivation import _project_verified
        result = self.result("not_applicable")
        result["data"]["intake"]["data"]["evaluation"]["retained_result"]["validity"]["valid"] = False
        output = _project_verified(result, None)
        self.assertIsNone(output["scores"])
        self.assertEqual(set(output["metrics"]), {"circuitbreaker_derivation"})
        self.assertNotIn("valid", output)
        result["data"]["decision"]["status"] = "derived"
        result["data"]["decision"]["payload"] = {"target_bound": self.target}
        with self.assertRaises(VerifierError):
            _project_verified(result, None)

    def test_unscored_lemma_stays_unscored(self):
        from score_derivation import _project_verified
        material = deepcopy(self.material)
        material["manifest"]["claims"] = ["Synthetic supporting lemma, not a competitive bound."]
        output = _project_verified(self.result("not_applicable"), material)
        self.assertIsNone(output["scores"])
        self.assertEqual(set(output["metrics"]), {"circuitbreaker_derivation"})
        with self.assertRaises(VerifierError):
            _project_verified(self.result("derived"), material)

    def test_legacy_review_cannot_admit_an_original_v2_declaration(self):
        import json
        from expressions import claims_for
        from score_derivation import _project_verified
        material = deepcopy(self.material)
        bound = deepcopy(self.target)
        material["intake"]["data"]["files"]["bound.json"]["text"] = json.dumps(bound)
        material["manifest"]["claims"] = claims_for(bound)
        with self.assertRaisesRegex(VerifierError, "legacy_v2"):
            _project_verified(self.result("unsupported"), material)

    def test_file_only_citation_cannot_enlarge_accepted_joint_claims(self):
        from score_derivation import _project_verified
        result = self.result()
        result["data"]["decision"]["source_refs"] = [{"source": "file:proof.md"}]
        with self.assertRaises(VerifierError):
            _project_verified(result, self.material)


class WorkflowRoutingTests(unittest.TestCase):
    def evaluation(self, scores=True, validity=False):
        from score_derivation import SCORING_POLICY
        return {"policies": {"scoring": deepcopy(SCORING_POLICY),
                             "validity": {"version": "fixed-family-full-b2-v1"}},
                "fields": {"scores": scores, "validity": validity}}

    def route(self, value):
        import json
        from select_scoring_route import select_route
        return select_route(json.dumps(value))

    def test_baseline_and_missing_validity_explicitly_request_ordinary_work(self):
        self.assertEqual(self.route(self.evaluation(validity=True)), "ordinary")
        self.assertEqual(self.route(self.evaluation(scores=False, validity=True)), "ordinary")

    def test_score_only_requires_descriptor_and_never_falls_back_to_full_review(self):
        value = self.evaluation()
        for descriptor in (None, {}, {"version": "unknown"}):
            with self.subTest(descriptor=descriptor):
                value["retained_result"] = descriptor
                with self.assertRaises(ValueError):
                    self.route(value)
        value.pop("retained_result")
        with self.assertRaises(ValueError):
            self.route(value)
        value["retained_result"] = {"version": "yukon-retained-result-v1"}
        self.assertEqual(self.route(value), "derive")

    def test_wrong_capability_fresh_validity_and_bad_flags_reject(self):
        value = self.evaluation()
        value["retained_result"] = {"version": "yukon-retained-result-v1"}
        value["policies"]["scoring"]["retained_result"] = "unknown"
        with self.assertRaises(ValueError):
            self.route(value)
        value = self.evaluation(validity=True)
        value["retained_result"] = {"version": "yukon-retained-result-v1"}
        with self.assertRaises(ValueError):
            self.route(value)
        for flag in (0, 1, "false", None):
            value = self.evaluation()
            value["fields"]["validity"] = flag
            with self.subTest(flag=flag), self.assertRaises(ValueError):
                self.route(value)

    def test_duplicate_keys_are_not_allowed_to_change_route(self):
        from select_scoring_route import select_route
        with self.assertRaises(ValueError):
            select_route('{"fields":{},"fields":{"validity":true,"scores":true},"policies":{}}')
