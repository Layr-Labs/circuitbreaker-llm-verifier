"""V2 producer integration through pinned, HTTP-stubbed verifier phase chains."""

import copy
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from export_metrics import (LEGACY_POLICIES, LEGACY_REVIEW, projection,
                            validate_contract)
from expressions import (CONVENTION, POLICIES, REVIEW, VERSION, VERSION_V2,
                         claims_for, polynomial, ratio, validate_bound)
from llm_verifier.errors import VerifierError
from tests.remainder_cases import big_o, explicit, term
from tests.synthetic import synthetic_dossier

ROOT = Path(__file__).resolve().parents[1]


class RemainderExportTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        (self.root / 'submission').mkdir()
        for name in ('verifier.json', 'rules.md'):
            shutil.copyfile(ROOT / name, self.root / name)
        shutil.copytree(ROOT / 'checkers', self.root / 'checkers')
        self.bound = {
            'version': VERSION_V2, 'convention': CONVENTION,
            'family': 'SYNTHETIC integration fixture; no circuit family asserted.',
            'theorem': 'SYNTHETIC transport test; no mathematical acceptance.',
            'guarantee': explicit(term(7)),
            'projection': {'kind': 'finite', 'coefficient': ratio(4), 'exponent': ratio(1)},
        }

    def dossier(self, *, manifest=None, **options):
        if manifest is None:
            manifest = {'version': 1, 'title': 'Synthetic v2 producer test',
                        'claims': claims_for(self.bound), 'experiments': [], 'checks': []}
        return synthetic_dossier(self.root, self.bound, manifest, **options)

    def legacy_config(self):
        for name in ('verifier.json', 'rules.md'):
            shutil.copyfile(ROOT / 'tests/fixtures/legacy-v1' / name, self.root / name)

    def test_v2_remainders_and_legacy_kinds_export_exact_versioned_tuple(self):
        for guarantee in (explicit(), explicit(term(7)), big_o(0, 1), polynomial(4, 1)):
            with self.subTest(guarantee=guarantee):
                self.bound['guarantee'] = guarantee
                record = self.dossier()
                output = projection(record)
                rich = output['metrics']['circuitbreaker']
                table = output['metrics']['circuitbreaker_frontier']
                resource = table['resources']['circuit_size']
                self.assertEqual(output['scores'], [4, 1])
                self.assertEqual(rich['version'], VERSION_V2)
                self.assertEqual(rich['bound'], self.bound)
                self.assertEqual((table['version'], table['comparator']),
                                 ('circuitbreaker-frontier-v2', 'circuitbreaker-v2'))
                self.assertEqual((table['scope'], table['model'], table['notation']),
                                 (CONVENTION, 'full-b2-v1', 'lower-bound'))
                self.assertEqual(resource['bound'], guarantee)
                self.assertEqual(resource['statement'], rich['statements'][0])
                self.assertEqual(resource['statement'], rich['statements'][1])
                self.assertEqual(resource['source_sha256'], rich['score_records'][0]['source_sha256'])
                self.assertEqual(resource['evidence'], rich['evidence'])
                self.assertEqual(table['dossier_sha256'], record['sha256'])
                self.assertEqual(table['intake_sha256'], rich['intake_sha256'])
                self.assertEqual(rich['selected_review_sha256'], record['data']['selected_review_sha256'])

    def test_legacy_declaration_under_current_review_stays_v1(self):
        self.bound.update(version=VERSION, guarantee=polynomial(4, 1))
        result = projection(self.dossier())
        self.assertEqual(result['metrics']['circuitbreaker']['version'], VERSION)
        table = result['metrics']['circuitbreaker_frontier']
        self.assertEqual((table['version'], table['comparator']),
                         ('circuitbreaker-frontier-v1', 'circuitbreaker-v1'))

    def test_legacy_review_rejects_all_v2_declarations(self):
        self.legacy_config()
        for guarantee in (explicit(term(7)), polynomial(4, 1)):
            with self.subTest(guarantee=guarantee):
                self.bound['guarantee'] = guarantee
                with self.assertRaisesRegex(VerifierError, 'legacy review cannot accept a v2'):
                    projection(self.dossier())

    def test_no_v1_eventual_or_unknown_version_can_export(self):
        manifest = {'version': 1, 'title': 'Synthetic stale claims',
                    'claims': claims_for(self.bound), 'experiments': [], 'checks': []}
        for version in (VERSION, 'circuitbreaker-expression-v3'):
            with self.subTest(version=version):
                self.bound['version'] = version
                with self.assertRaisesRegex(VerifierError, 'circuitbreaker_bound_contract'):
                    projection(self.dossier(manifest=manifest))

    def test_same_scores_do_not_authorize_different_remainder_or_leading_only_claim(self):
        manifest = {'version': 1, 'title': 'Synthetic stale claims',
                    'claims': claims_for(self.bound), 'experiments': [], 'checks': []}
        for guarantee in (explicit(term(6)), big_o(), polynomial(4, 1)):
            with self.subTest(guarantee=guarantee):
                self.bound['guarantee'] = guarantee
                with self.assertRaisesRegex(VerifierError, 'joint statement/score/evidence mismatch'):
                    projection(self.dossier(manifest=manifest))
        self.bound['guarantee'] = polynomial(4, 1)
        manifest['claims'] = claims_for(self.bound)
        self.bound['guarantee'] = explicit(term(7))
        with self.assertRaisesRegex(VerifierError, 'joint statement/score/evidence mismatch'):
            projection(self.dossier(manifest=manifest))

    def test_v2_endpoint_projections_preserve_null_scores_and_full_loss(self):
        for guarantee, projected in (
            (explicit(term(10, 1), leading=polynomial(1, 1, 1)),
             {'kind': 'lane_endpoint', 'exponent': ratio(1)}),
            ({'kind': 'superpolynomial'}, {'kind': 'global_endpoint'}),
        ):
            with self.subTest(guarantee=guarantee):
                self.bound.update(guarantee=guarantee, projection=projected)
                output = projection(self.dossier())
                self.assertIsNone(output['scores'])
                self.assertEqual(output['metrics']['circuitbreaker_frontier']['resources']['circuit_size']['bound'], guarantee)
                self.assertEqual(output['metrics']['circuitbreaker']['score_records'][0]['record'], 'positive_infinity')

    def test_rejection_and_simulation_never_export_v2_metrics(self):
        self.assertNotIn('metrics', projection(self.dossier(valid=False)))
        with self.assertRaisesRegex(VerifierError, 'simulation'):
            projection(self.dossier(simulation=True))

    def test_only_matched_review_policy_pairs_are_accepted(self):
        current = json.loads((ROOT / 'verifier.json').read_text())
        self.assertEqual(current['review_standard'], REVIEW)
        self.assertEqual(current['policies'], POLICIES)
        self.assertFalse(validate_contract(current))
        old = copy.deepcopy(current)
        old.update(review_standard=LEGACY_REVIEW, policies=LEGACY_POLICIES)
        self.assertTrue(validate_contract(old))
        for review, policies in ((REVIEW, LEGACY_POLICIES), (LEGACY_REVIEW, POLICIES)):
            bad = copy.deepcopy(current)
            bad.update(review_standard=review, policies=policies)
            with self.assertRaisesRegex(VerifierError, 'objective_contract'):
                validate_contract(bad)
        self.assertEqual(current['supported_policies']['scoring'][0]['selector'], POLICIES['scoring'])
        self.assertEqual(current['supported_policies']['validity'][0]['selector'], POLICIES['validity'])
        benchmark = json.loads((ROOT / 'benchmark.json').read_text())['tracks'][0]
        self.assertEqual(benchmark['paretoPolicies'], POLICIES)

    def test_all_documented_declarations_validate_and_bind_complete_claims(self):
        examples = json.loads((ROOT / 'examples/remainders-v2.json').read_text())
        self.assertIs(examples['synthetic'], True)
        for example in examples['examples']:
            with self.subTest(example=example['id']):
                bound = example['bound']
                validate_bound(bound)
                claims = claims_for(bound)
                self.assertEqual(claims[0]['statement'], claims[1]['statement'])
                self.assertEqual(claims[0]['evidence'], ['bound.json', 'proof.md'])
                self.assertIn('guarantee.remainder', claims[0]['statement'])

    def test_v2_schema_shape_and_closed_remainder_contract(self):
        bound = json.loads((ROOT / 'schemas/bound-v2.schema.json').read_text())
        rich = json.loads((ROOT / 'schemas/accepted-metric-v2.schema.json').read_text())
        table = json.loads((ROOT / 'schemas/frontier-v2.schema.json').read_text())
        self.assertEqual(bound['properties']['version']['const'], VERSION_V2)
        self.assertEqual(rich['properties']['version']['const'], VERSION_V2)
        self.assertEqual(rich['properties']['bound']['$ref'], 'bound-v2.schema.json')
        self.assertEqual(table['properties']['version']['const'], 'circuitbreaker-frontier-v2')
        self.assertEqual(table['properties']['comparator']['const'], 'circuitbreaker-v2')
        self.assertEqual(table['properties']['resources']['properties']['circuit_size']['properties']['bound']['$ref'],
                         'bound-v2.schema.json#/properties/guarantee')
        for schema in (bound, rich, table):
            self.assertFalse(schema['additionalProperties'])
            self.assertEqual(set(schema['required']), set(schema['properties']))
        eventual = bound['properties']['guarantee']['oneOf'][-1]
        self.assertEqual(eventual['properties']['kind']['const'], 'eventual_lower_bound')
        self.assertFalse(eventual['additionalProperties'])
        self.assertEqual(eventual['required'], ['kind', 'leading', 'remainder'])
        loss, unknown = eventual['properties']['remainder']['oneOf']
        self.assertFalse(loss['additionalProperties'])
        self.assertFalse(unknown['additionalProperties'])
        self.assertEqual(set(unknown['required']), {'kind', 'exponent', 'log_exponent', 'loglog_exponent'})
        self.assertNotIn('coefficient', unknown['properties'])
        self.assertEqual((loss['properties']['terms']['minItems'], loss['properties']['terms']['maxItems']), (0, 8))
        self.assertNotIn('kind', bound['$defs']['monomial']['properties'])
        self.assertFalse(bound['$defs']['monomial']['additionalProperties'])

    def test_historic_v1_schemas_fixtures_and_baseline_bytes_stay_frozen(self):
        # Original base 2cd3a17; portable to CI source archives without .git.
        expected = {
            "schemas/accepted-metric-v1.schema.json": "cec36afaa6b2f66709af71f1fac5c3d0bf3c95c9d3110b458c26c631d9883de3",
            "schemas/bound-v1.schema.json": "fec3907d91f19c1bfb6a29613f280bac5a2e0e0ecdfa70ee7fabf172179d0b33",
            "schemas/frontier-v1.schema.json": "d2a17ae11e3d9fc1df80f603d20881b2b33b2fe5ad91f7f319d464cf3d6c2f40",
            "submission/bound.json": "0b19620284317cbeb66782da4251fd852d2c0c810cf4a9e4763c31b4e56c4556",
            "submission/proof.md": "bcc793e84a7a44f5bbe1ed7c0d9057ea909346282793978d58743b63708e2f7e",
            "submission/submission.json": "92b3518ff6a94c984c6f4e52429e4537ea5934fc9f52cae7036caa11707ec647",
            "tests/fixtures/candidate-3_5-provenance.json": "6a80ebe02a67f1becd71bbc95da5bac2bdf9a19ee919a6bde8ac0ce0bb9e285f",
            "tests/fixtures/candidate-3_5/audit.py": "cb2d6f698eaf9547b27051115219e4b1626d47c30a92570616a300d65acf4bef",
            "tests/fixtures/candidate-3_5/bound.json": "1a9ea745225a27a228aff4e4183fc9824e7724ac97997a228a6ea84200d5e48c",
            "tests/fixtures/candidate-3_5/disperser.md": "5a82a17623db98cdf44730eaafa4bc5d02096d340702ff86d7ba16a39ae9ef2d",
            "tests/fixtures/candidate-3_5/elimination.md": "6b14cf339736d34aa934c8135444233b3da5e1c681e26cabec4c391c681ec24c",
            "tests/fixtures/candidate-3_5/family.py": "b51fc1ce74a62147e1e253f2b682d713d081197a7c18ae8c72332368fe2bceb4",
            "tests/fixtures/candidate-3_5/proof.md": "9b188eebcaa0d8dde054588eaa91fac10e9e526231f7df08fdf28fa67fb4806d",
            "tests/fixtures/candidate-3_5/submission.json": "59b8b179abffb36ab43af429145e8d2b369caf7360d6700fc94d0f62be17d6eb",
            "tests/fixtures/circuitbreaker-frontier-v1.json": "29896a4de110a34ac7c46d9f0db55406d9172b10f054e12d4c0ddab1d97be69f"
        }
        for name, digest in expected.items():
            with self.subTest(path=name):
                self.assertEqual(hashlib.sha256((ROOT / name).read_bytes()).hexdigest(), digest)
