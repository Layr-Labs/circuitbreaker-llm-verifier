"""Exact baseline/candidate transport regressions, never scientific judgments.

The seven-file fixture is byte-for-byte original candidate evidence. Native HTTP
and Docker boundaries are mocked; recorded finite computation is not a proof.
"""

import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
from threading import Barrier, Event
import unittest
from unittest.mock import patch

from llm_verifier import pipeline
from llm_verifier.artifacts import intake_snapshot
from llm_verifier.config import load_config
from llm_verifier.jsonutil import canonical
from llm_verifier.runner import DockerRunner, output_record

from check_submission import check_submission
from export_metrics import projection
from tests.synthetic import audit_review, request_parts, synthetic_review
from tools.generate_hosted import VERIFIER_REVISION

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'tests/fixtures/candidate-3_5'
PROVENANCE = json.loads((ROOT / 'tests/fixtures/candidate-3_5-provenance.json').read_text())


class EscalationTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        for name in ('verifier.json', 'rules.md', 'checker-tests.json'):
            shutil.copyfile(ROOT / name, self.root / name)
        shutil.copytree(ROOT / 'checkers', self.root / 'checkers')

    def prepare(self, candidate=False):
        source = FIXTURE if candidate else ROOT / 'submission'
        shutil.copytree(source, self.root / 'submission')
        before = {p.name: p.read_bytes() for p in source.iterdir()}
        if candidate:
            self.assertEqual(set(before), set(PROVENANCE['files']))
            for name, data in before.items():
                self.assertEqual({'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()},
                                 PROVENANCE['files'][name])
        check_submission(self.root)
        effective = load_config(self.root, 'verifier.json', 'submission')
        intake = intake_snapshot(effective, self.root / 'submission')
        self.assertEqual(len(intake['data']['files']), 7 if candidate else 3)
        self.assertEqual(before, {name: r['text'].encode() for name, r in intake['data']['files'].items()})
        fixture = json.loads((ROOT / 'checker-tests.json').read_text())[0]
        result = PROVENANCE['initial_result'] if candidate else {**fixture['expected'], 'message': 'Synthetic fixture'}
        with patch.object(DockerRunner, 'run', return_value=output_record('completed', 0, canonical(result), b'')) as docker:
            initial = pipeline.execute(intake)
            self.assertEqual(docker.call_count, 1)
        return intake, initial

    def followup(self, intake, initial, first):
        self.assertEqual(len(first['data']['decision']['requests']), 1)
        fixture = json.loads((ROOT / 'checker-tests.json').read_text())[0]
        output = output_record('completed', 0, canonical({**fixture['expected'], 'message': 'Synthetic fixture'}), b'')
        with patch.object(DockerRunner, 'run', return_value=output) as docker:
            result = pipeline.execute(intake, first, initial)
            self.assertEqual(docker.call_count, 1)
        return result

    def measure(self, name, record, requests):
        intake = record['data']['phases'][0]
        phases = [p for p in record['data']['phases'] if p['stage'] in ('initial_review', 'final_review')]
        sizes = [{'role': request_parts(w)[0]['id'],
                  'region': request_parts(w)[1].get('region', {}).get('id'),
                  'prompt_bytes': len(canonical(w['input']))} for w in requests]
        self.assertLessEqual(max(s['prompt_bytes'] for s in sizes), 2097152)
        self.assertLess(len(canonical(record)), 16 * 1024 * 1024)
        for wire in requests:
            self.assertEqual((wire['model'], wire['reasoning'], wire['background'], wire['max_output_tokens']),
                             ('gpt-6-astra', {'effort': 'high'}, True, 32768))
        executions = []
        for phase in phases:
            self.assertEqual(phase['data']['execution']['version'], 'parallel-children-v1')
            for child in ('primary', 'escalation'):
                value = phase['data'][child]
                if value['status'] != 'completed' or value.get('carried'):
                    continue
                trace = value['review']['execution']
                self.assertEqual(trace['version'], 'dependency-ready-v1')
                self.assertIsNone(trace['max_concurrent_calls'])
                executions.append({'phase': phase['stage'], 'child': child, 'trace': trace})
        receipt = {'executions': executions, 'synthetic_transport_only': True, 'verifier_revision': VERIFIER_REVISION,
                   'input_hashes': {p: r['sha256'] for p, r in intake['data']['files'].items()},
                   'calls': len(requests), 'prompts': sizes,
                   'reserved_completion_tokens': sum(p['data']['provenance']['allotted_completion_tokens'] for p in phases),
                   'dossier_bytes': len(canonical(record)), 'dossier_sha256': record['sha256'],
                   'intake_sha256': intake['sha256'],
                   'source_coverage': audit_review(phases[-1])['region_coverage'],
                   'score': projection(record)['scores']}
        if directory := os.environ.get('CB_TEST_RECEIPT_DIR'):
            path = Path(directory)
            path.mkdir(parents=True, exist_ok=True)
            (path / (name + '.json')).write_text(json.dumps(receipt, indent=2) + '\n')
        return receipt

    def test_exact_profile_and_shared_resource_contract(self):
        configured = json.loads((ROOT / 'verifier.json').read_text())
        self.assertEqual(configured['judging'], {
            'profile': 'builtin:probabilistic-escalation-v1', 'probability': 1, 'mode': 'override',
            'primary': {'profile': 'builtin:single-v1'},
            'escalation': {'profile': 'builtin:math-decomposed-v1'}})
        effective = load_config(ROOT, 'verifier.json', 'submission')
        self.assertNotIn('max_concurrent_calls', configured['limits'])
        self.assertNotIn('max_concurrent_calls', effective['config']['limits'])
        resolved = effective['judging']['committee']
        self.assertEqual(resolved['escalation']['max_regions'], 8)
        for child in ('primary', 'escalation'):
            for member in resolved[child]['members']:
                self.assertEqual(member['provider'], configured['provider'])
        limits = configured['limits']
        self.assertEqual((limits['total_provider_requests'], limits['total_review_tokens']), (23, 23 * 32768))
        self.assertEqual((limits['total_review_seconds'], limits['max_prompt_bytes'], limits['max_files'],
                          limits['max_followups']), (2400, 2097152, 32, 1))
        self.assertEqual((configured['provider']['max_response_bytes'], configured['provider']['timeout_seconds']),
                         (262144, 1200))

    def test_actual_baseline_primary_accepts_audit_rejects(self):
        intake, initial = self.prepare()
        requests = []
        first = synthetic_review(intake, initial, valid=True, audit_valid=False, blocking=True, requests=requests)
        record = pipeline.dossier(intake, initial, first)
        self.assertIsNone(projection(record)['scores'])
        self.assertFalse(projection(record)['valid'])
        self.assertEqual(first['data']['authority']['judge'], 'escalation')
        self.assertTrue(first['data']['comparison']['valid_disagreement'])
        self.assertEqual(set(audit_review(first)['region_plan']['sources']), set(intake['data']['files']))
        self.measure('baseline', record, requests)

    def test_exact_candidate_full_twenty_three_call_path(self):
        intake, initial = self.prepare(candidate=True)
        requests = []
        first = synthetic_review(intake, initial, request_followup=True, blocking=True,
                                 region_count=8, requests=requests)
        # Both children ask for the same checker, so one batch satisfies both.
        self.assertEqual(first['data']['primary']['review']['decision']['kind'], 'followup')
        self.assertEqual(audit_review(first)['decision']['kind'], 'followup')
        followup = self.followup(intake, initial, first)
        final = synthetic_review(intake, initial, first, followup, valid=False, audit_valid=True,
                                 blocking=True, region_count=8, requests=requests)
        record = pipeline.dossier(intake, initial, first, followup, final)
        self.assertEqual(audit_review(first)['region_plan'], audit_review(final)['region_plan'])
        self.assertEqual(len(audit_review(final)['region_coverage']), 8)
        self.assertEqual(set(audit_review(final)['region_plan']['sources']), set(PROVENANCE['files']))
        self.assertEqual([p['data']['provenance']['attempts'] for p in (first, final)], [11, 12])
        result = self.measure('candidate', record, requests)
        self.assertEqual((result['calls'], result['reserved_completion_tokens'], result['score']),
                         (23, 753664, [3.5, 1]))

    def test_actual_candidate_starts_all_ready_calls_without_four_worker_waves(self):
        intake, initial = self.prepare(candidate=True)
        primary = Event()
        proof = Event()
        # Primary + proof + all eight regions must be admitted concurrently.
        # A hidden four-worker limit or a proof-before-plan wait cannot pass.
        ready = Barrier(10, timeout=10)
        def observe(member, context):
            role = member['id']
            if role == 'reviewer':
                primary.set()
            elif role == 'proof':
                proof.set()
            elif role == 'planner':
                self.assertTrue(primary.wait(10))
                self.assertTrue(proof.wait(10))
            if role in ('reviewer', 'proof', 'region'):
                ready.wait()
        requests = []
        first = synthetic_review(intake, initial, region_count=8,
                                 requests=requests, on_request=observe)
        trace = audit_review(first)['execution']
        self.assertEqual([[call['member_id'] for call in group['calls']]
                          for group in trace['groups']],
                         [['proof'], ['planner'], ['region'] * 8, ['adjudicator']])
        record = pipeline.dossier(intake, initial, first)
        result = self.measure('candidate-all-ready', record, requests)
        self.assertEqual((result['calls'], result['score']), (12, [3.5, 1]))

    def test_only_audit_followup_carries_primary_without_another_primary_call(self):
        intake, initial = self.prepare()
        requests = []
        first = synthetic_review(intake, initial, primary_followup=False, audit_followup=True, requests=requests)
        followup = self.followup(intake, initial, first)
        final = synthetic_review(intake, initial, first, followup, audit_valid=False, requests=requests)
        self.assertTrue(final['data']['primary']['carried'])
        self.assertEqual(sum(request_parts(w)[0]['id'] == 'reviewer' for w in requests), 1)
        self.assertIsNone(projection(pipeline.dossier(intake, initial, first, followup, final))['scores'])

    def test_only_primary_followup_carries_audit_without_repeating_region_work(self):
        intake, initial = self.prepare()
        requests = []
        first = synthetic_review(intake, initial, primary_followup=True, audit_followup=False, requests=requests)
        followup = self.followup(intake, initial, first)
        final = synthetic_review(intake, initial, first, followup, requests=requests)
        self.assertTrue(final['data']['escalation']['carried'])
        self.assertEqual(sum(request_parts(w)[0]['id'] == 'region' for w in requests), 3)
        self.assertEqual(final['data']['provenance']['attempts'], 1)
        self.assertEqual(projection(pipeline.dossier(intake, initial, first, followup, final))['scores'], [3.1, 1])
