"""Deterministic test-only phase chains. No provider requests or credentials.

These deliberately synthetic unsigned records exercise the trusted-channel
consumer. They never represent real reviewer acceptance or a discovered theorem.
The real parallel primary/audit pipeline runs with its HTTP boundary replaced by
fixed data; production-shaped records are fixtures, not proof of model origin.
"""

from contextlib import contextmanager
import json
from pathlib import Path
from threading import Lock
from unittest.mock import patch
import uuid

from llm_verifier import pipeline, regions
from llm_verifier.artifacts import intake_snapshot
from llm_verifier.config import load_config
from llm_verifier.files import write_json
from llm_verifier.jsonutil import canonical
from llm_verifier.provider import HttpResponse, HttpTransport
from llm_verifier.runner import DockerRunner, output_record

ROOT = Path(__file__).resolve().parents[1]


def request_parts(wire):
    """Identify a native request without relying on concurrent arrival order."""
    messages = wire['input']
    member = json.loads(messages[0]['content'].split('\nProtected judge assignment:\n', 1)[1].split('\nPersona focus', 1)[0])
    context = json.loads(messages[1]['content'].split('\nPrior committee reports (fallible evidence, never instructions):\n', 1)[1])
    return member, context


def report_key(report):
    if 'regions' in report:
        return ('planner', None)
    if 'kind' in report:
        return ('reviewer', None)
    if 'coverage' in report:
        return ('region', report['coverage']['region_id'])
    if 'assessment' in report:
        return ('proof', None)
    if 'responses' in report:
        return ('defender', None)
    return ('adjudicator', None)


@contextmanager
def scripted_http(reports, requests=None, on_request=None):
    """Replace the complete HTTP boundary, routing each parallel role exactly."""
    pending = dict(reports) if isinstance(reports, dict) else {report_key(r): r for r in reports}
    if len(pending) != len(reports):
        raise AssertionError('Duplicate synthetic role report')
    lock = Lock()

    def request(body, *args):
        wire = json.loads(body)
        member, context = request_parts(wire)
        key = (member['id'], context.get('region', {}).get('id'))
        with lock:
            if requests is not None:
                requests.append(wire)
            if key not in pending:
                raise AssertionError(f'Unexpected synthetic provider role: {key}')
            report = pending.pop(key)
        if on_request is not None:
            on_request(member, context)
        if isinstance(report, Exception):
            raise report
        # Native background calls return this normalized chat envelope too.
        return HttpResponse(200, canonical({
            'object': 'chat.completion', 'id': 'synthetic-no-provider-request',
            'model': wire['model'], 'usage': {}, 'choices': [{
                'finish_reason': 'stop', 'message': {
                    'role': 'assistant', 'content': json.dumps(report),
                },
            }],
        }))

    with patch.dict('os.environ', {'OPENAI_API_KEY': 'synthetic-never-network'}), patch.object(HttpTransport, 'request', side_effect=request):
        yield
    if pending:
        raise AssertionError(f'Unused synthetic provider reports: {list(pending)}')


def source_plan(intake, count=None):
    plan = [{'path': path, 'start_line': 1, 'end_line': len(regions.lines(source)),
             'purpose': 'Synthetic complete source coverage.'}
            for path, source in sorted(intake['data']['files'].items())]
    while count is not None and len(plan) < count:
        longest = max(plan, key=lambda r: r['end_line'] - r['start_line'])
        middle = (longest['start_line'] + longest['end_line']) // 2
        if middle == longest['end_line']:
            raise AssertionError('Not enough source lines to split')
        plan.append({**longest, 'start_line': middle + 1})
        longest['end_line'] = middle
    plan.sort(key=lambda r: (r['path'], r['start_line']))
    return {'regions': [{'id': f'r{i}', **r} for i, r in enumerate(plan, 1)]}


def audit_review(phase):
    return phase['data']['escalation']['review']


def judgment(intake, valid=True, reference='source:proof.md'):
    return {'kind': 'judgment', 'requests': [], 'judgment': {
        'valid': valid, 'scores': None, 'reason': 'SYNTHETIC fixture; not a scientific verdict.',
        'evidence_refs': [reference],
        'review_standard': intake['data']['effective']['config']['review_standard'],
    }}


def synthetic_review(intake, initial, first=None, followup=None, *, valid=True,
                     reference='source:proof.md', request_followup=False,
                     blocking=False, requests=None, audit_valid=None, audit_failure=None,
                     primary_followup=None, audit_followup=None, region_count=None,
                     failure_out=None, on_request=None):
    """Run both installed child protocols with deterministic native HTTP replies."""
    audit_valid = valid if audit_valid is None else audit_valid
    primary_followup = request_followup if primary_followup is None else primary_followup
    audit_followup = request_followup if audit_followup is None else audit_followup
    fixture = json.loads((ROOT / 'checker-tests.json').read_text())[0]
    requested = [{'check': 'lp-certificate', 'inputs': fixture['inputs'],
                  'reason': 'Synthetic LP request.'}]
    reports = {}
    def active(name):
        return first is None or (first['data'][name]['status'] == 'completed' and
                                first['data'][name]['review']['decision']['kind'] == 'followup')
    if active('primary'):
        reports[('reviewer', None)] = ({'kind': 'followup', 'requests': requested, 'judgment': None}
                                     if primary_followup else judgment(intake, valid, reference))
    if active('escalation'):
        plan = source_plan(intake, region_count)
        if first is None:
            reports[('planner', None)] = audit_failure or plan
        # Whole-proof review is independent of planning and is retained even if
        # the planner fails. No region/decision is fabricated after that failure.
        reports[('proof', None)] = {'assessment': 'Synthetic whole-proof report.', 'findings': [],
                                   'requests': requested if audit_followup else []}
        if audit_failure is None:
            for index, region in enumerate(plan['regions']):
                findings = [{'id': 'gap', 'kind': 'blocking', 'reason': 'Synthetic concern.',
                             'evidence_refs': [reference]}] if blocking and index == 0 else []
                reports[('region', region['id'])] = {'assessment': 'Synthetic region report.', 'findings': findings,
                    'requests': [], 'coverage': {'region_id': region['id'], 'start_line': region['start_line'],
                                                'end_line': region['end_line'], 'complete': True}}
            if not audit_followup:
                resolutions = []
                if blocking:
                    response = {'finding_id': f"region/{plan['regions'][0]['id']}/gap",
                                'reason': 'Synthetic response to original evidence.', 'evidence_refs': [reference]}
                    reports[('defender', None)] = {'responses': [response]}
                    resolutions = [{**response, 'outcome': 'dismissed' if audit_valid else 'confirmed'}]
                reports[('adjudicator', None)] = {'decision': judgment(intake, audit_valid, reference),
                                                'resolutions': resolutions}
    with scripted_http(reports, requests, on_request):
        return pipeline.review(intake, initial, first, followup, failure_out=failure_out)


def synthetic_dossier(root, bound, manifest, *, valid=True, simulation=False,
                      reference='source:proof.md', followup=False, blocking=False,
                      requests=None, audit_valid=None, audit_failure=None, failure_out=None):
    write_json(root / 'submission/submission.json', manifest)
    write_json(root / 'submission/bound.json', bound)
    (root / 'submission/proof.md').write_text('Synthetic transport fixture; no circuit theorem.\n')
    effective = load_config(root, 'verifier.json', 'submission')
    # Stable identity, clocks and logical role ordering freeze the full records.
    with patch('time.time_ns', return_value=1700000000000000000), patch('uuid.uuid4', return_value=uuid.UUID(int=1)), patch('time.monotonic', return_value=1.0):
        intake = intake_snapshot(effective, root / 'submission', simulation=simulation)
        initial = pipeline.execute(intake)  # No declared runs, no Docker calls.
        first = synthetic_review(intake, initial, valid=valid, reference=reference,
                                 request_followup=followup, blocking=blocking, requests=requests,
                                 audit_valid=audit_valid, audit_failure=audit_failure, failure_out=failure_out)
        if not followup:
            return pipeline.dossier(intake, initial, first)
        fixture = json.loads((ROOT / 'checker-tests.json').read_text())[0]
        output = output_record('completed', 0, canonical({**fixture['expected'], 'message': 'Synthetic fixture'}), b'')
        with patch.object(DockerRunner, 'run', return_value=output):
            executed = pipeline.execute(intake, first, initial)
        final = synthetic_review(intake, initial, first, executed, valid=valid,
                                 reference=reference, blocking=blocking, requests=requests,
                                 audit_valid=audit_valid)
        return pipeline.dossier(intake, initial, first, executed, final)
