"""Generate v2/mixed-version parity fixtures through the actual exporter.

Every phase is explicitly synthetic with a stubbed HTTP boundary. Expected
ordering comes from hand-assigned tiers, never expressions.compare.
"""
import argparse
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from export_metrics import projection
from expressions import CONVENTION, claims_for, powers, ratio
from llm_verifier.errors import VerifierError
from tests.remainder_cases import expected_relation, ranked_cases, term
from tests.synthetic import synthetic_dossier
from tools.generate_hosted import VERIFIER_REVISION

FIXTURE = ROOT / 'tests/fixtures/circuitbreaker-frontier-v2.json'


def projected(g):
    if g['kind'] == 'superpolynomial':
        return {'kind': 'global_endpoint'}
    if g['kind'] == 'lane_endpoint':
        return g
    leading = g['leading'] if g['kind'] == 'eventual_lower_bound' else g
    p, q, r = powers(leading)
    if (q, r) > (0, 0):
        return {'kind': 'lane_endpoint', 'exponent': ratio(p)}
    if (q, r) < (0, 0):
        return {'kind': 'lane_endpoint', 'exponent': ratio((p + 1) / 2)}
    return {'kind': 'finite', 'coefficient': leading['coefficient'], 'exponent': leading['exponent']}


def generate():
    cases = []
    entries = ranked_cases()
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / 'submission').mkdir()
        for name in ('verifier.json', 'rules.md'):
            shutil.copyfile(ROOT / name, root / name)
        shutil.copytree(ROOT / 'checkers', root / 'checkers')

        def emit(entry, *, valid=True, lemma=False, simulation=False):
            bound = {'version': entry['version'], 'convention': CONVENTION,
                     'family': 'SYNTHETIC remainder transport family ' + entry['id'] + '; no circuit construction asserted.',
                     'theorem': 'SYNTHETIC test only; no mathematical review or discovery.',
                     'guarantee': entry['guarantee'], 'projection': projected(entry['guarantee'])}
            manifest = {'version': 1, 'title': 'SYNTHETIC remainder case ' + entry['id'],
                        'claims': ['SYNTHETIC unscored lemma.'] if lemma else claims_for(bound),
                        'experiments': [], 'checks': []}
            dossier = synthetic_dossier(root, bound, manifest, valid=valid, simulation=simulation)
            try:
                output = projection(dossier)
            except VerifierError as error:
                if not simulation:
                    raise
                assert str(error) == 'simulation_cannot_publish_score'
                return None
            if 'metrics' in output:
                output['metrics'] = {'circuitbreaker_frontier': output['metrics']['circuitbreaker_frontier']}
            return output

        for entry in entries:
            cases.append({'id': entry['id'], 'provenanceKind': 'exporter', 'consumerState': 'accepted',
                          'expectedEligible': True, 'strengthTier': entry['rank'], 'output': emit(entry)})
        for name, options in [('rejected', {'valid': False}), ('lemma', {'lemma': True}), ('simulation', {'simulation': True})]:
            cases.append({'id': name, 'provenanceKind': 'exporter', 'consumerState': name,
                          'expectedEligible': False, 'output': emit(entries[3], **options)})

    explicit_sample = next(c['output'] for c in cases if c['id'] == 'explicit-52-four-fifths')
    def bound(m):
        return m['resources']['circuit_size']['bound']
    mutations = {
        'v2-bound-in-v1-envelope': lambda m: m.update(version='circuitbreaker-frontier-v1', comparator='circuitbreaker-v1'),
        'v2-version-v1-comparator': lambda m: m.update(comparator='circuitbreaker-v1'),
        'v1-version-v2-comparator': lambda m: m.update(version='circuitbreaker-frontier-v1'),
        'unknown-version': lambda m: m.update(version='circuitbreaker-frontier-v3'),
        'wrong-scope': lambda m: m.update(scope='other-family-model'),
        'extra-envelope-field': lambda m: m.update(dominates_all=True),
        'wrong-source-hash': lambda m: m['resources']['circuit_size'].update(source_sha256='bad'),
        'missing-bound-evidence': lambda m: m['resources']['circuit_size'].update(evidence={'proof.md': '1'*64}),
        'unsafe-evidence-path': lambda m: m['resources']['circuit_size'].update(evidence={'bound.json': '0'*64, 'proof.md': '1'*64, '../bad': '2'*64}),
        'extra-resource': lambda m: m['resources'].update(other=copy.deepcopy(m['resources']['circuit_size'])),
        'non-subleading-remainder': lambda m: bound(m)['remainder'].update(terms=[term(1, 1)]),
        'zero-remainder-coefficient': lambda m: bound(m)['remainder'].update(terms=[term(0)]),
        'negative-remainder-coefficient': lambda m: bound(m)['remainder'].update(terms=[term(-1)]),
        'duplicate-remainder-powers': lambda m: bound(m)['remainder'].update(terms=[term(1), term(2)]),
        'ascending-remainder-powers': lambda m: bound(m)['remainder'].update(terms=[term(1), term(1, F(1,2))]),
        'too-many-remainder-terms': lambda m: bound(m)['remainder'].update(terms=[term(1,F(i,10)) for i in range(9,0,-1)]),
        'decaying-remainder': lambda m: bound(m)['remainder'].update(terms=[term(1,-1)]),
        'nonreduced-remainder': lambda m: bound(m)['remainder']['terms'][0].update(coefficient={'numerator':'2','denominator':'2'}),
        'big-o-with-constant': lambda m: bound(m).update(remainder={'kind':'big_o','coefficient':ratio(1),'exponent':ratio(0),'log_exponent':ratio(1),'loglog_exponent':ratio(0)}),
        'candidate-comparison-hint': lambda m: bound(m).update(compare_as='strongest'),
        'invalid-leading-kind': lambda m: bound(m).update(leading={'kind':'superpolynomial'}),
        'unspecified-remainder': lambda m: bound(m).update(remainder=None),
    }
    for name, mutate in mutations.items():
        output = copy.deepcopy(explicit_sample)
        mutate(output['metrics']['circuitbreaker_frontier'])
        cases.append({'id': name, 'provenanceKind': 'mutated-exporter', 'consumerState': 'accepted',
                      'expectedEligible': False, 'output': output})
    for state in ('unreviewed', 'revoked', 'stale'):
        cases.append({'id': state, 'provenanceKind': 'consumer-state', 'consumerState': state,
                      'expectedEligible': False, 'output': copy.deepcopy(explicit_sample)})
    eligible = [c for c in cases if c['expectedEligible']]
    def scenario(name, ids):
        selected = [c for c in eligible if c['id'] in ids]
        best = max(c['strengthTier'] for c in selected)
        return {'id': name, 'caseIds': ids, 'expectedFrontierIds': [c['id'] for c in selected if c['strengthTier'] == best]}
    contract = {'id':'complexity', 'kind':'metrics-frontier', 'label':'Complexity',
                'metricsPath':'circuitbreaker_frontier', 'metricsVersion':'circuitbreaker-frontier-v2',
                'comparator':'circuitbreaker-v2', 'scope':CONVENTION, 'model':'full-b2-v1', 'notation':'lower-bound',
                'columns':[{'key':'circuit_size','label':'Circuit size lower bound','direction':'+'}],
                'defaultSort':{'by':'metric','key':'circuit_size','direction':'desc'}}
    return {'contract':contract,
            'provenance':{'synthetic':True,'generator':'tools/generate_remainder_fixtures.py',
                          'verifierRevision':VERIFIER_REVISION,
                          'description':'Synthetic deterministic HTTP-stubbed dossiers through actual exporter. Hand-assigned real-envelope tiers; no scientific acceptance or integer-rounding equivalence asserted.'},
            'cases':cases, 'expectedEligibleIds':[c['id'] for c in eligible],
            'expectedFrontierIds':['global','v2-global'],
            'comparisons':[{'left':a['id'],'right':b['id'],'expected':expected_relation(a,b)} for a in entries for b in entries],
            'frontierScenarios':[
                scenario('historical-refinements', ['legacy-four','explicit-182-four-fifths','explicit-52-four-fifths','explicit-nine-hundredth','big-o-log']),
                scenario('same-rate-constant', ['big-o-log','explicit-ten-log','explicit-nine-log']),
                scenario('legacy-version-tie', ['legacy-four','v2-liminf-four']),
                scenario('explicit-multiterm', ['explicit-nine-hundredth','explicit-multiterm-23','explicit-multiterm-22']),
                scenario('equivalent-records-retained', ['eventual-exact','eventual-exact-another-family']),
                scenario('null-score-endpoints', ['eventual-nlog','eventual-quadratic-log-loss','lane-two']),
            ]}


def encoded():
    return (json.dumps(generate(), ensure_ascii=True, separators=(',', ':'), allow_nan=False) + '\n').encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    content = encoded()
    if args.check:
        if FIXTURE.read_bytes() != content:
            raise SystemExit('Remainder fixture drift; regenerate with the pinned verifier')
        print('V2 remainder fixtures reproduce exactly')
    else:
        FIXTURE.write_bytes(content)
        print(f'Wrote {FIXTURE.relative_to(ROOT)} ({len(content)} bytes)')


if __name__ == '__main__':
    main()
