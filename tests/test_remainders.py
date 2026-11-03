"""Exact remainder discrimination and compatibility, without provider calls."""
from fractions import Fraction as F
import itertools
import hashlib
import unittest

from expressions import (CONVENTION, VERSION, VERSION_V2, ExpressionError,
                         claims_for, compare, joint_statement, normalize,
                         polynomial, ratio, validate_bound, validate_projection)
from tests.remainder_cases import big_o, explicit, expected_relation, ranked_cases, term
from tests.test_expressions import finite, lane


class RemainderTests(unittest.TestCase):
    def test_hand_ordered_matrix(self):
        cases = ranked_cases()
        for a, b in itertools.product(cases, repeat=2):
            with self.subTest(left=a['id'], right=b['id']):
                self.assertEqual(compare(a['guarantee'], b['guarantee']), expected_relation(a, b))

    def test_smaller_losses_separate_identical_numeric_projections(self):
        for g in [explicit(term(182, F(4, 5))), explicit(term(52, F(4, 5))),
                  explicit(term(9, F(1, 100))), big_o(0, 1), explicit()]:
            self.assertEqual(validate_projection(g, finite(4, 1)), [ratio(4), ratio(1)])
        self.assertEqual(compare(big_o(0, 1), explicit(term(9, F(1, 100)))), 'stronger')

    def test_unknown_big_o_constant_is_never_an_explicit_constant(self):
        self.assertEqual(compare(big_o(0, 1), explicit(term(10**40, 0, 1))), 'weaker')
        self.assertEqual(compare(big_o(0, 1), explicit(term(F(1, 10**40), F(1, 100)))), 'stronger')
        self.assertEqual(compare(big_o(), explicit()), 'weaker')

    def test_lexicographic_remainder_domain(self):
        for g in [explicit(term(1, F(1, 2), -100)), big_o(F(1, 2), -100),
                  explicit(term(1, 0, 1, -100)), big_o(0, 0, 1)]:
            self.assertEqual(normalize(g), g)
        for g in [explicit(term(1, -1, 100)), big_o(0, -1, 100), big_o(0, 0, -1)]:
            with self.assertRaises(ExpressionError):
                normalize(g)

    def test_rejects_malformed_noncanonical_or_non_subleading_losses(self):
        bad = [explicit(term(0)), explicit(term(-1)), explicit(term(1, 1)), big_o(1),
               explicit(term(1), term(1, F(1, 2))), explicit(term(1), term(2)),
               explicit(*[term(1, F(i, 10)) for i in range(9, 0, -1)])]
        for value in [None, [], True, {'kind': 'O'}, {'kind': 'explicit', 'terms': {}},
                      {'kind': 'big_o', 'coefficient': ratio(1), 'exponent': ratio(0),
                       'log_exponent': ratio(1), 'loglog_exponent': ratio(0)}]:
            g = explicit(); g['remainder'] = value; bad.append(g)
        for value in [None, polynomial(0, 1), {'kind': 'superpolynomial'}, explicit()]:
            g = explicit(); g['leading'] = value; bad.append(g)
        g = explicit(term(1)); g['remainder']['terms'][0]['coefficient'] = {'numerator': '2', 'denominator': '2'}; bad.append(g)
        g = explicit(); g['comparison_hint'] = 'best'; bad.append(g)
        g = explicit(); g['remainder']['threshold'] = '10'; bad.append(g)
        g = explicit(term(1)); g['remainder']['terms'][0]['exponent'] = ratio(0); g['remainder']['terms'][0]['exponent']['numerator'] = '9' * 257; bad.append(g)
        for g in bad:
            with self.subTest(g=g):
                with self.assertRaises(ExpressionError):
                    normalize(g)
                self.assertEqual(compare(g, g), 'unknown')

    def test_version_gate_and_joint_quantifiers(self):
        bound = {'version': VERSION_V2, 'convention': CONVENTION, 'family': 'Synthetic family.',
                 'theorem': 'Synthetic quantified bound.', 'guarantee': explicit(term(9, F(1, 100))),
                 'projection': finite(4)}
        validate_bound(bound)
        text = joint_statement(bound)
        self.assertIn('There exists N', text)
        self.assertIn('exact declared constants', text)
        self.assertNotIn('liminf', text)
        self.assertEqual(claims_for(bound)[0]['statement'], claims_for(bound)[1]['statement'])
        bound['guarantee'] = big_o(0, 1)
        self.assertIn('fixed real K>0', joint_statement(bound))
        self.assertEqual(compare(bound['guarantee'], bound['guarantee'], left_version=VERSION), 'unknown')
        self.assertEqual(compare(bound['guarantee'], bound['guarantee'], right_version='unknown'), 'unknown')
        bound['version'] = VERSION
        with self.assertRaisesRegex(ExpressionError, 'require expression v2'):
            validate_bound(bound)

    def test_legacy_statement_bytes_unchanged(self):
        # SHA256 of statements independently produced by immutable base
        # 2cd3a17f3f296b7f0ac248a63a9b54a0fabdf851, without requiring Git at runtime.
        hashes = ['805587581336e7ac02182d7929426a4926887e5cb5d1cd2f9cfdb9873f42b6e4', '7bd4949238e68e6710b55f231e568eddd4ac56a7d6b416ae586b468f9286ed3c', 'e9e895e7f32f7610182a8abd8806fa234d26a0ebb0afe441ab0188166afe1414']
        for i, (guarantee, projection) in enumerate([(polynomial(4, 1), finite(4)), (lane(2), lane(2)),
                                       ({'kind': 'superpolynomial'}, {'kind': 'global_endpoint'})]):
            bound = {'version': VERSION, 'convention': CONVENTION, 'family': 'Synthetic family.',
                     'theorem': 'Synthetic theorem.', 'guarantee': guarantee, 'projection': projection}
            self.assertEqual(hashlib.sha256(joint_statement(bound).encode()).hexdigest(), hashes[i])

    def test_endpoint_and_unattained_projection_rules_survive(self):
        guarantee = big_o(1, leading=polynomial(1, 2, -1))
        validate_projection(guarantee, lane(F(3, 2)))
        for invalid in [finite(F(1, 10**40), 2), lane(2), {'kind': 'global_endpoint'}]:
            with self.assertRaises(ExpressionError):
                validate_projection(guarantee, invalid)
        self.assertEqual(compare(lane(2), guarantee), 'stronger')
        self.assertEqual(compare(guarantee, lane(F(3, 2))), 'stronger')

    def test_inverse_transitivity_and_detached_normalization(self):
        guarantees = [row['guarantee'] for row in ranked_cases()]
        guarantees += [explicit(term(c, p, q)) for c, p, q in
                       itertools.product([F(1, 3), 7], [F(1, 100), F(4, 5)], [-1, 0, 1])]
        relations = {(i, j): compare(a, b) for i, a in enumerate(guarantees) for j, b in enumerate(guarantees)}
        inverse = {'stronger': 'weaker', 'weaker': 'stronger', 'equivalent': 'equivalent', 'unknown': 'unknown'}
        for i, j in relations:
            self.assertEqual(relations[j, i], inverse[relations[i, j]])
        for i, j, k in itertools.product(range(len(guarantees)), repeat=3):
            if relations[i, j] in {'stronger', 'equivalent'} and relations[j, k] in {'stronger', 'equivalent'}:
                self.assertIn(relations[i, k], {'stronger', 'equivalent'})
        g = explicit(term(1)); detached = normalize(g)
        detached['remainder']['terms'].clear()
        self.assertTrue(g['remainder']['terms'])
        self.assertEqual(compare(g, g, right_convention='another-model'), 'unknown')


if __name__ == '__main__':
    unittest.main()
