"""Exact arithmetic and implication examples, never circuit discoveries."""

import copy
from fractions import Fraction as F
import itertools
import unittest

from expressions import (CONVENTION, ExpressionError, compare, normalize, polynomial,
                         ratio, validate_projection)


def finite(c=1, p=1):
    return {"kind": "finite", "coefficient": ratio(c), "exponent": ratio(p)}


def lane(p=1):
    return {"kind": "lane_endpoint", "exponent": ratio(p)}


GLOBAL = {"kind": "superpolynomial"}
GLOBAL_PROJECTION = {"kind": "global_endpoint"}


class ExpressionTests(unittest.TestCase):
    def test_exact_fractional_exponents_and_small_coefficients(self):
        g = polynomial(F(1, 10**40), F(3, 2))
        self.assertEqual(validate_projection(g, finite(F(1, 10**40), F(3, 2))),
                         [ratio(F(1, 10**40)), ratio(3, 2)])
        self.assertEqual(compare(g, polynomial(100, 1)), "stronger")
        self.assertEqual(compare(polynomial(F(31, 10), 1), polynomial(3, 1)), "stronger")

    def test_zero_negative_and_below_domain(self):
        for c, p in [(0, 10), (-1, 1), (1, F(1, 2))]:
            with self.subTest(c=c, p=p), self.assertRaises(ExpressionError):
                validate_projection(polynomial(c, p), finite(c, p))
        with self.assertRaises(ExpressionError):
            validate_projection(polynomial(1, 1), finite(0, 20))

    def test_reduced_bounded_exact_syntax(self):
        for bad in [{"numerator": "2", "denominator": "2"},
                    {"numerator": "01", "denominator": "1"},
                    {"numerator": "1" * 257, "denominator": "1"},
                    {"numerator": True, "denominator": "1"},
                    {"numerator": "1", "denominator": "0"}]:
            g = polynomial(1, 1)
            g["coefficient"] = bad
            with self.subTest(bad=bad), self.assertRaises(ExpressionError):
                normalize(g)

    def test_unattained_exponent_is_not_rounded_up(self):
        g = polynomial(1, 2, -1)
        self.assertEqual(validate_projection(g, lane(F(3, 2))), ["positive_infinity", ratio(3, 2)])
        validate_projection(g, finite(100, F(199, 100)))
        for projection in [finite(F(1, 10**40), 2), lane(2), finite(1, 3)]:
            with self.subTest(projection=projection), self.assertRaises(ExpressionError):
                validate_projection(g, projection)

    def test_leading_log_factor_and_base(self):
        self.assertEqual(compare(polynomial(1, 2, -1, 1000), polynomial(1, 2, 0, -1)), "weaker")
        self.assertEqual(compare(polynomial(1, 2, 0, 1), lane(2)), "stronger")
        self.assertEqual(compare(lane(2), polynomial(1, 2)), "stronger")
        self.assertEqual(compare(polynomial(1, 2, -1), lane(2)), "weaker")

    def test_sublinear_expressions_are_only_unscored_lemmas(self):
        for g in [polynomial(1, 1, -1), polynomial(1, 1, 0, -1)]:
            normalize(g)
            for projection in [finite(F(1, 10**40), 1), lane(1)]:
                with self.subTest(g=g, projection=projection), self.assertRaises(ExpressionError):
                    validate_projection(g, projection)

    def test_non_negligible_remainder_does_not_raise_coefficient(self):
        with self.assertRaises(ExpressionError):
            validate_projection(polynomial(F(31, 10), F(3, 2)), finite(F(16, 5), F(3, 2)))

    def test_endpoint_precedence(self):
        self.assertEqual(validate_projection(lane(F(3, 2)), lane(F(3, 2))),
                         ["positive_infinity", ratio(3, 2)])
        self.assertEqual(validate_projection(GLOBAL, GLOBAL_PROJECTION),
                         ["positive_infinity", "positive_infinity"])
        for g, p in [(GLOBAL, finite(100, 100)), (GLOBAL, lane(2)),
                     (lane(2), finite(100, 2)), (lane(2), lane(1))]:
            with self.subTest(g=g, p=p), self.assertRaises(ExpressionError):
                validate_projection(g, p)
        with self.assertRaises(ExpressionError):
            validate_projection(polynomial(1, 100), GLOBAL_PROJECTION)

    def test_no_rate_invention(self):
        self.assertEqual(compare(lane(1), polynomial(1, 1, 1)), "weaker")
        self.assertEqual(compare(GLOBAL, polynomial(1, 1000, 1000)), "stronger")
        self.assertEqual(compare(GLOBAL, GLOBAL), "equivalent")

    def test_unknown_is_neither_equality_nor_dominance(self):
        unsupported = {"kind": "exponential", "base": "2"}
        self.assertEqual(compare(unsupported, unsupported), "unknown")
        self.assertEqual(compare(polynomial(1, 1), polynomial(1, 1), right_convention="other"), "unknown")
        hinted = copy.deepcopy(polynomial(1, 1)); hinted["dominates_all"] = True
        self.assertEqual(compare(hinted, GLOBAL), "unknown")
        self.assertEqual(compare(polynomial(1, 1), polynomial(1, 1)), "equivalent")

    def test_transitivity_and_inverse_on_supported_fragment(self):
        guarantees = [polynomial(c, p, q, r) for c, p, q, r in
                      itertools.product([F(1, 10), F(31, 10)], [1, F(3, 2), 2], [-1, 0, 1], [0, 1])]
        guarantees += [lane(1), lane(F(3, 2)), lane(2), GLOBAL]
        relations = {(i, j): compare(a, b) for i, a in enumerate(guarantees) for j, b in enumerate(guarantees)}
        inverse = {"stronger": "weaker", "weaker": "stronger", "equivalent": "equivalent", "unknown": "unknown"}
        for i, j in relations:
            self.assertEqual(relations[j, i], inverse[relations[i, j]])
        for i, j, k in itertools.product(range(len(guarantees)), repeat=3):
            if relations[i, j] in {"stronger", "equivalent"} and relations[j, k] in {"stronger", "equivalent"}:
                self.assertIn(relations[i, k], {"stronger", "equivalent"})
