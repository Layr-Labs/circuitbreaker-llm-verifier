"""Hand-ordered synthetic envelopes. These are not accepted circuit results."""
from fractions import Fraction as F
from expressions import VERSION, VERSION_V2, polynomial, ratio


def term(c=1, p=0, q=0, r=0):
    return {k: v for k, v in polynomial(c, p, q, r).items() if k != 'kind'}


def explicit(*terms, leading=None):
    return {'kind': 'eventual_lower_bound', 'leading': leading or polynomial(4, 1),
            'remainder': {'kind': 'explicit', 'terms': list(terms)}}


def big_o(p=0, q=0, r=0, *, leading=None):
    return {'kind': 'eventual_lower_bound', 'leading': leading or polynomial(4, 1),
            'remainder': {'kind': 'big_o', 'exponent': ratio(p),
                          'log_exponent': ratio(q), 'loglog_exponent': ratio(r)}}


def ranked_cases():
    # Rankings are assigned by hand, independently of the implementation.
    rows = [
        ('legacy-31/10', 0, polynomial(F(31, 10), 1), VERSION),
        ('legacy-four', 1, polynomial(4, 1), VERSION),
        ('v2-liminf-four', 1, polynomial(4, 1), VERSION_V2),
        ('big-o-four-fifths', 2, big_o(F(4, 5)), VERSION_V2),
        ('explicit-182-four-fifths', 3, explicit(term(182, F(4, 5))), VERSION_V2),
        ('explicit-52-four-fifths', 4, explicit(term(52, F(4, 5))), VERSION_V2),
        ('big-o-sqrt', 5, big_o(F(1, 2)), VERSION_V2),
        ('explicit-nine-hundredth', 6, explicit(term(9, F(1, 100))), VERSION_V2),
        ('explicit-multiterm-23', 7, explicit(term(8, F(1, 100)), term(6, F(3, 50000)), term(23)), VERSION_V2),
        ('explicit-multiterm-22', 8, explicit(term(8, F(1, 100)), term(6, F(3, 50000)), term(22)), VERSION_V2),
        ('big-o-log', 9, big_o(0, 1), VERSION_V2),
        ('explicit-ten-log', 10, explicit(term(10, 0, 1)), VERSION_V2),
        ('explicit-nine-log-plus-loglog', 11, explicit(term(9, 0, 1), term(3, 0, 0, 1)), VERSION_V2),
        ('explicit-nine-log', 12, explicit(term(9, 0, 1)), VERSION_V2),
        ('big-o-loglog', 13, big_o(0, 0, 1), VERSION_V2),
        ('explicit-loglog', 14, explicit(term(1, 0, 0, 1)), VERSION_V2),
        ('big-o-constant', 15, big_o(), VERSION_V2),
        ('explicit-three', 16, explicit(term(3)), VERSION_V2),
        ('explicit-two', 17, explicit(term(2)), VERSION_V2),
        ('explicit-three-quarters', 18, explicit(term(F(3, 4))), VERSION_V2),
        ('explicit-quarter', 19, explicit(term(F(1, 4))), VERSION_V2),
        ('eventual-exact', 20, explicit(), VERSION_V2),
        ('eventual-exact-another-family', 20, explicit(), VERSION_V2),
        ('legacy-larger-coefficient', 21, polynomial(4 + F(1, 10**40), 1), VERSION),
        ('lane-one', 22, {'kind': 'lane_endpoint', 'exponent': ratio(1)}, VERSION),
        ('eventual-nlog', 23, explicit(term(100, 1), leading=polynomial(1, 1, 1)), VERSION_V2),
        ('legacy-quadratic-log-loss', 24, polynomial(1, 2, -1), VERSION),
        ('eventual-quadratic-log-loss', 25, big_o(1, leading=polynomial(1, 2, -1)), VERSION_V2),
        ('lane-two', 26, {'kind': 'lane_endpoint', 'exponent': ratio(2)}, VERSION),
        ('global', 27, {'kind': 'superpolynomial'}, VERSION),
        ('v2-global', 27, {'kind': 'superpolynomial'}, VERSION_V2),
    ]
    return [{'id': name, 'rank': rank, 'guarantee': guarantee, 'version': version}
            for name, rank, guarantee, version in rows]


def expected_relation(a, b):
    return 'stronger' if a['rank'] > b['rank'] else 'weaker' if a['rank'] < b['rank'] else 'equivalent'
