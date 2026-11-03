"""CircuitBreaker guarantee algebra. No claims of mathematical acceptance.

Legacy polynomial/log guarantees mean liminf C(n)/g(n) >= 1. V2 also
preserves explicit eventual inequalities and quantified Big-O losses. Comparisons
order the strength of guarantees, even for different eligible families.
They neither identify families nor infer upper bounds on circuit complexity.
"""

from fractions import Fraction
import json
import re


VERSION = "circuitbreaker-expression-v1"
VERSION_V2 = "circuitbreaker-expression-v2"
CONVENTION = "uniform-np-full-b2-actual-length-all-lengths-v1"
REVIEW = "circuitbreaker-joint-bounds-v3"
POLICIES = {
    "scoring": {"version": "coefficient-exponent-expression-v2"},
    "validity": {"version": "fixed-family-full-b2-v1"},
}


class ExpressionError(ValueError):
    pass


def fields(value, expected):
    if type(value) is not dict or set(value) != set(expected):
        raise ExpressionError("unexpected expression fields")


def rational(value):
    fields(value, ("numerator", "denominator"))
    n, d = value["numerator"], value["denominator"]
    if (type(n) is not str or type(d) is not str
            or not 1 <= len(n) <= 256 or not 1 <= len(d) <= 256
            or not re.fullmatch(r"0|-?[1-9][0-9]*", n)
            or not re.fullmatch(r"[1-9][0-9]*", d)):
        raise ExpressionError("bounded canonical rational required")
    result = Fraction(int(n), int(d))
    if str(result.numerator) != n or str(result.denominator) != d:
        raise ExpressionError("reduced rational required")
    return result


def ratio(n, d=1):
    f = Fraction(n, d)
    return {"numerator": str(f.numerator), "denominator": str(f.denominator)}


def normalize(guarantee, *, version=None):
    """Validate canonical AST; never silently rewrite a submitted claim."""
    if version not in (None, VERSION, VERSION_V2):
        raise ExpressionError("unsupported expression version")
    if type(guarantee) is not dict:
        raise ExpressionError("guarantee object required")
    kind = guarantee.get("kind")
    if version == VERSION and kind == "eventual_lower_bound":
        raise ExpressionError("eventual lower bounds require expression v2")
    if kind == "polynomial_log":
        fields(guarantee, ("kind", "coefficient", "exponent",
                           "log_exponent", "loglog_exponent"))
        if rational(guarantee["coefficient"]) <= 0:
            raise ExpressionError("positive coefficient required")
        if rational(guarantee["exponent"]) < 1:
            raise ExpressionError("exponent at least one required")
        rational(guarantee["log_exponent"])
        rational(guarantee["loglog_exponent"])
    elif kind == "lane_endpoint":
        fields(guarantee, ("kind", "exponent"))
        if rational(guarantee["exponent"]) < 1:
            raise ExpressionError("exponent at least one required")
    elif kind == "superpolynomial":
        fields(guarantee, ("kind",))
    elif kind == "eventual_lower_bound":
        fields(guarantee, ("kind", "leading", "remainder"))
        leading = guarantee["leading"]
        if type(leading) is not dict or leading.get("kind") != "polynomial_log":
            raise ExpressionError("polynomial/log leading term required")
        normalize(leading)
        remainder = guarantee["remainder"]
        if type(remainder) is not dict:
            raise ExpressionError("remainder object required")
        remainder_kind = remainder.get("kind")
        if remainder_kind == "explicit":
            fields(remainder, ("kind", "terms"))
            terms = remainder["terms"]
            if type(terms) is not list or len(terms) > 8:
                raise ExpressionError("at most eight explicit remainder terms required")
            previous = powers(leading)
            for term in terms:
                fields(term, ("coefficient", "exponent", "log_exponent", "loglog_exponent"))
                if rational(term["coefficient"]) <= 0:
                    raise ExpressionError("positive remainder coefficient required")
                term_powers = powers(term)
                if not (Fraction(0), Fraction(0), Fraction(0)) <= term_powers < previous:
                    raise ExpressionError("nondecaying strictly descending subleading terms required")
                previous = term_powers
        elif remainder_kind == "big_o":
            fields(remainder, ("kind", "exponent", "log_exponent", "loglog_exponent"))
            if not (Fraction(0), Fraction(0), Fraction(0)) <= powers(remainder) < powers(leading):
                raise ExpressionError("nondecaying subleading Big-O rate required")
        else:
            raise ExpressionError("unsupported remainder kind")
    else:
        raise ExpressionError("unsupported guarantee kind")
    # A detached canonical copy; no candidate comparison keys survive.
    return json.loads(json.dumps(guarantee, sort_keys=True, allow_nan=False))


def powers(g):
    return tuple(rational(g[k]) for k in
                 ("exponent", "log_exponent", "loglog_exponent"))


def implies(left, right):
    """Sound sufficient implication rules for these lower-bound guarantees."""
    a, b = left["kind"], right["kind"]
    if a == "eventual_lower_bound" or b == "eventual_lower_bound":
        left_leading = left["leading"] if a == "eventual_lower_bound" else left
        right_leading = right["leading"] if b == "eventual_lower_bound" else right
        if left_leading != right_leading:
            # Strict leading dominance supplies enough slack for every
            # supported lower-order loss. A lane theorem has unbounded
            # coefficient at its lane, so its polynomial implications also
            # supply strict slack.
            return implies(left_leading, right_leading)
        if b != "eventual_lower_bound":
            return True
        if a != "eventual_lower_bound":
            # Liminf supplies no uniform subleading error rate.
            return False
        return remainder_implies(left["remainder"], right["remainder"])
    if a == "superpolynomial":
        return True
    if b == "superpolynomial":
        return False
    if a == b == "polynomial_log":
        return (powers(left) > powers(right) or
                (powers(left) == powers(right) and
                 rational(left["coefficient"]) >= rational(right["coefficient"])))
    if a == b == "lane_endpoint":
        return rational(left["exponent"]) >= rational(right["exponent"])
    if a == "polynomial_log":
        return powers(left) > (rational(right["exponent"]), Fraction(0), Fraction(0))
    # Divergence at p implies a polynomial/log bound below n^p, and
    # every finite coefficient at exactly n^p, but no rate of divergence.
    return (rational(left["exponent"]), Fraction(0), Fraction(0)) >= powers(right)


def remainder_implies(left, right):
    """Whether the left loss guarantees at least the right eventual bound."""
    a, b = left["kind"], right["kind"]
    if a == b == "big_o":
        return powers(left) <= powers(right)
    if a == "explicit" and b == "big_o":
        return not left["terms"] or powers(left["terms"][0]) <= powers(right)
    if a == "big_o" and b == "explicit":
        return bool(right["terms"]) and powers(left) < powers(right["terms"][0])
    # The eventual sign of a finite sum is the sign of its largest nonzero
    # coefficient in the polynomial/log/loglog asymptotic scale.
    coefficients = {}
    for sign, remainder in ((1, left), (-1, right)):
        for term in remainder["terms"]:
            key = powers(term)
            coefficients[key] = coefficients.get(key, Fraction(0)) + sign * rational(term["coefficient"])
    for key in sorted(coefficients, reverse=True):
        if coefficients[key]:
            return coefficients[key] < 0
    return True


def compare(left, right, *, left_convention=CONVENTION, right_convention=CONVENTION,
            left_version=None, right_version=None):
    """Return stronger/weaker/equivalent/unknown. No candidate ordering hints."""
    if left_convention != CONVENTION or right_convention != CONVENTION:
        return "unknown"
    try:
        left, right = normalize(left, version=left_version), normalize(right, version=right_version)
    except (ExpressionError, TypeError, ValueError):
        return "unknown"
    ab, ba = implies(left, right), implies(right, left)
    if ab and ba:
        return "equivalent"
    if ab:
        return "stronger"
    if ba:
        return "weaker"
    return "unknown"


def polynomial(c, p, q=0, r=0):
    return {"kind": "polynomial_log", "coefficient": ratio(c), "exponent": ratio(p),
            "log_exponent": ratio(q), "loglog_exponent": ratio(r)}


def validate_projection(guarantee, projection):
    guarantee = normalize(guarantee)
    if type(projection) is not dict:
        raise ExpressionError("projection required")
    kind = projection.get("kind")
    if guarantee["kind"] == "superpolynomial" and kind != "global_endpoint":
        raise ExpressionError("explicit global endpoint requires global projection")
    if guarantee["kind"] == "lane_endpoint" and projection != guarantee:
        raise ExpressionError("explicit lane endpoint requires its exact lane projection")
    if kind == "finite":
        fields(projection, ("kind", "coefficient", "exponent"))
        c, p = rational(projection["coefficient"]), rational(projection["exponent"])
        if c <= 0 or p < 1:
            raise ExpressionError("finite projection requires c > 0 and p >= 1")
        target = polynomial(c, p)
        scores = [projection["coefficient"], projection["exponent"]]
    elif kind == "lane_endpoint":
        target = normalize(projection)
        scores = ["positive_infinity", projection["exponent"]]
    elif kind == "global_endpoint":
        fields(projection, ("kind",))
        target = {"kind": "superpolynomial"}
        # Transport placeholders for ONE theorem, not n^infinity.
        scores = ["positive_infinity", "positive_infinity"]
    else:
        raise ExpressionError("unsupported projection kind")
    if not implies(guarantee, target):
        raise ExpressionError("projection is not implied by the full guarantee")
    return scores


def validate_bound(bound):
    fields(bound, ("version", "convention", "family", "theorem", "guarantee", "projection"))
    if bound["version"] not in (VERSION, VERSION_V2) or bound["convention"] != CONVENTION:
        raise ExpressionError("unsupported bound convention/version")
    if (bound["version"] == VERSION and type(bound["guarantee"]) is dict
            and bound["guarantee"].get("kind") == "eventual_lower_bound"):
        raise ExpressionError("eventual lower bounds require expression v2")
    for key, limit in (("family", 2000), ("theorem", 2000)):
        if type(bound[key]) is not str or not bound[key].strip() or len(bound[key]) > limit:
            raise ExpressionError("bounded family and theorem text required")
    return validate_projection(bound["guarantee"], bound["projection"])


def joint_statement(bound):
    validate_bound(bound)
    kind = bound["guarantee"]["kind"]
    meanings = {
        "polynomial_log": (
            "For g(n)=c*n^p*(log2 n)^q*(log2 log2 n)^r with the exact parameters "
            "in guarantee, liminf over all integer n of C_B2(f_n)/g(n) >= 1. "
            "Equivalently, for every epsilon>0 there is N_epsilon such that "
            "C_B2(f_n)>=(1-epsilon)*g(n) for every integer n>=N_epsilon. "
        ),
        "lane_endpoint": (
            "For the fixed exponent p in guarantee, for every real M>0 there "
            "is N_M such that C_B2(f_n)>=M*n^p for every integer n>=N_M. "
        ),
        "superpolynomial": (
            "For every fixed real p>=1 and every real M>0 there is N_(p,M) "
            "such that C_B2(f_n)>=M*n^p for every integer n>=N_(p,M). "
            "The family and its uniform NP verifier do not depend on p or M. "
        ),
    }
    if kind == "eventual_lower_bound":
        remainder = bound["guarantee"]["remainder"]
        meaning = (
            "Let L(n)=c*n^p*(log2 n)^q*(log2 log2 n)^r with the exact parameters "
            "in guarantee.leading. "
        )
        if remainder["kind"] == "explicit":
            meaning += (
                "Let R(n) be the sum of the monomials in guarantee.remainder.terms, "
                "each coefficient times n^exponent*(log2 n)^log_exponent*"
                "(log2 log2 n)^loglog_exponent; an empty sum is zero. "
                "There exists N such that C_B2(f_n)>=L(n)-R(n) for every integer "
                "n>=N. All remainder coefficients are the exact declared constants. "
            )
        else:
            meaning += (
                "Let h(n)=n^exponent*(log2 n)^log_exponent*"
                "(log2 log2 n)^loglog_exponent with the exact parameters in "
                "guarantee.remainder. There exist a fixed real K>0 and N such "
                "that C_B2(f_n)>=L(n)-K*h(n) for every integer n>=N. "
                "K is existential, not a specified coefficient. "
            )
    else:
        meaning = meanings[kind]
    return (
        "Joint CircuitBreaker theorem under rules.md: the following declaration "
        "defines ONE fixed eligible total uniform NP family in the pinned full B2 "
        "model at actual input length n. Its full guarantee and its declared "
        "projection both hold unconditionally, for all sufficiently large integer "
        "lengths with the quantifiers specified by the theorem kinds. Both score "
        "source claims assert this ENTIRE declaration, not separate lemmas. "
        + meaning + "Declaration: " + json.dumps(bound, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    )


def claims_for(bound, evidence=None):
    """Authoring convenience. Review and export independently check its output."""
    scores = validate_bound(bound)
    evidence = evidence if evidence is not None else ["bound.json", "proof.md"]
    statement = joint_statement(bound)
    if len(statement) > 8000:
        raise ExpressionError("joint statement exceeds verifier claim limit")
    return [{"statement": statement, "score": score, "evidence": list(evidence)}
            for score in scores]
