"""Exact certificate for the UNCHANGED Li–Yang LP; no theorem acceptance."""
from fractions import Fraction as F

POINT = {"alpha_I": F(48, 5), "alpha_Q": F(9, 5), "alpha_phi": F(1, 5), "delta": F(64, 5)}


def terms(point):
    i, q, phi = (point[x] for x in ("alpha_I", "alpha_Q", "alpha_phi"))
    return (i / 3, 2 - 2 * phi + q, 4 - 4 * phi, 3 + phi, 5 - q, (5 - 2 * phi + q) / 2)


def feasible(point):
    if set(point) != set(POINT) or any(type(v) is not F for v in point.values()):
        return False
    i, q, phi, delta = (point[x] for x in ("alpha_I", "alpha_Q", "alpha_phi", "delta"))
    return 0 <= phi <= F(1, 2) and 0 <= q <= i and delta >= 0 and all(delta - i <= t for t in terms(point))


def coefficient(point):
    if not feasible(point):
        raise ValueError("not a feasible exact point of the unchanged LP")
    return point["delta"] - point["alpha_I"] - point["alpha_phi"] / 2


def certificate():
    # Rows represent u <= constant + phi_coeff*phi + q_coeff*Q.
    rows = ((F(3), F(1), F(0)), (F(5), F(0), F(-1)), (F(5, 2), F(-1), F(1, 2)))
    weights = (F(7, 10), F(1, 10), F(1, 5))
    weighted = tuple(sum(w * row[j] for w, row in zip(weights, rows)) for j in range(3))
    assert sum(weights) == 1 and all(w >= 0 for w in weights)
    assert weighted == (F(31, 10), F(1, 2), F(0))
    assert feasible(POINT) and coefficient(POINT) == F(31, 10)
    return {"weights": weights, "weighted_rhs": weighted, "optimum": F(31, 10), "terms": terms(POINT)}


def run(inputs):
    """Optional protected checker: exact point/certificate, never a circuit proof."""
    import re
    if not isinstance(inputs, dict) or set(inputs) != set(POINT):
        raise ValueError("provide exactly alpha_I, alpha_Q, alpha_phi, delta")
    point={}
    for key,value in inputs.items():
        if not isinstance(value,str) or len(value)>64 or not re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?",value):
            raise ValueError("parameters must be bounded canonical rational strings")
        point[key]=F(value)
        if str(point[key])!=value:
            raise ValueError("parameters must be reduced canonical rational strings")
    cert=certificate()
    u=point['delta']-point['alpha_I']
    good=feasible(point)
    return {'passed':good,'coefficient':str(u-point['alpha_phi']/2),
        'constraint_slacks':[str(t-u) for t in terms(point)],
        'dual_weights':[str(w) for w in cert['weights']], 'dual_bound':'31/10',
        'message':('Feasible exact point and unchanged-LP upper certificate; no circuit theorem certified.' if good else 'Infeasible point of the unchanged LP; no theorem or score certified.')}


if __name__=='__main__':
    import json,sys
    try:
        result=run(json.load(sys.stdin))
    except (ValueError,TypeError,ZeroDivisionError,RecursionError) as exc:
        result={'passed':False,'coefficient':'0','constraint_slacks':[],
          'dual_weights':['7/10','1/10','1/5'],'dual_bound':'31/10','message':str(exc)[:1000]}
    print(json.dumps(result,separators=(',',':'),allow_nan=False))
