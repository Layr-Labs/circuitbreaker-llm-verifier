#!/usr/bin/env python3
"""Offline numerical witness. No network, credentials or paid services.
Default verifies embedded Git blobs and replays arithmetic using frozen transcripts.
--full regenerates every numerical transcript used by the shifted-cut checker.
Requires g++/OpenMP, time and memory beyond the hosted Python-slim runtime.
These computations certify finite numerical inequalities, not the analytic reduction.
"""
import argparse, hashlib, os, pathlib, shutil, subprocess, sys, tempfile
# The verifier executes this file in isolated mode beside its evidence mount.
script_directory = pathlib.Path(__file__).resolve().parent
evidence_directory = script_directory / "files"
if not (evidence_directory / "data_a.py").is_file():
    evidence_directory = script_directory
sys.path.insert(0, str(evidence_directory))
from data_a import DATA as A
from data_b import DATA as B

def run_cone(root):
    import runpy

    inputs = root / "cone-inputs"
    inputs.mkdir(exist_ok=True)
    sources = {
        "U5_REFINED.stdout": root / "cubic_certified_two_pass_20261124" / "U5_REFINED.stdout",
        "GAIN.stdout": root / "fair_shift35_20261124" / "GAIN.stdout",
        "MOTIF.stdout": root / "fair_shift35_20261124" / "MOTIF.stdout",
    }
    for name, source in sources.items():
        shutil.copyfile(source, inputs / name)
    previous = sys.argv
    sys.argv = [str(root / "cone-check.py"), str(inputs)]
    try:
        result = runpy.run_path(str(root / "cone-check.py"), run_name="__main__")["result"]
    finally:
        sys.argv = previous
    if result["score_exact"] != "62405759/12405759":
        raise ValueError("incorrect cone score")
    return result

def run():
    parser=argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true")
    args=parser.parse_args()
    if A.keys() & B.keys():
        raise ValueError("duplicate witness path")
    with tempfile.TemporaryDirectory(prefix="circuitbreaker-") as tmp:
        root=pathlib.Path(tmp)
        for name,rec in {**A,**B}.items():
            p=root/name
            if p.resolve().parent!=root.resolve() and root.resolve() not in p.resolve().parents:
                raise ValueError("unsafe path")
            p.parent.mkdir(parents=True,exist_ok=True)
            raw=rec["content"].encode("utf-8")
            got=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
            if got!=rec["git_blob"]:
                raise ValueError("Git blob mismatch: "+name)
            p.write_bytes(raw)
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE="1",OMP_NUM_THREADS="2")
        base=root/"cubic_certified_two_pass_20261124"
        fair=root/"cubic_fair_patch_erasure_20261124"
        motif=root/"fair_edge_motif_20261124"
        shift=root/"fair_shift35_20261124"
        def capture(cmd):
            return subprocess.check_output([str(x) for x in cmd],env=env)
        def same(cmd,expected):
            if capture(cmd)!=expected.read_bytes():
                raise ValueError("Output mismatch: "+str(expected))
        for opt in ([],["-O"]):
            same([sys.executable,*opt,shift/"check.py"],shift/"EXACT.stdout")
        print("PASS: embedded Git blob identities and final rational arithmetic in both Python modes; recorded numerical transcripts used.")
        if not args.full:
            run_cone(root)
            print("Directed C++ computations NOT replayed in this invocation.")
            return
        def build(src,dst,more=()):
            subprocess.run(["g++","-O3","-std=c++17",*more,str(src),"-o",str(dst)],check=True,env=env)
        build(root/"cubic_crossing_line_20261123"/"certify.cpp",root/"moment")
        same([root/"moment",2048,1024],base/"U5_REFINED.stdout")
        build(base/"cert_ray.cpp",root/"gain",["-fopenmp"])
        same([root/"gain",6,-13600,-13600,200,280,100,200,280],fair/"GAIN.stdout")
        same([root/"gain",6,-14000,-14000,200,280,100,200,280],shift/"GAIN.stdout")
        build(fair/"survive.cpp",root/"survive",["-fopenmp"])
        same([root/"survive",-13600,200,300,40,76000,128000],fair/"SURVIVAL.stdout")
        same([root/"survive",-14000,160,240,45,76000,128000],shift/"SURVIVAL.stdout")
        src=(motif/"motif.cpp").read_text()
        if src.count("T=-13600")!=1:
            raise ValueError("unexpected motif threshold")
        shifted=root/"motif_shift.cpp"
        shifted.write_text(src.replace("T=-13600","T=-14000"))
        build(shifted,root/"motif",["-I"+str(motif),"-Wall","-Wextra","-Werror"])
        same([root/"motif"],shift/"MOTIF.stdout")
        for opt in ([],["-O"]):
            same([sys.executable,*opt,shift/"check.py"],shift/"EXACT.stdout")
        run_cone(root)
        print("PASS: full directed moment, gain, survival and motif replays, both final Python modes.")

def hosted(request, data):
    """Small, exact audit. No C++ regeneration or universal proof certification."""
    import contextlib, io, json, runpy
    from itertools import product
    if request != {}:
        raise ValueError('expected empty request object')
    required = {
        'cone-check.py',
        'cubic_certified_two_pass_20261124/mathcheck.py',
        'cubic_certified_two_pass_20261124/U5_REFINED.stdout',
        'cubic_fair_patch_erasure_20261124/bias.py',
        'cubic_fair_patch_erasure_20261124/GAIN.stdout',
        'cubic_fair_patch_erasure_20261124/SURVIVAL.stdout',
        'fair_shift35_20261124/check.py',
        'fair_shift35_20261124/EXACT.stdout',
        'fair_shift35_20261124/GAIN.stdout',
        'fair_shift35_20261124/MOTIF.stdout',
        'fair_shift35_20261124/SURVIVAL.stdout',
    }
    with tempfile.TemporaryDirectory(prefix='cb-audit-') as tmp:
        root = pathlib.Path(tmp)
        for name in sorted(required):
            record = data[name]
            p = root / name
            raw = record['content'].encode('utf-8')
            blob = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
            if blob != record['git_blob']:
                raise ValueError('original Git blob mismatch: ' + name)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(raw)
        shift = root / 'fair_shift35_20261124'
        def replay_arithmetic():
            saved_path = sys.path[:]
            saved_module = sys.modules.pop('mathcheck', None)
            try:
                out = io.StringIO()
                with contextlib.redirect_stdout(out):
                    runpy.run_path(str(shift / 'check.py'), run_name='__main__')
                if out.getvalue().encode('utf-8') != (shift / 'EXACT.stdout').read_bytes():
                    raise ValueError('arithmetic output mismatch')
            finally:
                sys.path[:] = saved_path
                sys.modules.pop('mathcheck', None)
                if saved_module is not None:
                    sys.modules['mathcheck'] = saved_module
        replay_arithmetic()
        mutants = [('GAIN.stdout', '851968', '851967'),
                   ('MOTIF.stdout', '300311352', '300311351'),
                   ('SURVIVAL.stdout', '33422386', '33554432'),
                   ('GAIN.stdout', '-14000', '-13600')]
        rejected = 0
        for name, old, new in mutants:
            f = shift / name
            original = f.read_bytes()
            if old.encode() not in original:
                raise ValueError('mutation target missing')
            try:
                f.write_bytes(original.replace(old.encode(), new.encode(), 1))
                try:
                    replay_arithmetic()
                except ValueError:
                    rejected += 1
                else:
                    raise RuntimeError('receipt mutant passed')
            finally:
                f.write_bytes(original)
        replay_arithmetic()
        with contextlib.redirect_stdout(io.StringIO()):
            cone = run_cone(root)
    def sg(a):
        return (a > 0) - (a < 0)
    def nae(*bits):
        return int(min(bits) != max(bits))
    def least_minimum(values):
        best = min(e for e, _ in values)
        opts = [a for e, a in values if e == best]
        least = tuple(min(a[i] for a in opts) for i in range(len(opts[0])))
        if least not in opts:
            raise ValueError('least minimizer absent')
        return least
    root_cases = 0
    for b in (0, 1):
        for z in product((-1, 0, 1), repeat=4):
            vals = []
            for a in product((0, 1), repeat=5):
                if a[4] > b:
                    continue
                cost = -sum(z[i] * a[i] for i in range(4))
                cost += nae(a[0], a[1], a[4]) + nae(a[2], a[3], a[4])
                vals.append((cost, a))
            m = int(b and sg(z[0] + z[1]) + sg(z[2] + z[3]) > 0)
            e = tuple(int(z[i] + sg(z[i ^ 1] + (sg(z[i ^ 2] + z[i ^ 3]) if b else -1)) > 0) for i in range(4))
            if (*e, m) != least_minimum(vals):
                raise ValueError('root message table mismatch')
            root_cases += 1
    paths = 0
    for word, target in [('0011100', '001100'), ('00100', '0000')]:
        pins, expected = tuple(map(int, word)), tuple(map(int, target))
        for l, r in product((0, 1), repeat=2):
            vals = []
            for a in product((0, 1), repeat=len(pins)-1):
                q = (l, *a, r)
                vals.append((sum(nae(q[i], q[i+1], b) for i, b in enumerate(pins)), a))
            if least_minimum(vals) != expected:
                raise ValueError('motif guard mismatch')
            paths += 1
    return {'all_passed': True, 'source_blobs_checked': len(required),
            'root_least_minimizer_cases': root_cases, 'guarded_motif_cases': paths,
            'receipt_mutants_rejected': rejected,
            'coefficient': '62405759/12405759',
            'cone_faces': cone['A3']['face_cases'] + cone['A5']['face_cases'],
            'cone_cost_slack_lower': cone['cost_slack_lower'],
            'scope': 'Exact arithmetic from frozen numerical transcripts and local finite identities only. Does not regenerate the directed numerical transcripts or verify the universal mathematical proof.'}

if __name__=="__main__":
    import json
    # Official experiment protocol supplies a JSON request on stdin.
    # Offline invocation without input or with --full retains the original CLI.
    raw = sys.stdin.read() if len(sys.argv) == 1 and not sys.stdin.isatty() else ""
    if raw.strip():
        print(json.dumps(hosted(json.loads(raw), {**A, **B}), sort_keys=True))
    else:
        run()
