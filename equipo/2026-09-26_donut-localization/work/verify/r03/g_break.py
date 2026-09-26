# -*- coding: utf-8 -*-
"""Try to break scripts/check_provenance.py with scratch copies (--root); real paper/ untouched."""
import os, shutil, subprocess, sys, json
H = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(H, *['..']*5))
S = os.path.join(H, 'scratch'); shutil.rmtree(S, ignore_errors=True)
GOOD = {"crb_center_lg_L50_N100_nm": {"statement": "centre CRB", "type": "script", "reproduce": "scripts/compute_paper_numbers.py",
        "number_keys": ["crb_center_lg_L50_N100_nm"]}}
def mk(name, body, prov=GOOD, extra=None, numbers_tex_mod=None):
    r = os.path.join(S, name)
    shutil.copytree(os.path.join(ROOT, 'paper'), os.path.join(r, 'paper'), ignore=shutil.ignore_patterns('figures'))
    os.makedirs(os.path.join(r, 'data')); shutil.copy(os.path.join(ROOT, 'data', 'paper_numbers.json'), os.path.join(r, 'data'))
    shutil.copytree(os.path.join(ROOT, 'structure'), os.path.join(r, 'structure'))
    shutil.copytree(os.path.join(ROOT, 'scripts'), os.path.join(r, 'scripts'), ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copytree(os.path.join(ROOT, 'tests'), os.path.join(r, 'tests'), ignore=shutil.ignore_patterns('__pycache__'))
    os.makedirs(os.path.join(r, 'paper', 'sections'), exist_ok=True)
    open(os.path.join(r, 'paper', 'sections', 'body.tex'), 'w', encoding='utf-8').write(body)
    m = open(os.path.join(r, 'paper', 'main.tex'), encoding='utf-8').read().replace('% \input{sections/introduction}', '\input{sections/body}')
    open(os.path.join(r, 'paper', 'main.tex'), 'w', encoding='utf-8').write(m)
    json.dump(prov, open(os.path.join(r, 'paper', 'provenance.json'), 'w'))
    if extra: extra(r)
    if numbers_tex_mod:
        p = os.path.join(r, 'paper', 'generated', 'numbers.tex'); t = open(p, encoding='utf-8').read()
        open(p, 'w', encoding='utf-8').write(numbers_tex_mod(t))
    return r
def run(name, expect_fail, **kw):
    r = mk(name, **kw)
    pr = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'check_provenance.py'), '--root', r], capture_output=True, text=True, cwd=r)
    failed = pr.returncode != 0
    verdict = 'as expected' if failed == expect_fail else '*** CHECKER HOLE / FALSE ALARM ***'
    print('[%s] exit=%d expect_fail=%s -> %s' % (name, pr.returncode, expect_fail, verdict))
    print('   ' + '\n   '.join((pr.stdout + pr.stderr).strip().splitlines()[-4:]))
K = 'crb_center_lg_L50_N100_nm'
run('01_good', False, body='CRB \pnum{%s}\src{%s}.\n' % (K, K))
run('02_unknown_pnum', True, body='CRB \pnum{no_such_key}\src{%s}.\n' % K)
run('03_unresolved_src', True, body='CRB \pnum{%s}\src{no_such_entry}.\n' % K)
run('04_pnumse_without_se', True, body='\pnumse{%s}\src{%s}\n' % (K, K))
run('05_bad_number_keys', True, body='\pnum{%s}\src{%s}\n' % (K, K), prov={K: dict(GOOD[K], number_keys=['bogus_key'])})
run('06_missing_reproduce_file', True, body='\pnum{%s}\src{%s}\n' % (K, K), prov={K: dict(GOOD[K], reproduce='scripts/nope.py')})
run('07_numbers_tex_out_of_sync', True, body='\pnum{%s}\src{%s}\n' % (K, K), numbers_tex_mod=lambda t: t.replace('{1.605}', '{1.705}', 1))
run('08_commented_bad_pnum_ignored', False, body='%% \pnum{no_such_key}\n\pnum{%s}\src{%s}\n' % (K, K))
run('09_escaped_percent_then_bad', True, body='50\%% \pnum{no_such_key}\src{%s}\n' % K)
run('10_input_without_braces', True, body='\pnum{%s}\src{%s}\n\input sections/hidden\n' % (K, K),
    extra=lambda r: open(os.path.join(r, 'paper', 'sections', 'hidden.tex'), 'w').write('\pnum{no_such_key}\src{no_such_entry}\n'))
run('11_src_with_space', True, body='\pnum{%s}\src {no_such_entry}\n' % K)
run('12_pnum_via_macro', True, body='\newcommand{\nn}[1]{\pnum{#1}}\nn{no_such_key}\src{%s}\n' % K)
run('13_missing_input_file', True, body='\input{sections/doesnotexist}\n')
run('14_claim_bad_number', True, body='\pnum{%s}\src{%s}\n' % (K, K),
    extra=lambda r: json.dump(json.load(open(os.path.join(r, 'structure', 'claims.json'), encoding='utf-8')) + [{"key": "x", "numbers": ["bogus"], "script": "scripts/compute_paper_numbers.py"}],
                              open(os.path.join(r, 'structure', 'claims.json'), 'w', encoding='utf-8')))
run('15_input_subdir_relative_nested', True, body='\input{sections/sub}\n',
    extra=lambda r: (open(os.path.join(r, 'paper', 'sections', 'sub.tex'), 'w').write('\input{sections/leaf}\n'),
                     open(os.path.join(r, 'paper', 'sections', 'leaf.tex'), 'w').write('\pnum{bogus2}\n')))
run('16_circular_input', True, body='\input{sections/body}\n')
