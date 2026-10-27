'''Run one flex ddG trajectory and classify existing ones.'''
import json, os, re, shutil, subprocess, sys

PKG = os.path.dirname(os.path.abspath(__file__))
XML = os.path.join(PKG, 'data', 'ddG-backrub.xml')
ROSETTA = 'rosetta_scripts.default.linuxgccrelease'
_trials_re = re.compile(r'number_backrub_trials=(\d+)')


def traj_dir(cfg, mut, i, base=None):
    return os.path.join(base or cfg['dir'], 'output', mut, '%02d' % i)


def traj_state(cfg, d):
    '''"done", "missing", "incomplete" (started, no success line), or "mismatch" (ran other ntrials).'''
    out = os.path.join(d, 'rosetta.out')
    if not os.path.exists(out):
        return 'missing'
    text = open(out, errors='replace').read()
    if 'reported success' not in text or not os.path.exists(os.path.join(d, 'ddG.db3')):
        return 'incomplete'
    m = _trials_re.search(text)
    if m and int(m.group(1)) != cfg['protocol']['ntrials']:
        return 'mismatch'
    return 'done'


def command(cfg, input_dir, i, protocol):
    manifest = json.load(open(os.path.join(input_dir, 'input.json')))
    return [ROSETTA,
            '-s', os.path.join(input_dir, manifest['pdb']),
            '-parser:protocol', XML,
            '-parser:script_vars', 'chainstomove=' + manifest['chains_to_move'],
            'mutate_resfile_relpath=' + os.path.join(input_dir, manifest['resfile']),
            'number_backrub_trials=%d' % protocol['ntrials'], 'max_minimization_iter=%d' % protocol['maxmin'],
            'abs_score_convergence_thresh=%s' % protocol['conv'],
            'backrub_trajectory_stride=%d' % protocol['stride'],
            '-restore_talaris_behavior', '-in:file:fullatom', '-ignore_unrecognized_res',
            '-ignore_zero_occupancy', 'false', '-ex1', '-ex2', '-constant_seed', '-jran', str(1000 + i)
            ] + list(protocol['extra_flags'])


def run_trajectory(cfg, mut, i, protocol=None, base=None):
    '''Run trajectory i of mut into <base>/output/<mut>/<ii>; skip if already done. Returns exit code.'''
    protocol = protocol or cfg['protocol']
    input_dir = os.path.join(cfg['dir'], 'inputs', mut)
    out = traj_dir(cfg, mut, i, base)
    if traj_state(dict(cfg, protocol=protocol), out) == 'done':
        print('done: %s' % out)
        return 0
    if shutil.which(ROSETTA) is None:
        sys.exit('%s not on PATH; load it first: ml %s' % (ROSETTA, cfg['rosetta_module']))
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    with open(os.path.join(out, 'rosetta.out'), 'w') as f:
        rc = subprocess.call(command(cfg, input_dir, i, protocol), cwd=out, stdout=f, stderr=subprocess.STDOUT)
    print('%s %d: exit %d -> %s' % (mut, i, rc, traj_state(dict(cfg, protocol=protocol), out)))
    return rc
