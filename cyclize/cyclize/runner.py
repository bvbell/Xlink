'''Run one variant's closure job and classify finished ones.

One array task = one variant, running all nstruct in a single Rosetta process. Rosetta's own nstruct
loop is much cheaper than a process per structure, and it keeps each task comfortably above the
10-minute floor that makes a Slurm job worth submitting.

Partial failure is normal and is NOT an error: GenKIC closure yield is well under 100% for constrained
rings, so Rosetta exits non-zero whenever any single structure failed to close. The number that
matters is how many PDBs landed, so `state` counts output rather than trusting the exit code.
'''
import glob, json, os, shutil, subprocess, sys

PKG = os.path.dirname(os.path.abspath(__file__))
ROSETTA = 'rosetta_scripts.default.linuxgccrelease'


def variant_dir(cfg, vid, base=None):
    return os.path.join(base or cfg['dir'], 'output', vid)


def closed_count(d):
    return len(glob.glob(os.path.join(d, '*.pdb')))


def state(cfg, vid, base=None):
    '''"done" (nstruct reached), "partial", "empty" (ran, nothing closed), or "missing".'''
    d = variant_dir(cfg, vid, base)
    log = os.path.join(d, 'rosetta.out')
    if not os.path.exists(log):
        return 'missing'
    n = closed_count(d)
    if n >= cfg['protocol']['nstruct']:
        return 'done'
    return 'partial' if n else 'empty'


def command(cfg, input_dir, out_dir, seed):
    manifest = json.load(open(os.path.join(input_dir, 'input.json')))
    p = cfg['protocol']
    return [ROSETTA,
            '-s', os.path.normpath(os.path.join(input_dir, manifest['pdb'])),
            '-parser:protocol', os.path.join(input_dir, 'protocol.xml'),
            '-nstruct', str(p['nstruct']),
            '-out:path:pdb', out_dir, '-out:prefix', '',
            '-overwrite',
            # Mirror-symmetric scoring, required for a fair comparison of D and L turn residues.
            '-score:symmetric_gly_tables',
            # Without this the S-S CONECT record is dropped and the ring looks open downstream.
            '-write_all_connect_info',
            # SG atoms are far apart before closure, so autodetection would mis-pair or miss them;
            # DeclareBond/ForceDisulfides in the protocol establish the bond explicitly.
            '-in:detect_disulf', 'false',
            '-ignore_unrecognized_res', '-ignore_zero_occupancy', 'false',
            '-constant_seed', '-jran', str(1000 + seed),
            ] + list(p['extra_flags'])


def run_variant(cfg, vid, seed=1, base=None):
    '''Run all nstruct of one variant into <base>/output/<vid>/. Returns 0 unless nothing closed.'''
    input_dir = os.path.join(cfg['dir'], 'inputs', vid)
    if not os.path.isdir(input_dir):
        sys.exit('%s has no inputs; run `cyclize setup` first' % vid)
    out = variant_dir(cfg, vid, base)
    if state(cfg, vid, base) == 'done':
        print('done: %s' % out)
        return 0
    if shutil.which(ROSETTA) is None:
        sys.exit('%s not on PATH; load it first: ml %s' % (ROSETTA, cfg['rosetta_module']))
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    with open(os.path.join(out, 'rosetta.out'), 'w') as f:
        rc = subprocess.call(command(cfg, input_dir, out, seed), cwd=out, stdout=f,
                             stderr=subprocess.STDOUT)
    n = closed_count(out)
    print('%s: rosetta exit %d, %d/%d structures closed -> %s'
          % (vid, rc, n, cfg['protocol']['nstruct'], state(cfg, vid, base)))
    # Non-zero exit just means some structures did not close, which is expected; only a total failure
    # is worth reporting as one.
    return 0 if n else 1
