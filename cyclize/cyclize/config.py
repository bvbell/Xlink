'''cyclize.json: one file per target directory, validated on load, defaults for everything optional.

    {"title": "APJ d006 CDR3 macrocycles",
     "complex": "../../apj_complex.pdb",
     "motif_chain": "A",
     "receptor_chain": "B",
     "scaffolds": {"short": "../../apj_scaffolds/short.pdb"},
     "scan": "../../flexddg/apj_nb/analysis_output/alascan_summary.csv"}

complex         the parent complex the motif came from; validation superimposes designs back onto it
motif_chain     chain the scaffolds were excised from (numbering must still match the complex)
receptor_chain  the chain a closed design must not clash with
scaffolds       name -> excised motif PDB. Names become the prefix of every variant id.
scan            flex ddG alascan_summary.csv, read for per-residue ddG_GAM; hotspots are protected
                from being released to GenKIC and from being mutated to Cys
'''
import copy, json, os

DEFAULTS = {
    'title': None,
    'complex': None,
    'motif_chain': 'A',
    'receptor_chain': None,
    'scaffolds': None,
    'scan': None,
    'hotspot': 1.0,
    'routes': ['amide', 'ssapp', 'ssmut'],
    # Sweep every legal identity at i+1 of the turn instead of using the five hand-picked
    # TURN_SEEDS. Only the amide route has an i+1 to sweep; ssapp and ssmut are unaffected.
    'redesign_i1': False,
    'protocol': {'nstruct': 50, 'closure_attempts': 2000, 'min_solutions': 1,
                 'scorefxn': 'ref2015', 'relax_rounds': 1, 'design_i2': False, 'rama_cutoff': None, 'seed_rama_cutoff': None, 'coord_dev': 0.5,
                 'omega_trans': True, 'omega_weight': None,
                 # target_loop was read by inputs.py but missing here, so _merge rejected any
                 # cyclize.json that set it -- the option was unreachable. See variants.TARGET_LOOP.
                 'target_loop': 5,
                 'extra_flags': []},
    # clash 2.2 A, not 3.0: the native APJ interface packs motif residues to 2.44-2.80 A of the
    # receptor (F104 2.44, E108 2.65, K110 2.70), so a 3.0 A cutoff rejects the native structure
    # itself. 2.2 A is below anything observed natively, so it flags only real overlap.
    'validate': {'clash': 2.2, 'motif_rmsd': 1.0},
    # per_task: variants per array task. One variant runs in about a minute, so a task per variant
    # is mostly queue latency and is poor scheduler etiquette; 20 gives tasks of roughly 20 minutes.
    'slurm': {'partition': 'normal,possu', 'time': '02:00:00', 'mem': '2G', 'max_concurrent': 100,
              'per_task': 20, 'analysis_time': '00:30:00', 'analysis_mem': '4G'},
    'rosetta_module': 'biology rosetta/3.14',
}
REQUIRED = ['complex', 'scaffolds']


class ConfigError(ValueError):
    pass


def _merge(defaults, given, where):
    out = copy.deepcopy(defaults)
    for k, v in given.items():
        if k not in defaults:
            raise ConfigError('%s: unknown key %r (allowed: %s)' % (where, k, ', '.join(defaults)))
        if isinstance(defaults[k], dict) and defaults[k] and k != 'scaffolds':
            if not isinstance(v, dict):
                raise ConfigError('%s.%s must be an object' % (where, k))
            out[k] = _merge(defaults[k], v, '%s.%s' % (where, k))
        else:
            out[k] = v
    return out


def validate(cfg):
    for k in REQUIRED:
        if cfg[k] in (None, '', {}, []):
            raise ConfigError('missing required key %r' % k)
    if not isinstance(cfg['scaffolds'], dict):
        raise ConfigError('scaffolds must map a name to a PDB path')
    for name in cfg['scaffolds']:
        if '__' in name or not name:
            raise ConfigError('scaffold name %r: "__" separates fields in variant ids, so it is '
                              'reserved' % name)
    from .variants import ROUTES
    for r in cfg['routes']:
        if r not in ROUTES:
            raise ConfigError('routes: %r is not one of %s' % (r, ', '.join(ROUTES)))
    for c in ('motif_chain', 'receptor_chain'):
        if cfg[c] is not None and len(cfg[c]) != 1:
            raise ConfigError('%s must be a single character, got %r' % (c, cfg[c]))
    p = cfg['protocol']
    for k in ('nstruct', 'closure_attempts', 'min_solutions'):
        if not (isinstance(p[k], int) and p[k] > 0):
            raise ConfigError('protocol.%s must be a positive integer' % k)
    if not isinstance(cfg['hotspot'], (int, float)) or cfg['hotspot'] <= 0:
        raise ConfigError('hotspot must be a positive ddG cutoff in kcal/mol')
    if cfg['title'] is None:
        cfg['title'] = os.path.basename(cfg['dir'])
    return cfg


def load(target_dir):
    target_dir = os.path.abspath(target_dir)
    path = os.path.join(target_dir, 'cyclize.json')
    if not os.path.exists(path):
        raise ConfigError('%s not found' % path)
    try:
        given = json.load(open(path))
    except json.JSONDecodeError as e:
        raise ConfigError('%s: %s' % (path, e))
    cfg = _merge(DEFAULTS, given, 'cyclize.json')
    cfg['dir'] = target_dir
    cfg = validate(cfg)
    for k in ('complex', 'scan'):
        if cfg[k]:
            cfg[k] = os.path.normpath(os.path.join(target_dir, cfg[k]))
    cfg['scaffolds'] = {n: os.path.normpath(os.path.join(target_dir, p))
                        for n, p in cfg['scaffolds'].items()}
    return cfg


def read_scan(path, column='ddG_GAM'):
    '''alascan_summary.csv -> {resnum: ddG}. Case ids look like "A_F104A".'''
    import csv
    if not path:
        return {}
    if not os.path.exists(path):
        raise ConfigError('scan file %s not found' % path)
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            case = row['case']
            try:
                num = int(case.split('_', 1)[1][1:-1])
            except (IndexError, ValueError):
                raise ConfigError('%s: cannot read a residue number out of case %r' % (path, case))
            out[num] = float(row[column])
    return out
