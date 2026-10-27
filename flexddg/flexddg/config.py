'''target.json: one file per target directory, validated on load, with defaults for everything optional.

    {"title": "CLDN18.2 (9v32)", "pdb": "../../9v32_complex.pdb",
     "groups": {"fab": ["H", "L"], "antigen": ["A"]},
     "move": "antigen",
     "chain_labels": {"H": "Heavy chain", "L": "Light chain"},
     "scan": {"side": ["fab"]}}

groups   named sets of chains; every chain of interest must be in exactly one group
move     the group pulled away to make the unbound state (Rosetta's chains_to_move)
scan     interface alanine scan: residues of the `side` groups within `cutoff` A (heavy atoms) of
         the other side; mode "none" disables it so that only `mutations` are run
mutations  explicit mutants, e.g. "H:Y33F" or "H:Y33A+L:Y100A" (insertion code before the new
         residue: "H:S52AA" is S52A -> A); run in addition to the scan
'''
import copy, json, os, re

DEFAULTS = {
    'title': None,
    'pdb': None,
    'groups': None,
    'move': None,
    'chain_labels': {},
    'scan': {'mode': 'interface', 'side': None, 'cutoff': 5.0, 'to': 'A',
             'skip': ['ALA', 'GLY', 'PCA', 'disulfide']},
    'mutations': [],
    'protocol': {'nstruct': 35, 'ntrials': 35000, 'stride': 35000, 'maxmin': 5000, 'conv': 1.0,
                 'extra_flags': []},
    'slurm': {'partition': 'normal,possu', 'time': '04:00:00', 'mem': '3G', 'max_concurrent': 100,
              'analysis_time': '01:00:00', 'analysis_mem': '8G'},
    'analysis': {'bootstrap': 1000, 'seed': 0, 'hotspot': 1.0, 'min_nstruct': None},
    'rosetta_module': 'biology rosetta/3.14',
}
REQUIRED = ['pdb', 'groups', 'move']
AA1 = set('ACDEFGHIKLMNPQRSTVWY')
MUT_RE = re.compile(r'^([A-Za-z0-9]):([A-Z])(-?\d+)([A-Z]?)([A-Z])$')


class ConfigError(ValueError):
    pass


def _merge(defaults, given, where):
    out = copy.deepcopy(defaults)
    for k, v in given.items():
        if k not in defaults:
            raise ConfigError('%s: unknown key %r (allowed: %s)' % (where, k, ', '.join(defaults)))
        if isinstance(defaults[k], dict) and defaults[k] and k != 'chain_labels':
            if not isinstance(v, dict):
                raise ConfigError('%s.%s must be an object' % (where, k))
            out[k] = _merge(defaults[k], v, '%s.%s' % (where, k))
        else:
            out[k] = v
    return out


def parse_mutation(spec):
    '''"H:Y33A+L:Y100A" -> [(chain, wt, resnum, icode, mut), ...]'''
    out = []
    for part in spec.split('+'):
        m = MUT_RE.match(part.strip())
        if not m:
            raise ConfigError('bad mutation %r: expected CHAIN:WT<num>[icode]MUT, e.g. H:Y33A' % part)
        chain, wt, num, icode, mut = m.groups()
        if wt not in AA1 or mut not in AA1:
            raise ConfigError('bad mutation %r: unknown amino acid' % part)
        out.append((chain, wt, int(num), icode or ' ', mut))
    if len({(c, n, i) for c, _, n, i, _ in out}) != len(out):
        raise ConfigError('bad mutation %r: a position appears twice' % spec)
    return out


def mutation_name(muts):
    '''[(H, Y, 33, ' ', A)] -> "H_Y33A"; several positions are joined with "-".'''
    return '-'.join('%s_%s%d%s%s' % (c, wt, n, i.strip(), mut) for c, wt, n, i, mut in muts)


def validate(cfg):
    for k in REQUIRED:
        if cfg[k] in (None, '', {}, []):
            raise ConfigError('missing required key %r' % k)
    groups = cfg['groups']
    if not isinstance(groups, dict) or len(groups) < 2:
        raise ConfigError('groups must map at least two names to chain lists')
    seen = {}
    for g, chains in groups.items():
        if not isinstance(chains, list) or not chains:
            raise ConfigError('groups.%s must be a non-empty list of chain IDs' % g)
        for c in chains:
            if not (isinstance(c, str) and len(c) == 1):
                raise ConfigError('groups.%s: chain IDs are single characters, got %r' % (g, c))
            if c in seen:
                raise ConfigError('chain %s is in both groups %s and %s' % (c, seen[c], g))
            seen[c] = g
    if cfg['move'] not in groups:
        raise ConfigError('move=%r is not a group (groups: %s)' % (cfg['move'], ', '.join(groups)))
    scan = cfg['scan']
    if scan['mode'] not in ('interface', 'none'):
        raise ConfigError('scan.mode must be "interface" or "none"')
    if scan['side'] is None:
        scan['side'] = list(groups)
    if isinstance(scan['side'], str):
        scan['side'] = [scan['side']]
    for g in scan['side']:
        if g not in groups:
            raise ConfigError('scan.side: %r is not a group' % g)
    if scan['to'] not in AA1:
        raise ConfigError('scan.to must be a one-letter amino acid')
    for c in cfg['chain_labels']:
        if c not in seen:
            raise ConfigError('chain_labels: chain %s is not in any group' % c)
    for spec in cfg['mutations']:
        for c, *_ in parse_mutation(spec):
            if c not in seen:
                raise ConfigError('mutation %s: chain %s is not in any group' % (spec, c))
    p = cfg['protocol']
    for k in ('nstruct', 'ntrials', 'stride', 'maxmin'):
        if not (isinstance(p[k], int) and p[k] > 0):
            raise ConfigError('protocol.%s must be a positive integer' % k)
    if p['ntrials'] % p['stride']:
        raise ConfigError('protocol.ntrials must be a multiple of protocol.stride')
    if cfg['analysis']['min_nstruct'] is None:
        cfg['analysis']['min_nstruct'] = p['nstruct']
    if cfg['title'] is None:
        cfg['title'] = os.path.basename(cfg['dir'])
    return cfg


def load(target_dir):
    target_dir = os.path.abspath(target_dir)
    path = os.path.join(target_dir, 'target.json')
    if not os.path.exists(path):
        raise ConfigError('%s not found' % path)
    try:
        given = json.load(open(path))
    except json.JSONDecodeError as e:
        raise ConfigError('%s: %s' % (path, e))
    cfg = _merge(DEFAULTS, given, 'target.json')
    cfg['dir'] = target_dir
    cfg = validate(cfg)
    cfg['pdb'] = os.path.normpath(os.path.join(target_dir, cfg['pdb']))
    return cfg


def chain_order(cfg):
    return [c for chains in cfg['groups'].values() for c in chains]


def group_of(cfg, chain):
    for g, chains in cfg['groups'].items():
        if chain in chains:
            return g
