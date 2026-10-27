'''Set up a target directory: clean PDB, interface list, one input directory per mutation.

Writes into the target dir:
    <id>_clean.pdb          (pdbprep)
    interface_residues.txt  interface scan residues kept, in chain order
    chain_breaks.txt        backbone gaps Rosetta will treat as bonded (informational)
    inputs/<MUT>/           <id>_clean.pdb symlink, chains_to_move.txt, nataa_mutations.resfile, input.json
    mutations.txt           every mutation set up; mutations_<group>.txt per group (scan only);
                            mutations_custom.txt for explicit `mutations` from target.json
Files are only rewritten when their content changes, so rerunning setup is safe while jobs run.
'''
import json, os

from . import config, interface, pdbprep

THREE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
         'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
         'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def clean_name(pdb):
    return os.path.basename(pdb).split('_')[0].split('.')[0] + '_clean.pdb'


def write_if_changed(path, text):
    if os.path.exists(path) and open(path).read() == text:
        return
    tmp = path + '.tmp'
    open(tmp, 'w').write(text)
    os.replace(tmp, path)


def symlink(target, path):
    if os.path.islink(path) and os.readlink(path) == target:
        return
    tmp = path + '.tmp'
    if os.path.lexists(tmp):
        os.remove(tmp)
    os.symlink(target, tmp)
    os.replace(tmp, path)


def resfile(muts):
    return 'NATAA\nstart\n' + ''.join('%d%s %s PIKAA %s\n' % (n, i.strip(), c, m) for c, _, n, i, m in muts)


def setup(target_dir, log=print):
    cfg = config.load(target_dir)
    os.chdir(cfg['dir'])
    atoms = pdbprep.read_atoms(cfg['pdb'])
    present = {l[21] for l in atoms}
    order = config.chain_order(cfg)
    missing = [c for c in order if c not in present]
    if missing:
        raise config.ConfigError('chains %s are not in %s (it has %s)'
                                 % (''.join(missing), cfg['pdb'], ''.join(sorted(present))))
    clean = clean_name(cfg['pdb'])
    seq = pdbprep.sequence_by_residue(atoms)

    breaks = pdbprep.chain_breaks(atoms)
    write_if_changed('chain_breaks.txt', ''.join('%s %d%s -> %d%s  C-N %.1f A\n' % (
        a[0], a[1], a[2].strip(), b[1], b[2].strip(), d) for a, b, d in breaks))
    if breaks:
        log('chain breaks (kept as bonded by Rosetta): ' + ', '.join(
            '%s%d-%d' % (a[0], a[1], b[1]) for a, b, _ in breaks))

    moving = cfg['groups'][cfg['move']]
    partner = [c for c in order if c not in moving]
    scan = cfg['scan']
    sets = []    # (name, muts, list_name)
    keep = []
    if scan['mode'] == 'interface':
        iface = interface.interface_residues(atoms, moving, partner, scan['cutoff'])
        pca = [k for k in iface if k[3] == 'PCA']
        for k in pca:
            log('PCA %s%d at %.1f A from the partner' % (k[0], k[1], iface[k]))
        disulf = interface.disulfide_cys(atoms)
        keep, skipped = interface.select(iface, scan['skip'], disulf, order)
        keep = [k for k in keep if config.group_of(cfg, k[0]) in scan['side']]
        log('interface residues (<%.1f A heavy atom): %d' % (scan['cutoff'], len(iface)))
        log('skipped (%s): %s' % ('/'.join(scan['skip']), ' '.join('%s%d%s%s' % (c, n, i.strip(), r)
                                                                    for c, n, i, r in skipped)))
        log('to scan (%s side): %d' % ('+'.join(scan['side']), len(keep)))
        for c, n, i, r in keep:
            log('  %s %4d%s %s  %.1f A' % (c, n, i.strip(), r, iface[(c, n, i, r)]))
        for c, n, i, r in keep:
            if r not in THREE:
                raise config.ConfigError('cannot mutate nonstandard residue %s %s%d' % (r, c, n))
            if THREE[r] == scan['to']:
                continue
            muts = [(c, THREE[r], n, i, scan['to'])]
            sets.append((config.mutation_name(muts), muts, 'mutations_%s.txt' % config.group_of(cfg, c)))
    for spec in cfg['mutations']:
        muts = config.parse_mutation(spec)
        for c, wt, n, i, _ in muts:
            have = seq.get((c, n, i))
            if have is None:
                raise config.ConfigError('mutation %s: no residue %s%d%s in the PDB' % (spec, c, n, i.strip()))
            if THREE.get(have) != wt:
                raise config.ConfigError('mutation %s: residue %s%d%s is %s, not %s' % (spec, c, n, i.strip(), have, wt))
        sets.append((config.mutation_name(muts), muts, 'mutations_custom.txt'))

    names = [s[0] for s in sets]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        raise config.ConfigError('mutations set up twice: %s' % ', '.join(sorted(dup)))

    pdbprep.write_clean(atoms, clean + '.new')
    if os.path.exists(clean) and open(clean).read() == open(clean + '.new').read():
        os.remove(clean + '.new')
    else:
        os.replace(clean + '.new', clean)
    write_if_changed('interface_residues.txt', ''.join('%s %d%s %s\n' % (c, n, i.strip(), r) for c, n, i, r in keep))

    lists = {}
    for name, muts, list_name in sets:
        d = os.path.join('inputs', name)
        os.makedirs(d, exist_ok=True)
        symlink(os.path.join('..', '..', clean), os.path.join(d, clean))
        write_if_changed(os.path.join(d, 'chains_to_move.txt'), ''.join(moving) + '\n')
        write_if_changed(os.path.join(d, 'nataa_mutations.resfile'), resfile(muts))
        write_if_changed(os.path.join(d, 'input.json'), json.dumps({
            'pdb': clean, 'chains_to_move': ''.join(moving), 'resfile': 'nataa_mutations.resfile',
            'mutations': ['%s:%s%d%s%s' % (c, wt, n, i.strip(), m) for c, wt, n, i, m in muts]},
            indent=1) + '\n')
        lists.setdefault(list_name, []).append(name)
    write_if_changed('mutations.txt', ''.join(n + '\n' for n in names))
    for list_name, ns in lists.items():
        write_if_changed(list_name, ''.join(n + '\n' for n in ns))

    stale = sorted(set(os.listdir('inputs')) - set(names)) if os.path.isdir('inputs') else []
    if stale:
        log('WARNING: inputs/ has directories not in the current config (left alone): ' + ' '.join(stale))
    log('%d mutations: ' % len(names) + ', '.join('%d in %s' % (len(v), k) for k, v in lists.items()))
    return cfg, names
