'''Build design:receptor complexes and score the interface.

The designs Rosetta writes are peptide-only -- chain A, no receptor, because the scaffolds are
excised motifs and chain B is never in the pose. So there is nothing to score an interface on until
one is built: this applies `validate.superimpose_on_motif`'s transform (the same one `validate` uses
to measure motif RMSD and receptor clash) and concatenates the receptor chain, giving a two-chain
PDB with the design as chain P.

Scoring then reuses the route rule from rescore.py -- PeptideCyclizeMover for amide, nothing for the
disulfide routes -- so the peptide is scored as the macrocycle it is rather than as a linear chain
with charged termini sitting in the binding site.

What the numbers are, and are not:

  dG_separated  the complex's score minus the two chains pulled apart and repacked separately. A
                binding score in REU with no reference-state correction, so it is ONLY a relative
                ranking across structures prepared identically. It is not a ddG and will not
                reproduce the flex ddG scans in ../flexddg/.
  ddg           the Ddg filter, which repacks both states. repeats=1 is single-shot and noisy.

The repack is the part that decides whether any of this means anything. On a raw superposed complex
both numbers are dominated by fa_rep, so `pack_input` is on: without it you rank clash severity
rather than binding. It is a SIDECHAIN repack only -- backbone clash from the superposition is not
relieved -- which is the intended behaviour here, since `validate` has already rejected anything
closer than 2.2 A. The one rule that matters is that every design gets the identical treatment.
'''
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from cyclize import config, struct, validate            # noqa: E402
from rescore import CYCLIC_ROUTES, ROSETTA, run         # noqa: E402

PEP_CHAIN = 'P'


def transformed_lines(path, transform, chain=PEP_CHAIN):
    '''The design's ATOM/HETATM lines moved into the complex frame and relabelled.'''
    R, pc, qc = transform
    import numpy as np
    out = []
    for line in open(path):
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        xyz = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])]])
        x, y, z = ((xyz - pc).dot(R.T) + qc)[0]
        out.append('%s%s%s%8.3f%8.3f%8.3f%s' % (line[:21], chain, line[22:30], x, y, z,
                                                line[54:].rstrip('\n')))
    return out


def build_complex(design_path, native, receptor_lines, offset, dest):
    '''-> motif RMSD. Writes a two-chain PDB: design as chain P, receptor as it was.'''
    design = struct.chain_residues(struct.read(design_path), 'A')
    transform, rmsd, _ = validate.superimpose_on_motif(design, native, offset=offset)
    with open(dest, 'w') as fh:
        fh.write('\n'.join(transformed_lines(design_path, transform)) + '\nTER\n')
        fh.write(''.join(receptor_lines))
        fh.write('TER\nEND\n')
    return rmsd


def xml(cyclize_it, receptor_chain):
    mover = ('    <PeptideCyclizeMover name="cyc" residue_selector="pep"/>\n' if cyclize_it else '')
    step = '    <Add mover="cyc"/>\n' if cyclize_it else ''
    return ('<ROSETTASCRIPTS>\n'
            '  <SCOREFXNS><ScoreFunction name="sfxn" weights="ref2015"/></SCOREFXNS>\n'
            '  <RESIDUE_SELECTORS><Chain name="pep" chains="%s"/></RESIDUE_SELECTORS>\n'
            '  <MOVERS>\n' % PEP_CHAIN + mover +
            '    <InterfaceAnalyzerMover name="ia" scorefxn="sfxn" interface="%s_%s"\n'
            '                            pack_separated="true" pack_input="true"\n'
            '                            packstat="false" interface_sc="true"\n'
            '                            scorefile_reporting_prefix="if"/>\n' % (PEP_CHAIN, receptor_chain) +
            '  </MOVERS>\n'
            '  <FILTERS>\n'
            '    <Ddg name="ddg" scorefxn="sfxn" jump="1" repeats="1" repack="1" confidence="0"/>\n'
            '  </FILTERS>\n'
            '  <PROTOCOLS>\n' + step +
            '    <Add mover="ia"/>\n'
            '    <Add filter="ddg"/>\n'
            '  </PROTOCOLS>\n'
            '  <OUTPUT scorefxn="sfxn"/>\n'
            '</ROSETTASCRIPTS>\n')


def read_scorefile(path):
    out, hdr = {}, None
    for line in open(path):
        f = line.split()
        if not f or f[0] != 'SCORE:':
            continue
        if hdr is None:
            hdr = f[1:]
            continue
        row = dict(zip(hdr, f[1:]))
        out[row['description']] = row
    return out


def score_target(target, use_all=False, keep=False):
    cfg = config.load(os.path.join(ROOT, target))
    if not cfg['receptor_chain']:
        sys.exit('%s: no receptor_chain in cyclize.json; there is no interface to score' % target)
    work = os.path.join(cfg['dir'], 'interface')
    cdir = os.path.join(work, 'complexes')
    os.makedirs(cdir, exist_ok=True)

    receptor_lines = [l for l in open(cfg['complex'])
                      if l.startswith(('ATOM', 'HETATM')) and l[21] == cfg['receptor_chain']]
    if not receptor_lines:
        sys.exit('%s: no chain %s in %s' % (target, cfg['receptor_chain'], cfg['complex']))

    sys.path.insert(0, HERE)
    import rescore
    jobs = rescore.collect(target, use_all)
    natives, manifests, built = {}, {}, []
    for path, variant, route in jobs:
        if variant not in manifests:
            manifests[variant] = json.load(
                open(os.path.join(cfg['dir'], 'inputs', variant, 'input.json')))
        man = manifests[variant]
        if man['scaffold'] not in natives:
            natives[man['scaffold']] = struct.chain_residues(
                struct.read(cfg['scaffolds'][man['scaffold']], cfg['motif_chain']), cfg['motif_chain'])
        tag = os.path.basename(path)[:-4]
        dest = os.path.join(cdir, tag + '.pdb')
        # ssapp prepends a Cys, so the motif starts at design position 2 -- the same offset validate
        # uses. Getting it wrong misaligns every pair and the superposition is silently garbage.
        offset = man.get('motif_range', [1])[0] - 1
        rmsd = build_complex(path, natives[man['scaffold']], receptor_lines, offset, dest)
        built.append({'tag': tag, 'pdb': dest, 'variant': variant, 'route': route,
                      'motif_rmsd': round(rmsd, 3)})
    print('%s: built %d complexes -> %s' % (target, len(built), cdir))

    rows = []
    for cyclic, label in ((True, 'cyclic'), (False, 'linear')):
        group = [b for b in built if (b['route'] in CYCLIC_ROUTES) == cyclic]
        if not group:
            continue
        open(os.path.join(work, '%s.xml' % label), 'w').write(xml(cyclic, cfg['receptor_chain']))
        listfile = os.path.join(work, '%s.list' % label)
        open(listfile, 'w').write(''.join(b['pdb'] + '\n' for b in group))
        print('  %s: %d structures (%s)' % (label, len(group),
                                            'PeptideCyclizeMover' if cyclic else 'read as-is'))
        run([ROSETTA, '-l', listfile, '-parser:protocol', '%s.xml' % label,
             '-out:file:score_only', '%s.sc' % label, '-no_nstruct_label',
             '-overwrite', '-mute', 'all', '-out:nooutput'], cwd=work)
        scores = read_scorefile(os.path.join(work, '%s.sc' % label))
        for b in group:
            s = scores.get(b['tag'])
            if s is None:
                print('    WARNING: no score row for %s' % b['tag'])
                continue
            row = {'target': target, 'structure': b['tag'], 'variant': b['variant'],
                   'route': b['route'], 'scaffold': b['variant'].split('__')[0],
                   'turn': b['variant'].split('__')[2], 'motif_rmsd': b['motif_rmsd'],
                   'cyclized': int(cyclic)}
            for key, col in (('dG_separated', 'if_dG_separated'), ('dSASA_int', 'if_dSASA_int'),
                             ('dG_dSASA_norm', 'if_dG_separated/dSASAx100'),
                             ('hbonds_int', 'if_hbonds_int'), ('nres_int', 'if_nres_int'),
                             ('sc_value', 'if_sc_value'), ('ddg', 'ddg'),
                             ('total_score', 'total_score')):
                row[key] = round(float(s[col]), 3) if col in s else ''
            rows.append(row)

    if not keep:
        for f in os.listdir(work):
            if f.endswith('.list'):
                os.remove(os.path.join(work, f))
    out = os.path.join(work, 'interface.csv')
    rows.sort(key=lambda r: (r['scaffold'], r['dG_separated'] if r['dG_separated'] != '' else 1e9))
    with open(out, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print('  -> %s (%d rows)' % (out, len(rows)))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('targets', nargs='*', default=['apj', 'uor'])
    ap.add_argument('--all', action='store_true', help='every passing structure, not just reps')
    ap.add_argument('--keep', action='store_true', help='keep the -l lists')
    args = ap.parse_args(argv)
    if not shutil.which(ROSETTA):
        sys.exit('%s not on PATH; ml biology rosetta/3.14 first' % ROSETTA)
    for target in args.targets:
        rows = score_target(target, args.all, args.keep)
        for r in sorted([r for r in rows if r['dG_separated'] != ''],
                        key=lambda r: r['dG_separated'])[:5]:
            print('      %-34s dG_sep %9.2f  ddg %9.2f  dSASA %7.1f  motif %.2f'
                  % (r['variant'], r['dG_separated'], r['ddg'], r['dSASA_int'], r['motif_rmsd']))


if __name__ == '__main__':
    main()
