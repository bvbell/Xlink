'''Validate every closed design and write analysis_output/.

results.csv     one row per closed structure: motif RMSD, seed phi, closest receptor approach, verdict
summary.csv     one row per variant: closure yield, pass rate, best motif RMSD
'''
import csv, glob, json, os

from . import config, runner, struct, validate


def analyse(cfg):
    scan = config.read_scan(cfg['scan'])
    comp = struct.read(cfg['complex'])
    rec = (struct.chain_residues(comp, cfg['receptor_chain']) if cfg['receptor_chain'] else [])
    if not rec:
        print('no receptor_chain set: clash filtering is skipped, only motif RMSD and seed phi apply')
    outdir = os.path.join(cfg['dir'], 'analysis_output')
    os.makedirs(outdir, exist_ok=True)
    rows, per_variant = [], []
    for vid in [l.strip() for l in open(os.path.join(cfg['dir'], 'variants.txt')) if l.strip()]:
        idir = os.path.join(cfg['dir'], 'inputs', vid)
        man = json.load(open(os.path.join(idir, 'input.json')))
        native = struct.chain_residues(struct.read(cfg['scaffolds'][man['scaffold']],
                                                  cfg['motif_chain']), cfg['motif_chain'])
        seed_num = man['n_motif'] + 1 if man['route'] == 'amide' and man['added'] else None
        # ssapp prepends a Cys, so the motif starts at design position 2, not 1.
        offset = man.get('motif_range', [1])[0] - 1
        motif_nums = range(offset + 1, offset + man['n_motif'] + 1)
        got = sorted(glob.glob(os.path.join(runner.variant_dir(cfg, vid), '*.pdb')))
        passes, best = 0, None
        for f in got:
            try:
                r = validate.check(f, native, rec, motif_nums, cfg['validate']['clash'],
                                   cfg['validate']['motif_rmsd'], seed_num, offset=offset)
            except (struct.StructError, ValueError) as e:
                rows.append({'variant': vid, 'structure': os.path.basename(f), 'pass': False,
                             'reasons': 'unreadable: %s' % e})
                continue
            r['variant'], r['structure'] = vid, os.path.basename(f)
            r['route'], r['scaffold'] = man['route'], man['scaffold']
            rows.append(r)
            passes += r['pass']
            if best is None or r['motif_rmsd'] < best:
                best = r['motif_rmsd']
        per_variant.append({'variant': vid, 'scaffold': man['scaffold'], 'route': man['route'],
                            'added': '-'.join(man['added']), 'turn': man['turn'] or '',
                            'nstruct': cfg['protocol']['nstruct'], 'closed': len(got),
                            'closure_yield': round(len(got) / cfg['protocol']['nstruct'], 3),
                            'passed': passes,
                            'pass_rate_of_closed': round(passes / len(got), 3) if got else 0.0,
                            'best_motif_rmsd': round(best, 2) if best is not None else ''})
    _write(os.path.join(outdir, 'results.csv'), rows,
           ['variant', 'scaffold', 'route', 'structure', 'motif_rmsd', 'seed_phi',
            'min_receptor_dist', 'closest_design_res', 'closest_receptor_res', 'pass', 'reasons'])
    per_variant.sort(key=lambda r: (-r['passed'], r['variant']))
    _write(os.path.join(outdir, 'summary.csv'), per_variant,
           ['variant', 'scaffold', 'route', 'added', 'turn', 'nstruct', 'closed', 'closure_yield',
            'passed', 'pass_rate_of_closed', 'best_motif_rmsd'])
    print('%-46s %6s %7s %7s %9s' % ('variant', 'closed', 'passed', 'yield', 'bestRMSD'))
    for r in per_variant:
        print('%-46s %6d %7d %7.2f %9s' % (r['variant'], r['closed'], r['passed'],
                                           r['closure_yield'], r['best_motif_rmsd']))
    print('\n%d structures over %d variants; %d passed'
          % (len(rows), len(per_variant), sum(r['passed'] for r in per_variant)))
    print('wrote %s/{results,summary}.csv' % outdir)
    return per_variant


def _write(path, rows, fields):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)
