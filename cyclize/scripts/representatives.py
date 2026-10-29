'''Reduce each condition's passing structures to its distinct conformations.

Every variant ran 50 GenKIC trajectories on one fixed sequence and most of them landed in the same
basin (see scripts/degeneracy.py), so the passing set is far smaller than it looks. This picks one
representative per conformational cluster per variant: cluster the passers by backbone RMSD, then
take the member with the lowest motif RMSD, since that is the only structure-quality number the
pipeline produces. The representatives are what any downstream rescoring should run on -- scoring all
50 members of a cluster is 50x the cost for one answer.

Writes <target>/representatives/<variant>__c<NN>.pdb plus a manifest CSV.
'''
import argparse
import csv
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from degeneracy import backbone, cluster, pairwise_rmsd   # noqa: E402

ROOT = os.path.dirname(HERE)


def representatives(target, cutoff):
    '''-> list of manifest rows, one per conformational cluster.'''
    root = os.path.join(ROOT, target)
    rows = [r for r in csv.DictReader(open(os.path.join(root, 'analysis_output/results.csv')))
            if r['pass'] == 'True']
    by_variant = {}
    for r in rows:
        by_variant.setdefault(r['variant'], {})[r['structure']] = r

    out = []
    for variant in sorted(by_variant):
        members = by_variant[variant]
        names = sorted(members)
        coords = [backbone(os.path.join(root, 'output', variant, n))[1] for n in names]
        d = pairwise_rmsd(coords)
        for rank, (_, idx) in enumerate(cluster(d, cutoff)):
            best = min(idx, key=lambda i: float(members[names[i]]['motif_rmsd']))
            row = members[names[best]]
            spread = max((d[best, i] for i in idx), default=0.0)
            out.append({
                'target': target, 'variant': variant, 'cluster': rank, 'size': len(idx),
                'fraction_of_passing': round(len(idx) / float(len(names)), 3),
                'representative': names[best], 'max_rmsd_to_rep': round(float(spread), 3),
                'motif_rmsd': round(float(row['motif_rmsd']), 3),
                'min_receptor_dist': row['min_receptor_dist'],
                'scaffold': row['scaffold'], 'route': row['route'],
                'turn': variant.split('__')[2],
                'pdb': '%s__c%02d.pdb' % (variant, rank),
            })
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('targets', nargs='*', default=['apj', 'uor'])
    ap.add_argument('--cutoff', type=float, default=1.0,
                    help='backbone RMSD (A) within which two structures are the same conformation')
    ap.add_argument('--dry-run', action='store_true', help='report only, copy nothing')
    args = ap.parse_args(argv)

    for target in args.targets:
        rows = representatives(target, args.cutoff)
        dest = os.path.join(ROOT, target, 'representatives')
        if not args.dry_run:
            if os.path.isdir(dest):
                shutil.rmtree(dest)       # stale reps from a different cutoff would be misleading
            os.makedirs(dest)
            for r in rows:
                shutil.copyfile(os.path.join(ROOT, target, 'output', r['variant'], r['representative']),
                                os.path.join(dest, r['pdb']))
            with open(os.path.join(dest, 'representatives.csv'), 'w') as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)

        passing = sum(r['size'] for r in rows)
        singles = sum(1 for r in rows if r['size'] == 1)
        print('%-5s %4d passing structures -> %3d distinct conformations at %.2f A  '
              '(%d are singletons)%s'
              % (target, passing, len(rows), args.cutoff, singles,
                 '' if args.dry_run else '  -> %s' % dest))
        for r in sorted(rows, key=lambda r: -r['size'])[:5]:
            print('      %-34s c%02d  %3d/%-3d of passing  motif %.2f A'
                  % (r['variant'], r['cluster'], r['size'],
                     int(round(r['size'] / r['fraction_of_passing'])), r['motif_rmsd']))


if __name__ == '__main__':
    main()
