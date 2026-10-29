'''setup: read the scaffolds, analyse the register, enumerate variants, write inputs/ and variants.txt.

Everything here is derived, so setup is safe to rerun. It refuses outright on an off-register scaffold
rather than quietly producing designs whose H-bond ladder cannot survive closure.
'''
import os

from . import build, config, hairpin, struct, variants


def load_scaffold(cfg, name):
    res = struct.chain_residues(struct.read(cfg['scaffolds'][name], cfg['motif_chain']),
                                cfg['motif_chain'])
    if len(res) < 4:
        raise variants.VariantError('%s: %d residues is below the 4-residue floor for a cyclic peptide'
                                    % (name, len(res)))
    return res


def setup(target_dir):
    cfg = config.load(target_dir)
    scan = config.read_scan(cfg['scan'])
    report, all_variants, problems = [], [], []
    for name in sorted(cfg['scaffolds']):
        res = load_scaffold(cfg, name)
        report.append('=== %s ===' % name)
        report.append(hairpin.describe(res, scan=scan))
        try:
            vs = variants.enumerate_variants(name, res, scan, cfg['hotspot'], cfg['routes'],
                                             target_loop=cfg['protocol'].get('target_loop',
                                                                             variants.TARGET_LOOP))
        except variants.VariantError as e:
            problems.append(str(e))
            report.append('  SKIPPED: %s' % e)
            continue
        if vs and vs[0].get('underdetermined'):
            report.append('  WARNING: only %d loop residues available (needs > %d for GenKIC to have '
                          'any sampling freedom); closure will probably fail for this scaffold'
                          % (vs[0]['loop_size'], variants.MIN_LOOP))
        report.append('  %d variants, loop %d residues, releasing N%s C%s'
                      % (len(vs), vs[0]['loop_size'], vs[0]['release_n'], vs[0]['release_c']))
        for v in vs:
            build.write_variant(cfg, v, res, os.path.join(cfg['dir'], 'inputs', v['id']))
        all_variants.extend(vs)
    open(os.path.join(cfg['dir'], 'variants.txt'), 'w').write(
        ''.join(v['id'] + '\n' for v in all_variants))
    text = '\n'.join(report) + '\n\n' + variants.summary(all_variants) + '\n'
    open(os.path.join(cfg['dir'], 'setup.log'), 'w').write(text)
    print(text)
    print('%d variants written to %s/inputs/' % (len(all_variants), cfg['dir']))
    if problems:
        print('\n%d scaffold(s) skipped:' % len(problems))
        for p in problems:
            print('  ' + p)
    return all_variants
