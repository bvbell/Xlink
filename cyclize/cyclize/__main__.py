'''cyclize <command> <target_dir> ...   (run as `python3 -m cyclize` or bin/cyclize)

    register   hairpin register report for the scaffolds, without writing anything
    setup      enumerate variants, write inputs/ and variants.txt from cyclize.json
    status     per-variant closure state, plus active jobs
    submit     queue only variants not yet run, then a validation job after them
    validate   analysis_output/{results,summary}.csv from the closed structures
    smoke      one variant at low nstruct in a scratch dir, as a check before a full submit
    run-task   (used by the array job) run line SLURM_ARRAY_TASK_ID of a task file
'''
import argparse, os, sys

from . import config


def main(argv=None):
    p = argparse.ArgumentParser(prog='cyclize', description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    for name in ['register', 'setup', 'status', 'validate']:
        sub.add_parser(name).add_argument('target')
    s = sub.add_parser('submit')
    s.add_argument('target')
    s.add_argument('--retry', action='store_true', help='also rerun partial/empty variants (find out why first)')
    s.add_argument('--no-analysis', action='store_true')
    s.add_argument('--dry-run', action='store_true')
    s.add_argument('--force', action='store_true', help='submit even though this target has active jobs')
    s = sub.add_parser('smoke')
    s.add_argument('target')
    s.add_argument('variant', nargs='?', help='default: first in variants.txt')
    s.add_argument('--nstruct', type=int, default=4)
    s = sub.add_parser('run-task')
    s.add_argument('target')
    s.add_argument('taskfile')
    s.add_argument('index', type=int)
    a = p.parse_args(argv)

    try:
        if a.cmd == 'setup':
            from . import inputs
            inputs.setup(a.target)
            return 0
        cfg = config.load(a.target)
    except config.ConfigError as e:
        sys.exit('error: %s' % e)

    if a.cmd == 'register':
        from . import hairpin, inputs
        scan = config.read_scan(cfg['scan'])
        for name in sorted(cfg['scaffolds']):
            print('=== %s ===' % name)
            print(hairpin.describe(inputs.load_scaffold(cfg, name), scan=scan))
    elif a.cmd == 'status':
        from . import slurm
        slurm.status(cfg)
    elif a.cmd == 'submit':
        from . import slurm
        slurm.submit(cfg, retry=a.retry, then_analyze=not a.no_analysis, dry_run=a.dry_run,
                     force=a.force)
    elif a.cmd == 'validate':
        from . import analysis
        analysis.analyse(cfg)
    elif a.cmd == 'smoke':
        from . import runner, slurm
        vid = a.variant or slurm.selected_variants(cfg)[0]
        base = os.path.join(os.environ.get('SCRATCH', '/tmp'), 'cyclize_smoke',
                            os.path.basename(cfg['dir']))
        os.makedirs(base, exist_ok=True)
        cfg = dict(cfg, protocol=dict(cfg['protocol'], nstruct=a.nstruct))
        print('smoke: %s, %d structures -> %s' % (vid, a.nstruct, base))
        return runner.run_variant(cfg, vid, base=base)
    elif a.cmd == 'run-task':
        from . import slurm
        return slurm.run_task(cfg, a.taskfile, a.index)
    return 0


if __name__ == '__main__':
    sys.exit(main())
