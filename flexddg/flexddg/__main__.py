'''flexddg <command> <target_dir> ...   (run as `python3 -m flexddg` or bin/flexddg)

    setup      clean PDB, interface scan, inputs/ and mutations*.txt from target.json
    status     done / incomplete / missing trajectories per mutation, plus active jobs
    submit     queue only unfinished trajectories, then an analysis job after them
    analyze    queue the analysis job alone, optionally --after JOBID
    summarize  analysis_output/{results,struct_scores,alascan_summary}.csv
    plot       plots/ from analysis_output/
    compare    side-by-side hotspot figure for several targets
    smoke      one short trajectory in a scratch dir, as a check before a full submit
    run-task   (used by the array job) run line SLURM_ARRAY_TASK_ID of a task file
'''
import argparse, os, sys

from . import config


def main(argv=None):
    p = argparse.ArgumentParser(prog='flexddg', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    for name in ['setup', 'status', 'plot']:
        sub.add_parser(name).add_argument('target')
    s = sub.add_parser('submit')
    s.add_argument('target')
    s.add_argument('--retry-failed', action='store_true', help='also rerun incomplete trajectories (find out why first)')
    s.add_argument('--no-analysis', action='store_true', help='do not queue summarize+plot after the array')
    s.add_argument('--dry-run', action='store_true')
    s.add_argument('--force', action='store_true', help='submit even though this target has active jobs')
    s = sub.add_parser('analyze')
    s.add_argument('target')
    s.add_argument('--after', help='job id to wait for (afterany)')
    s = sub.add_parser('summarize')
    s.add_argument('target')
    s.add_argument('--min-nstruct', type=int, help='default: analysis.min_nstruct, else protocol.nstruct')
    s.add_argument('--all-output', action='store_true', help='every case in output/, not just mutations.txt')
    s = sub.add_parser('compare')
    s.add_argument('targets', nargs='+')
    s.add_argument('-o', '--out', default='compare')
    s = sub.add_parser('smoke')
    s.add_argument('target')
    s.add_argument('mutation', nargs='?', help='default: first in mutations.txt')
    s.add_argument('--ntrials', type=int, default=100)
    s.add_argument('--sbatch', action='store_true', help='submit as a job instead of running here')
    s = sub.add_parser('run-task')
    s.add_argument('target')
    s.add_argument('taskfile')
    s.add_argument('index', type=int)
    a = p.parse_args(argv)

    try:
        if a.cmd == 'compare':
            from . import plots
            plots.compare([config.load(t) for t in a.targets], os.path.abspath(a.out))
            return 0
        if a.cmd == 'setup':
            from . import inputs
            inputs.setup(a.target)
            return 0
        cfg = config.load(a.target)
    except config.ConfigError as e:
        sys.exit('error: %s' % e)

    if a.cmd == 'status':
        from . import slurm
        slurm.status(cfg)
    elif a.cmd == 'submit':
        from . import slurm
        slurm.submit(cfg, retry_failed=a.retry_failed, then_analyze=not a.no_analysis, dry_run=a.dry_run, force=a.force)
    elif a.cmd == 'analyze':
        from . import slurm
        slurm.submit_analysis(cfg, after=a.after)
    elif a.cmd == 'summarize':
        from . import analysis
        if a.min_nstruct:
            cfg['analysis']['min_nstruct'] = a.min_nstruct
        analysis.summarize(cfg, all_output=a.all_output)
    elif a.cmd == 'plot':
        from . import plots
        plots.plot_target(cfg)
    elif a.cmd == 'smoke':
        from . import smoke
        return smoke.smoke(cfg, a.mutation, a.ntrials, a.sbatch)
    elif a.cmd == 'run-task':
        from . import slurm
        return slurm.run_task(cfg, a.taskfile, a.index)
    return 0


if __name__ == '__main__':
    sys.exit(main())
