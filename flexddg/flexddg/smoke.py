'''One short trajectory in $SCRATCH/flexddg_smoke/, never in the target's output/, so a test run can
never be mistaken for a finished production trajectory. Keeps protocol.maxmin: a small maxmin makes
the minimizer restart until it converges, which is slower, not faster.
'''
import datetime, os, shlex, sys

from . import analysis, runner, slurm


def smoke(cfg, mutation=None, ntrials=100, use_sbatch=False):
    mutation = mutation or slurm.selected_mutations(cfg)[0]
    if not os.path.isdir(os.path.join(cfg['dir'], 'inputs', mutation)):
        sys.exit('no inputs/%s; run setup first' % mutation)
    base = os.path.join(os.environ.get('SCRATCH', '/tmp'), 'flexddg_smoke',
                        '%s-%s' % (os.path.basename(cfg['dir']), datetime.datetime.now().strftime('%Y%m%d-%H%M%S')))
    if not use_sbatch and not os.environ.get('SLURM_JOB_ID'):
        sys.exit('not inside a Slurm allocation; rerun with --sbatch (or inside sh_dev)')
    os.makedirs(base)
    if use_sbatch:
        script = slurm.write_script(os.path.join(base, 'smoke.sbatch'), [
            'ml purge; ml %s' % cfg['rosetta_module'],
            'exec %s -m flexddg smoke %s %s --ntrials %d' % (slurm.PYTHON, shlex.quote(cfg['dir']), mutation, ntrials)])
        jid = slurm.sbatch(['-J', slurm.job_name(cfg) + '_smoke', '-p', cfg['slurm']['partition'], '-t', '01:00:00',
                            '-c', '1', '--mem', cfg['slurm']['mem'], '-o', os.path.join(base, 'smoke_%j.out'),
                            '--export', slurm.env_export(), script])
        print('submitted smoke job %s; log in %s' % (jid, base))
        return 0
    protocol = dict(cfg['protocol'], ntrials=ntrials, stride=ntrials)
    rc = runner.run_trajectory(cfg, mutation, 1, protocol=protocol, base=base)
    test_cfg = dict(cfg, protocol=protocol, analysis=dict(cfg['analysis'], min_nstruct=1, bootstrap=10))
    d = runner.traj_dir(cfg, mutation, 1, base)
    if rc != 0 or runner.traj_state(test_cfg, d) != 'done':
        print('SMOKE FAILED: see %s/rosetta.out' % d)
        return 1
    r = analysis.summarize(test_cfg, out_dir=os.path.join(base, 'analysis_output'), finished={mutation: [d]})
    print('SMOKE OK: %s ddG_GAM %.2f from one %d-trial trajectory (not meaningful, only parsed); %s'
          % (mutation, r.ddG_GAM.iloc[0], ntrials, base))
    return 0
