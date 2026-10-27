'''Submit only the trajectories still to do, report progress, chain the analysis after the runs.

Each submission writes tasks/<timestamp>.txt ("MUT I" per line, one array task each) and appends
"<jobid> <taskfile> <ntasks>" to tasks/submissions.tsv, so a log can be traced back to its task.
'''
import datetime, os, shlex, subprocess, sys

from . import runner

PYTHON = sys.executable
PKG_PARENT = os.path.dirname(runner.PKG)


def job_name(cfg):
    return 'flexddg_' + os.path.basename(cfg['dir'])


def selected_mutations(cfg):
    return [l.strip() for l in open(os.path.join(cfg['dir'], 'mutations.txt')) if l.strip()]


def active_jobs(cfg):
    '''{jobid: "N running, M pending"} for this target's run jobs (not the analysis job).'''
    try:
        out = subprocess.run(['squeue', '--me', '-h', '-r', '-n', job_name(cfg), '-o', '%F %T'],
                             capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError) as e:
        sys.exit('squeue failed (%s); not guessing at job state' % e)
    jobs = {}
    for line in out.splitlines():
        jid, state = line.split()
        jobs.setdefault(jid, {}).setdefault(state, 0)
        jobs[jid][state] += 1
    return jobs


def scan_states(cfg, muts=None):
    '''{mut: {state: [i, ...]}} over nstruct trajectories.'''
    out = {}
    for m in muts or selected_mutations(cfg):
        states = {}
        for i in range(1, cfg['protocol']['nstruct'] + 1):
            states.setdefault(runner.traj_state(cfg, runner.traj_dir(cfg, m, i)), []).append(i)
        out[m] = states
    return out


def status(cfg):
    jobs = active_jobs(cfg)
    states = scan_states(cfg)
    total = {}
    print('%-18s %5s %5s %5s %5s' % ('mutation', 'done', 'incmp', 'miss', 'other'))
    for m, s in states.items():
        print('%-18s %5d %5d %5d %5d' % (m, len(s.get('done', [])), len(s.get('incomplete', [])),
                                        len(s.get('missing', [])), len(s.get('mismatch', []))))
        for k, v in s.items():
            total[k] = total.get(k, 0) + len(v)
    print('total:', ', '.join('%s %d' % kv for kv in sorted(total.items())))
    if jobs:
        for jid, st in jobs.items():
            print('active job %s: %s' % (jid, ', '.join('%d %s' % (n, s.lower()) for s, n in st.items())))
        print('"incomplete" includes trajectories still running.')
    elif total.get('incomplete'):
        print('no active jobs, so incomplete = failed. Last lines of each:')
        for m, s in states.items():
            for i in s.get('incomplete', []):
                out = os.path.join(runner.traj_dir(cfg, m, i), 'rosetta.out')
                tail = open(out, errors='replace').read().strip().splitlines()[-3:]
                print('--- %s %d\n%s' % (m, i, '\n'.join(tail)))
    if total.get('mismatch'):
        print('"other": finished with a different number_backrub_trials than target.json; never resubmitted.')
    return states, jobs


def sbatch(args):
    out = subprocess.run(['sbatch', '--parsable'] + args, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit('sbatch failed: ' + out.stderr.strip())
    return out.stdout.strip().split(';')[0]


def write_script(path, lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, 'w').write('#!/bin/bash\nset -uo pipefail\n' + '\n'.join(lines) + '\n')
    return path


def env_export():
    path = PKG_PARENT + (':' + os.environ['PYTHONPATH'] if os.environ.get('PYTHONPATH') else '')
    return 'ALL,PYTHONPATH=' + path


def submit(cfg, retry_failed=False, then_analyze=True, dry_run=False, force=False):
    jobs = active_jobs(cfg)
    if jobs and not force:
        sys.exit('target already has active jobs (%s); wait for them, or --force if they cover other '
                 'mutations. Pending is not a reason to resubmit.' % ', '.join(jobs))
    states = scan_states(cfg)
    want = ['missing'] + (['incomplete'] if retry_failed else [])
    tasks = [(m, i) for m, s in states.items() for k in want for i in s.get(k, [])]
    failed = sum(len(s.get('incomplete', [])) for s in states.values())
    if failed and not retry_failed:
        print('%d incomplete trajectories not resubmitted: check why with `flexddg status`, then pass '
              '--retry-failed.' % failed)
    if not tasks:
        print('nothing to submit')
        return None
    tasks.sort(key=lambda t: (t[1], t[0]))   # spread mutations across the array, so partial results cover all of them
    print('%d trajectories to run (%d mutations)' % (len(tasks), len({m for m, _ in tasks})))
    if dry_run:
        return None
    os.makedirs(os.path.join(cfg['dir'], 'tasks'), exist_ok=True)
    os.makedirs(os.path.join(cfg['dir'], 'logs'), exist_ok=True)
    taskfile = os.path.join(cfg['dir'], 'tasks', datetime.datetime.now().strftime('%Y%m%d-%H%M%S') + '.txt')
    open(taskfile, 'w').write(''.join('%s %d\n' % t for t in tasks))
    s = cfg['slurm']
    script = write_script(taskfile[:-4] + '.sbatch', [
        'ml purge; ml %s' % cfg['rosetta_module'],
        'echo "flexddg task $SLURM_ARRAY_TASK_ID on $(hostname); python %s"' % PYTHON,
        'exec %s -m flexddg run-task %s %s "$SLURM_ARRAY_TASK_ID"' % (
            PYTHON, shlex.quote(cfg['dir']), shlex.quote(taskfile))])
    jid = sbatch(['-D', cfg['dir'], '-J', job_name(cfg), '-p', s['partition'], '-t', s['time'],
                  '-c', '1', '--mem', s['mem'], '--array=0-%d%%%d' % (len(tasks) - 1, s['max_concurrent']),
                  '-o', os.path.join(cfg['dir'], 'logs', '%A_%a.out'), '--export', env_export(), script])
    with open(os.path.join(cfg['dir'], 'tasks', 'submissions.tsv'), 'a') as f:
        f.write('%s\t%s\t%d\n' % (jid, os.path.basename(taskfile), len(tasks)))
    print('submitted array job %s (%d tasks, %d at a time)' % (jid, len(tasks), s['max_concurrent']))
    if then_analyze:
        submit_analysis(cfg, after=jid)
    return jid


def submit_analysis(cfg, after=None):
    '''summarize + plot as a job, optionally after job `after` ends (afterany: partial results still count).'''
    s = cfg['slurm']
    os.makedirs(os.path.join(cfg['dir'], 'logs'), exist_ok=True)
    script = write_script(os.path.join(cfg['dir'], 'tasks', 'analysis.sbatch'), [
        '%s -m flexddg summarize %s && %s -m flexddg plot %s' % (
            PYTHON, shlex.quote(cfg['dir']), PYTHON, shlex.quote(cfg['dir']))])
    args = ['-D', cfg['dir'], '-J', job_name(cfg) + '_analysis', '-p', s['partition'],
            '-t', s['analysis_time'], '-c', '1', '--mem', s['analysis_mem'],
            '-o', os.path.join(cfg['dir'], 'logs', 'analysis_%j.out'), '--export', env_export(), script]
    if after:
        args = ['--dependency=afterany:%s' % after] + args
    jid = sbatch(args)
    print('submitted analysis job %s%s' % (jid, ' (after %s)' % after if after else ''))
    return jid


def run_task(cfg, taskfile, index):
    lines = [l.split() for l in open(taskfile) if l.strip()]
    mut, i = lines[int(index)]
    print('%s nstruct %s; target %s' % (mut, i, cfg['dir']), flush=True)
    return runner.run_trajectory(cfg, mut, int(i))
