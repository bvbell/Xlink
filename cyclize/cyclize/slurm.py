'''Submit one array task per variant, report progress, chain the analysis after the runs.

Follows the same discipline as the flexddg package: each submission writes tasks/<timestamp>.txt (one
variant id per line) and appends "<jobid> <taskfile> <ntasks>" to tasks/submissions.tsv, so a log can
always be traced back to the variant it ran.
'''
import datetime, os, shlex, subprocess, sys

from . import runner

PYTHON = sys.executable
PKG_PARENT = os.path.dirname(runner.PKG)


def job_name(cfg):
    return 'cyclize_' + os.path.basename(cfg['dir'])


def selected_variants(cfg):
    path = os.path.join(cfg['dir'], 'variants.txt')
    if not os.path.exists(path):
        sys.exit('%s not found; run `cyclize setup` first' % path)
    return [l.strip() for l in open(path) if l.strip()]


def active_jobs(cfg):
    try:
        out = subprocess.run(['squeue', '--me', '-h', '-r', '-n', job_name(cfg), '-o', '%F %T'],
                             capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError) as e:
        sys.exit('squeue failed (%s); not guessing at job state' % e)
    jobs = {}
    for line in out.splitlines():
        jid, st = line.split()
        jobs.setdefault(jid, {}).setdefault(st, 0)
        jobs[jid][st] += 1
    return jobs


def scan_states(cfg, vids=None):
    return {v: runner.state(cfg, v) for v in (vids or selected_variants(cfg))}


def status(cfg):
    jobs = active_jobs(cfg)
    states = scan_states(cfg)
    total = {}
    print('%-46s %-8s %s' % ('variant', 'state', 'closed'))
    for v, st in states.items():
        print('%-46s %-8s %d' % (v, st, runner.closed_count(runner.variant_dir(cfg, v))))
        total[st] = total.get(st, 0) + 1
    print('total:', ', '.join('%s %d' % kv for kv in sorted(total.items())))
    for jid, st in jobs.items():
        print('active job %s: %s' % (jid, ', '.join('%d %s' % (n, s.lower()) for s, n in st.items())))
    if not jobs and total.get('empty'):
        print('"empty" with no active jobs = ran but nothing closed. Check inputs/<variant>/protocol.xml '
              'and the loop size reported in setup.log.')
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


def submit(cfg, retry=False, then_analyze=True, dry_run=False, force=False):
    jobs = active_jobs(cfg)
    if jobs and not force:
        sys.exit('target already has active jobs (%s); wait for them, or --force if they cover other '
                 'variants. Pending is not a reason to resubmit.' % ', '.join(jobs))
    states = scan_states(cfg)
    want = ['missing'] + (['partial', 'empty'] if retry else [])
    tasks = [v for v, st in states.items() if st in want]
    stale = sum(1 for st in states.values() if st in ('partial', 'empty'))
    if stale and not retry:
        print('%d variants are partial/empty and were not resubmitted; check `cyclize status` first, '
              'then pass --retry.' % stale)
    if not tasks:
        print('nothing to submit')
        return None
    print('%d variants to run (%d structures each)' % (len(tasks), cfg['protocol']['nstruct']))
    if dry_run:
        for t in tasks:
            print('  ' + t)
        return None
    os.makedirs(os.path.join(cfg['dir'], 'logs'), exist_ok=True)
    taskfile = os.path.join(cfg['dir'], 'tasks',
                            datetime.datetime.now().strftime('%Y%m%d-%H%M%S') + '.txt')
    os.makedirs(os.path.dirname(taskfile), exist_ok=True)
    s = cfg['slurm']
    # One variant takes about a minute, and a Slurm array of one-minute tasks is poor scheduler
    # etiquette (and mostly queue latency). Pack several variants onto each array task; one line of
    # the task file is one task. A line with a single id is still valid, so old task files replay.
    per = max(1, int(s['per_task']))
    lines = [tasks[i:i + per] for i in range(0, len(tasks), per)]
    open(taskfile, 'w').write(''.join(' '.join(chunk) + '\n' for chunk in lines))
    script = write_script(taskfile[:-4] + '.sbatch', [
        'ml purge; ml %s' % cfg['rosetta_module'],
        'echo "cyclize task $SLURM_ARRAY_TASK_ID on $(hostname)"',
        'exec %s -m cyclize run-task %s %s "$SLURM_ARRAY_TASK_ID"' % (
            PYTHON, shlex.quote(cfg['dir']), shlex.quote(taskfile))])
    jid = sbatch(['-D', cfg['dir'], '-J', job_name(cfg), '-p', s['partition'], '-t', s['time'],
                  '-c', '1', '--mem', s['mem'],
                  '--array=0-%d%%%d' % (len(lines) - 1, s['max_concurrent']),
                  '-o', os.path.join(cfg['dir'], 'logs', '%A_%a.out'), '--export', env_export(), script])
    with open(os.path.join(cfg['dir'], 'tasks', 'submissions.tsv'), 'a') as f:
        f.write('%s\t%s\t%d\n' % (jid, os.path.basename(taskfile), len(tasks)))
    print('submitted array job %s (%d variants in %d tasks of up to %d, %d tasks at a time)'
          % (jid, len(tasks), len(lines), per, s['max_concurrent']))
    if then_analyze:
        submit_analysis(cfg, after=jid)
    return jid


def submit_analysis(cfg, after=None):
    s = cfg['slurm']
    os.makedirs(os.path.join(cfg['dir'], 'logs'), exist_ok=True)
    script = write_script(os.path.join(cfg['dir'], 'tasks', 'analysis.sbatch'),
                          ['%s -m cyclize validate %s' % (PYTHON, shlex.quote(cfg['dir']))])
    args = ['-D', cfg['dir'], '-J', job_name(cfg) + '_analysis', '-p', s['partition'],
            '-t', s['analysis_time'], '-c', '1', '--mem', s['analysis_mem'],
            '-o', os.path.join(cfg['dir'], 'logs', 'analysis_%j.out'), '--export', env_export(), script]
    if after:
        args = ['--dependency=afterany:%s' % after] + args
    jid = sbatch(args)
    print('submitted analysis job %s%s' % (jid, ' (after %s)' % after if after else ''))
    return jid


def run_task(cfg, taskfile, index):
    '''Run every variant on line `index` of the task file. One line = one array task.'''
    lines = [l.split() for l in open(taskfile) if l.strip()]
    vids = lines[int(index)]
    # Seed from the variant's position in the whole file, not its position in the line, so every
    # variant keeps a distinct RNG stream however the file happens to be chunked.
    first = sum(len(l) for l in lines[:int(index)])
    rc = 0
    for i, vid in enumerate(vids):
        print('[%d/%d] %s; target %s' % (i + 1, len(vids), vid, cfg['dir']), flush=True)
        try:
            rc |= runner.run_variant(cfg, vid, seed=first + i + 1)
        except Exception as e:       # one bad variant must not take the other 19 down with it
            print('%s FAILED: %s' % (vid, e), flush=True)
            rc = 1
    return rc
