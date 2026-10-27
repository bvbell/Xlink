'''Analyze a target's trajectories with the tutorial's analyze_flex_ddG.py and summarize per mutation.

Writes into <target>/analysis_output/ (or out_dir):
    results.csv          ddG, wt_dG, mut_dG per mutation at nstruct checkpoints 10, 20, 30, ..., n
                         for talaris2014 and its GAM reweighting (analyze_flex_ddG.py format)
    struct_scores.csv    per-trajectory ddG and per-term breakdown (talaris)
    alascan_summary.csv  per mutation: n, talaris mean and SEM, GAM ddG with a bootstrap 95% CI over
                         trajectories, GAM at 20 and 30 nstruct (convergence), fraction of ddG > 0
Only mutations with at least analysis.min_nstruct finished trajectories (default: protocol.nstruct)
are included, and only trajectories run with the configured number_backrub_trials.
'''
import os, sys

import numpy as np
import pandas as pd

from . import runner, slurm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vendor'))
import analyze_flex_ddG as vendor  # noqa: E402


def finished_trajectories(cfg, muts=None):
    '''{mut: [trajectory dirs]} for trajectories in state "done".'''
    out = {}
    for m in muts or sorted(os.listdir(os.path.join(cfg['dir'], 'output'))):
        dirs = [runner.traj_dir(cfg, m, i) for i in range(1, cfg['protocol']['nstruct'] + 1)]
        out[m] = [d for d in dirs if runner.traj_state(cfg, d) == 'done']
    return out


def analyze(finished, out_dir, min_nstruct, stride, log=print):
    '''Run the tutorial analysis over {case: [struct dirs]}; returns (results, struct_scores).

    stride labels the checkpoints (backrub_steps). Rosetta 3.14 does not record it in ddG.db3, so it
    comes from target.json rather than the vendor fallback of 5.
    '''
    ddg_dfs, struct_dfs, skipped = [], [], []
    for case, structs in sorted(finished.items()):
        if len(structs) < min_nstruct:
            skipped.append('%s (%d)' % (case, len(structs)))
            continue
        scores = pd.concat([vendor.process_finished_struct(s, case, stride)[0] for s in structs])
        ddg, struct_scores = vendor.calc_ddg(scores)
        struct_dfs.append(struct_scores)
        ddg_dfs += [ddg, vendor.apply_zemu_gam(ddg)] + vendor.calc_dgs(scores)
    if skipped:
        log('skipped, fewer than %d finished trajectories: %s' % (min_nstruct, ', '.join(skipped)))
    if not struct_dfs:
        raise SystemExit('no mutation has %d finished trajectories yet' % min_nstruct)
    os.makedirs(out_dir, exist_ok=True)
    results, struct_scores = pd.concat(ddg_dfs), pd.concat(struct_dfs)
    results.to_csv(os.path.join(out_dir, 'results.csv'))
    struct_scores.to_csv(os.path.join(out_dir, 'struct_scores.csv'))
    return (pd.read_csv(os.path.join(out_dir, 'results.csv')),
            pd.read_csv(os.path.join(out_dir, 'struct_scores.csv')))


def gam_of_mean(df, terms):
    m = df[terms].mean()
    return sum(vendor.gam_function(m[t], t) for t in terms)


def summary(results, struct_scores, bootstrap=1000, seed=0):
    terms = list(vendor.zemu_gam_params)
    dd = results[(results.scored_state == 'ddG') & (results.score_function_name == 'fa_talaris2014-gam')]
    gam = dd.set_index(['case_name', 'nstruct']).total_score
    rng = np.random.default_rng(seed)
    rows = []
    for c, g in struct_scores.groupby('case_name'):
        tal = g.total_score.values
        boots = [gam_of_mean(g.iloc[rng.integers(0, len(g), len(g))], terms) for _ in range(bootstrap)]
        lo, hi = np.percentile(boots, [2.5, 97.5])
        rows.append(dict(case=c, n=len(g), ddG_tal=tal.mean(), tal_sem=tal.std(ddof=1) / np.sqrt(len(tal)),
                         ddG_GAM=gam[(c, g.struct_num.max())], GAM_lo=lo, GAM_hi=hi, GAM_sd_boot=np.std(boots),
                         GAM20=gam.get((c, 20), np.nan), GAM30=gam.get((c, 30), np.nan),
                         frac_pos=(tal > 0).mean()))
    return pd.DataFrame(rows).sort_values('ddG_GAM', ascending=False)


def summarize(cfg, out_dir=None, finished=None, all_output=False, log=print):
    out_dir = out_dir or os.path.join(cfg['dir'], 'analysis_output')
    if finished is None:
        muts = None if all_output else slurm.selected_mutations(cfg)
        finished = finished_trajectories(cfg, muts)
    a = cfg['analysis']
    results, structs = analyze(finished, out_dir, a['min_nstruct'], cfg['protocol']['stride'], log)
    r = summary(results, structs, a['bootstrap'], a['seed'])
    r.round(3).to_csv(os.path.join(out_dir, 'alascan_summary.csv'), index=False)
    with pd.option_context('display.width', 200):
        log(r.round(2).to_string(index=False))
    return r
