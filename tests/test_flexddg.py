'''Regression and unit tests. Run with `python3 tests/test_flexddg.py` (pytest also works if installed).

The golden files in tests/golden/ were frozen on 2026-10-26 from the original one-off scripts' output
(TROP2 setup and summary) and from the inputs the first CLDN array job ran with.
The summary regression reads all 875 TROP2 trajectories and takes about a minute.
'''
import glob, hashlib, json, os, shutil, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import pandas as pd  # noqa: E402

from flexddg import analysis, config, inputs, runner  # noqa: E402

GOLD = os.path.join(HERE, 'golden')
TARGETS = {'trop2_9pi9': '9pi9', 'cldn18_9v32': '9v32'}


def temp_target(name, **overrides):
    '''A temp target dir with the real target.json (pdb made absolute), scanning both sides.'''
    src = os.path.join(ROOT, name)
    cfg = json.load(open(os.path.join(src, 'target.json')))
    cfg['pdb'] = os.path.normpath(os.path.join(src, cfg['pdb']))
    cfg['scan'] = dict(cfg.get('scan', {}), side=list(cfg['groups']))
    cfg.update(overrides)
    d = tempfile.mkdtemp(prefix='flexddg_test_', dir=os.environ.get('SCRATCH'))
    json.dump(cfg, open(os.path.join(d, 'target.json'), 'w'))
    return d


def inputs_listing(d, muts):
    out = []
    for m in muts:
        out.append('== %s\n' % m)
        for f in ['chains_to_move.txt', 'nataa_mutations.resfile']:
            out.append(open(os.path.join(d, 'inputs', m, f)).read())
        out.append(os.path.basename(glob.glob(os.path.join(d, 'inputs', m, '*.pdb'))[0]) + '\n')
    return ''.join(out)


def check_setup_matches_golden(name):
    d = temp_target(name)
    try:
        inputs.setup(d, log=lambda *a: None)
        g = os.path.join(GOLD, name)
        for f in ['interface_residues.txt', 'mutations.txt', 'mutations_fab.txt']:
            assert open(os.path.join(d, f)).read() == open(os.path.join(g, f)).read(), f
        muts = open(os.path.join(d, 'mutations.txt')).read().split()
        assert inputs_listing(d, muts) == open(os.path.join(g, 'inputs.txt')).read()
        clean = os.path.join(d, TARGETS[name] + '_clean.pdb')
        assert hashlib.md5(open(clean, 'rb').read()).hexdigest() == open(os.path.join(g, 'clean_pdb.md5')).read().strip()
        for m in muts:
            manifest = json.load(open(os.path.join(d, 'inputs', m, 'input.json')))
            assert os.path.exists(os.path.join(d, 'inputs', m, manifest['pdb']))
    finally:
        shutil.rmtree(d)


def test_setup_trop2_matches_original_scripts():
    check_setup_matches_golden('trop2_9pi9')


def test_setup_cldn_matches_running_inputs():
    check_setup_matches_golden('cldn18_9v32')


def test_side_selection_and_custom_mutations():
    d = temp_target('cldn18_9v32', mutations=['H:W33F', 'H:W33A+L:Y100A'])
    try:
        cfg = json.load(open(os.path.join(d, 'target.json')))
        cfg['scan']['side'] = ['fab']
        json.dump(cfg, open(os.path.join(d, 'target.json'), 'w'))
        _, names = inputs.setup(d, log=lambda *a: None)
        fab = open(os.path.join(GOLD, 'cldn18_9v32', 'mutations_fab.txt')).read().split()
        assert names == fab + ['H_W33F', 'H_W33A-L_Y100A']
        assert open(os.path.join(d, 'inputs', 'H_W33A-L_Y100A', 'nataa_mutations.resfile')).read() == \
            'NATAA\nstart\n33 H PIKAA A\n100 L PIKAA A\n'
        assert open(os.path.join(d, 'mutations_custom.txt')).read() == 'H_W33F\nH_W33A-L_Y100A\n'
    finally:
        shutil.rmtree(d)


def test_wrong_wildtype_is_rejected():
    d = temp_target('cldn18_9v32', mutations=['H:Y33A'])   # H33 is TRP
    try:
        inputs.setup(d, log=lambda *a: None)
        raise AssertionError('expected ConfigError')
    except config.ConfigError as e:
        assert 'is TRP' in str(e)
    finally:
        shutil.rmtree(d)


def test_config_validation():
    base = {'pdb': 'x.pdb', 'groups': {'fab': ['H', 'L'], 'antigen': ['A']}, 'move': 'antigen'}
    bad = [({'move': 'nope'}, 'not a group'),
           ({'groups': {'fab': ['H'], 'antigen': ['H']}}, 'both groups'),
           ({'colour': 1}, 'unknown key'),
           ({'protocol': {'nstruct': 0}}, 'positive integer'),
           ({'protocol': {'ntrials': 1000, 'stride': 300}}, 'multiple'),
           ({'scan': {'side': ['x']}}, 'not a group'),
           ({'mutations': ['H33A']}, 'bad mutation'),
           ({'mutations': ['Q:Y33A']}, 'not in any group')]
    for override, msg in bad:
        d = tempfile.mkdtemp()
        try:
            json.dump(dict(base, **override), open(os.path.join(d, 'target.json'), 'w'))
            try:
                config.load(d)
                raise AssertionError('accepted %r' % override)
            except config.ConfigError as e:
                assert msg in str(e), (override, e)
        finally:
            shutil.rmtree(d)
    d = tempfile.mkdtemp()
    json.dump(base, open(os.path.join(d, 'target.json'), 'w'))
    cfg = config.load(d)
    assert cfg['scan']['side'] == ['fab', 'antigen'] and cfg['analysis']['min_nstruct'] == 35
    shutil.rmtree(d)


def test_mutation_names():
    assert config.mutation_name(config.parse_mutation('H:Y33A')) == 'H_Y33A'
    assert config.mutation_name(config.parse_mutation('H:S52AA+L:Y100F')) == 'H_S52AA-L_Y100F'
    assert config.parse_mutation('H:S52AA') == [('H', 'S', 52, 'A', 'A')]
    assert inputs.resfile(config.parse_mutation('H:S52AA')) == 'NATAA\nstart\n52A H PIKAA A\n'


def test_rosetta_command_matches_legacy_run_one():
    '''The flags the old pipeline/run_one.sh passed, so new and old trajectories are interchangeable.'''
    cfg = config.load(os.path.join(ROOT, 'cldn18_9v32'))
    ind = os.path.join(cfg['dir'], 'inputs', 'L_Y100A')
    cmd = runner.command(cfg, ind, 7, cfg['protocol'])
    assert cmd == ['rosetta_scripts.default.linuxgccrelease', '-s', os.path.join(ind, '9v32_clean.pdb'),
                   '-parser:protocol', runner.XML, '-parser:script_vars', 'chainstomove=A',
                   'mutate_resfile_relpath=' + os.path.join(ind, 'nataa_mutations.resfile'),
                   'number_backrub_trials=35000', 'max_minimization_iter=5000', 'abs_score_convergence_thresh=1.0',
                   'backrub_trajectory_stride=35000', '-restore_talaris_behavior', '-in:file:fullatom',
                   '-ignore_unrecognized_res', '-ignore_zero_occupancy', 'false', '-ex1', '-ex2', '-constant_seed',
                   '-jran', '1007']


def test_summary_reproduces_trop2():
    cfg = config.load(os.path.join(ROOT, 'trop2_9pi9'))
    out = tempfile.mkdtemp(prefix='flexddg_test_', dir=os.environ.get('SCRATCH'))
    try:
        fab = open(os.path.join(GOLD, 'trop2_9pi9', 'mutations_fab.txt')).read().split()
        r = analysis.summarize(cfg, out_dir=out, finished=analysis.finished_trajectories(cfg, fab), log=lambda *a: None)
        new = pd.read_csv(os.path.join(out, 'alascan_summary.csv')).set_index('case').sort_index()
        old = pd.read_csv(os.path.join(GOLD, 'trop2_9pi9', 'fab_alascan_summary.csv')).set_index('case').sort_index()
        assert list(new.columns) == list(old.columns) and list(new.index) == list(old.index)
        assert (new - old).abs().max().max() < 1e-9
        assert len(r) == 25
        cols = ['case_name', 'backrub_steps', 'score_function_name', 'scored_state', 'nstruct']
        key = lambda df: df.drop(columns='Unnamed: 0').sort_values(cols).reset_index(drop=True)
        new_res = key(pd.read_csv(os.path.join(out, 'results.csv')))
        old_res = key(pd.read_csv(os.path.join(GOLD, 'trop2_9pi9', 'fab_complete-results.csv')))
        pd.testing.assert_frame_equal(new_res, old_res, check_exact=False, rtol=0, atol=1e-9)
    finally:
        shutil.rmtree(out)


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted((k, v) for k, v in dict(globals()).items() if k.startswith('test_')):
        try:
            fn()
            print('PASS', name, flush=True)
        except Exception as e:  # noqa: BLE001
            failed += 1
            print('FAIL', name, repr(e), flush=True)
    sys.exit(1 if failed else 0)
