"""Locate the published, unchanged ElPiGraph engine; no result files are read."""
from pathlib import Path
import importlib.util
import sys

HERE = Path(__file__).resolve().parent
SITE = HERE.parent
ROOT = SITE / 'chronology-stability/replication'
if not ROOT.is_dir():
    ROOT = HERE.parents[2]
TEAM = ROOT / 'solve/botanical_enrichment/team_golf'
STUDY = TEAM / 'golf_published_stability_20261016'
RF = TEAM / 'golf_rf1b_stability_20261016'
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(STUDY))


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def preprocessing():
    return module('text_order_preprocessing', RF / 'prepare_data.py')


def engine():
    import ensemble
    return ensemble
