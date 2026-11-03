"""Repository-level syntax and version checks; no data files needed."""
import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def test_ctd_script_compiles():
    p = ROOT / 'src' / 'ctd_processor_v2_6.py'
    ast.parse(p.read_text(encoding='utf-8'), filename=str(p))

def test_water_mass_script_compiles():
    p = ROOT / 'src' / 'water_masses_v3_4.py'
    ast.parse(p.read_text(encoding='utf-8'), filename=str(p))

def test_cli_compiles():
    p = ROOT / 'scripts' / 'run_water_masses.py'
    ast.parse(p.read_text(encoding='utf-8'), filename=str(p))

def test_versions_present():
    assert '2.6-DISCOVERY-VF-ESTATICOS' in (ROOT / 'src' / 'ctd_processor_v2_6.py').read_text(encoding='utf-8')
    assert 'VERSION = "3.4"' in (ROOT / 'src' / 'water_masses_v3_4.py').read_text(encoding='utf-8')

def test_no_private_old_windows_path():
    assert 'D:\\lorenzo' not in (ROOT / 'src' / 'water_masses_v3_4.py').read_text(encoding='utf-8')
