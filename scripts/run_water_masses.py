#!/usr/bin/env python3
"""Ejecuta el analizador V3.4 con rutas independientes del equipo/campaña.

Uso:
    python scripts/run_water_masses.py --input <ruta/CTD_BASE_FINAL.xlsx> --output <carpeta>
"""
import argparse
import importlib.util
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Análisis de masas de agua v3.4 desde la base CTD validada')
    parser.add_argument('--input', required=True, type=Path, help='Excel CTD_BASE_FINAL.xlsx de v2.6')
    parser.add_argument('--output', default=Path('outputs/analisis_masas_agua'), type=Path,
                        help='Carpeta de resultados')
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f'No existe el Excel de entrada: {args.input}')
    script = Path(__file__).resolve().parents[1] / 'src' / 'water_masses_v3_4.py'
    spec = importlib.util.spec_from_file_location('ctd_water_masses_v3_4', script)
    module = importlib.util.module_from_spec(spec)
    import sys
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.CFG['archivo_entrada'] = args.input.resolve()
    module.CFG['carpeta_salida'] = args.output.resolve()
    module.main()

if __name__ == '__main__':
    main()
