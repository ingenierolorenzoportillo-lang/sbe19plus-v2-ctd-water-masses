# -*- coding: utf-8 -*-
"""
ANÁLISIS DE MASAS DE AGUA — CARIBE COLOMBIANO
==============================================
Versión 3.4 (historial de versiones integrado); diseñada para CTD_BASE_FINAL.xlsx.

Objetivo
--------
Leer CTD_BASE_FINAL.xlsx y producir una rutina reproducible para:
  1) describir Caribbean Surface Water (CSW) mediante MLD/ILD;
  2) detectar el núcleo de North Atlantic Subtropical Underwater (SUW)
     como máximo de Salinidad Absoluta por debajo de 70 m;
  3) detectar el núcleo de Antarctic Intermediate Water (AAIW)
     como mínimo de Salinidad Absoluta por debajo del núcleo SUW;
  4) estimar fracciones mSUW–mAAIW–mNADW mediante OMP clásico
     usando Θ, SA, O2 y conservación de masa;
  5) generar tablas y figuras orientadas a un manuscrito científico.

Referencia metodológica principal
---------------------------------
Torres, R. R., Latandret, S., Salon, J., & Dagua, C. (2023).
Water masses in the Caribbean Sea and sub-annual variability in the
Guajira upwelling region. Ocean Dynamics, 73, 39–57.
https://doi.org/10.1007/s10236-022-01529-5

Criterios reproducidos de Torres et al. (2023):
  - MLD: incremento de sigma0 = 0.03 kg m-3 respecto a 10 m.
  - ILD: diferencia de Θ = 0.2 °C respecto a 10 m.
  - SUW core: máximo de SA por debajo de 70 m.
  - OMP debajo del máximo de salinidad.
  - SWT modificados del Caribe y pesos OMP de su Tabla 1.
  - contribución >= 30 % como umbral descriptivo para capas de masas de agua.
  - residual de conservación de masa <= 10 % como control del ajuste OMP.

IMPORTANTE
----------
Este script NO convierte los porcentajes obtenidos en "porcentaje del volumen
del Caribe". Los perfiles del crucero no constituyen una malla volumétrica
basin-wide. Las fracciones OMP se interpretan localmente, por estación y
profundidad, dentro del dominio muestreado.

El OMP clásico normaliza cada propiedad usando media/desviación estándar de
los source water types y aplica los pesos publicados antes de resolver el
sistema por mínimos cuadrados no negativos. Se conserva el residual de masa
sin forzar artificialmente una suma exacta de 1. Para visualización se crean
fracciones normalizadas únicamente cuando el residual de masa es aceptable.

Requisitos
----------
    numpy pandas scipy matplotlib openpyxl gsw

Ejecutar completo con F5 en Spyder.
"""

from __future__ import annotations

import json
import math
import re
import sys
import warnings
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.optimize import lsq_linear

try:
    import gsw
except ImportError as exc:
    raise SystemExit(
        "Falta gsw. En la consola IPython de Spyder ejecuta: %pip install gsw"
    ) from exc


# =============================================================================
# CONFIGURACIÓN
# =============================================================================

VERSION = "1.1"

DEFAULT_SESSION = Path.cwd()

CFG = {
    "archivo_entrada": None,  # Elegir por consola o ejecutar scripts/run_water_masses.py
    "hoja_base": "Base_CTD",
    "hoja_decisiones": "Decisiones_auto",
    "carpeta_salida": DEFAULT_SESSION / "ANALISIS_MASAS_AGUA_VF",

    # QC del producto CTD v2.2. Se conservan PASS y SUSPECT porque un valor
    # sospechoso puede corresponder a estructura oceanográfica real.
    "usar_qartod": True,
    "qartod_aceptados": (1, 3),

    # Para OMP con oxígeno, por defecto se usan solamente lances cuyo estado
    # dinámico O2 haya sido AUTO_OK. Puede cambiarse después de revisión.
    "omp_requiere_O2_AUTO_OK": True,

    # Torres et al. (2023)
    "ref_depth_m": 10.0,
    "mld_delta_sigma0": 0.03,
    "ild_delta_CT_C": 0.20,
    "suw_min_depth_m": 70.0,
    "omp_mass_residual_max_pct": 10.0,
    "fraccion_significativa": 0.30,

    # Si el mínimo AAIW aparece muy cerca del fondo muestreado no se interpreta
    # como núcleo plenamente resuelto.
    "aaiw_min_depth_m": 350.0,
    "aaiw_bottom_margin_m": 60.0,

    # Secciones. None = ordenar por número de estación extraído del nombre.
    # Ejemplo manual: [1, 2, 3, 4, 5, 6, 7, 9]
    "orden_estaciones": None,
    "seccion_dz_m": 10.0,

    # Profundidades de interés comparables con el trabajo de Torres et al.
    "profundidades_estandar_m": [150, 500, 700, 1200, 1450],
    "tolerancia_prof_estandar_m": 15.0,

    "dpi": 300,
    "mostrar_figuras": False,

    # Figuras VF: PASS es la capa principal; SUSPECT se conserva como
    # diagnóstico visual, pero no controla los límites de ejes.
    "qartod_pass": 1,
    "qartod_suspect": 3,

    # Secciones espaciales propuestas según la numeración de estaciones.
    # IMPORTANTE: confirmar contra el diseño/bitácora del crucero antes de
    # denominarlas transectos definitivos en el manuscrito.
    "transectos": {
        "Norte": [1, 2, 3],
        "Centro": [4, 5, 6],
        "Sur": [7, 8, 9],
    },

    # Diagnóstico OMP. Estos NO son umbrales universales; son alarmas
    # operativas para detectar desacople sistemático de propiedades.
    "omp_o2_bias_warn_umol_kg": 10.0,
    "omp_same_sign_warn_fraction": 0.80,
    "omp_sensitivity_warn_fraction": 0.15,

    # Guardar también PDF vectorial además de PNG.
    "guardar_pdf": True,

    # Contexto cartográfico opcional. Si Basemap está instalado se dibuja
    # costa y relieve ETOPO de contexto. No es requisito del análisis.
    "mapa_basemap": True,
}


# Source Water Types modificados para el Caribe y pesos, Tabla 1 de
# Torres et al. (2023), publicación corregida.
SWT = pd.DataFrame(
    {
        "mSUW": {"CT_C": 23.24, "SA_g_kg": 37.30, "O2_umol_kg": 182.8},
        "mAAIW": {"CT_C": 5.96, "SA_g_kg": 34.84, "O2_umol_kg": 135.3},
        "mNADW": {"CT_C": 3.87, "SA_g_kg": 35.15, "O2_umol_kg": 241.6},
    }
)

OMP_WEIGHTS = {
    "CT_C": 4.142,
    "SA_g_kg": 6.175,
    "O2_umol_kg": 1.678,
    "mass": 6.0,
}

SOURCE_REFERENCE = {
    "title": "Water masses in the Caribbean Sea and sub-annual variability in the Guajira upwelling region",
    "authors": "Torres, Latandret, Salon & Dagua",
    "journal": "Ocean Dynamics 73, 39-57",
    "year": 2023,
    "doi": "10.1007/s10236-022-01529-5",
}


# =============================================================================
# UTILIDADES
# =============================================================================


def slug(text: str) -> str:
    s = str(text).strip().lower()
    s = re.sub(r"[^a-z0-9áéíóúüñ]+", "_", s)
    return s.strip("_") or "perfil"


def periodo_from_station(name: str) -> str:
    s = str(name).lower()
    if "noct" in s or "noche" in s:
        return "Nocturno"
    if "diurn" in s or re.search(r"\bdia\b", s):
        return "Diurno"
    return "Sin_clasificar"


def station_number(name: str) -> float:
    m = re.search(r"(?:estaci[oó]n|e)?\s*0*(\d+)", str(name), flags=re.I)
    if not m:
        nums = re.findall(r"\d+", str(name))
        return float(nums[0]) if nums else np.nan
    return float(m.group(1))


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    p1 = np.deg2rad(lat1)
    p2 = np.deg2rad(lat2)
    dphi = np.deg2rad(lat2 - lat1)
    dlambda = np.deg2rad(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlambda / 2) ** 2
    return 2 * r * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def safe_interp(x, y, x0):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 2:
        return np.nan
    xx, yy = x[m], y[m]
    o = np.argsort(xx)
    xx, yy = xx[o], yy[o]
    if x0 < xx.min() or x0 > xx.max():
        return np.nan
    return float(np.interp(x0, xx, yy))


def first_threshold_depth(depth, values, threshold, mode="ge"):
    """Profundidad interpolada del primer cruce desde superficie."""
    z = np.asarray(depth, float)
    v = np.asarray(values, float)
    m = np.isfinite(z) & np.isfinite(v)
    if m.sum() < 2:
        return np.nan
    z, v = z[m], v[m]
    o = np.argsort(z)
    z, v = z[o], v[o]
    hit = v >= threshold if mode == "ge" else v <= threshold
    idx = np.flatnonzero(hit)
    if len(idx) == 0:
        return np.nan
    j = int(idx[0])
    if j == 0:
        return float(z[0])
    z1, z2 = z[j - 1], z[j]
    v1, v2 = v[j - 1], v[j]
    if not np.isfinite(v2 - v1) or abs(v2 - v1) < 1e-12:
        return float(z2)
    f = (threshold - v1) / (v2 - v1)
    f = float(np.clip(f, 0, 1))
    return float(z1 + f * (z2 - z1))


def choose_input_file(path: Path) -> Path:
    if path.exists():
        return path
    print(f"No se encontró la ruta por defecto:\n  {path}")
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        sel = filedialog.askopenfilename(
            title="Selecciona CTD_BASE_FINAL.xlsx",
            filetypes=[("Excel", "*.xlsx"), ("Todos", "*.*")],
        )
        root.destroy()
        if sel:
            return Path(sel)
    except Exception:
        pass
    raise FileNotFoundError(
        "No se encontró CTD_BASE_FINAL.xlsx. Edita CFG['archivo_entrada'] con la ruta correcta."
    )


def savefig(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=CFG["dpi"], bbox_inches="tight", facecolor="white")
    if CFG["mostrar_figuras"]:
        plt.show()
    else:
        plt.close(fig)


# =============================================================================
# CARGA Y NORMALIZACIÓN DEL DATASET
# =============================================================================


def find_o2_column(df: pd.DataFrame) -> str | None:
    cand = [c for c in df.columns if str(c).startswith("oxigeno_") and str(c).endswith("_umol_kg")]
    if not cand:
        return None
    if len(cand) == 1:
        return cand[0]
    counts = {c: int(pd.to_numeric(df[c], errors="coerce").notna().sum()) for c in cand}
    chosen = max(counts, key=counts.get)
    warnings.warn(f"Se detectaron varios O2 en µmol/kg: {cand}. Se usará {chosen} por mayor cobertura.")
    return chosen


def cargar_dataset(path: Path):
    xls = pd.ExcelFile(path)
    if CFG["hoja_base"] not in xls.sheet_names:
        raise ValueError(f"No existe la hoja {CFG['hoja_base']} en {path.name}")

    base = pd.read_excel(path, sheet_name=CFG["hoja_base"])
    decisiones = (
        pd.read_excel(path, sheet_name=CFG["hoja_decisiones"])
        if CFG["hoja_decisiones"] in xls.sheet_names
        else pd.DataFrame()
    )

    req = {"estacion", "latitud", "longitud", "salinidad_absoluta_g_kg", "temperatura_conservativa_C", "sigma0_kg_m3"}
    faltan = req - set(base.columns)
    if faltan:
        raise ValueError("Faltan columnas esenciales en Base_CTD: " + ", ".join(sorted(faltan)))

    pcol = "bin_centro_dbar" if "bin_centro_dbar" in base else "presion_dbar"
    if pcol not in base:
        raise ValueError("No se encontró presión ni bin_centro_dbar.")

    if "profundidad_m" not in base:
        if not {"latitud", pcol}.issubset(base.columns):
            raise ValueError("No hay profundidad_m y tampoco latitud para reconstruirla con TEOS-10.")
        base["profundidad_m"] = -gsw.z_from_p(
            pd.to_numeric(base[pcol], errors="coerce").to_numpy(float),
            pd.to_numeric(base["latitud"], errors="coerce").to_numpy(float),
        )

    if "direccion" in base:
        base = base[base["direccion"].astype(str).str.lower().eq("bajada")].copy()

    # Elimina bins sin T/S/SA/CT físicamente analizables.
    for c in [pcol, "profundidad_m", "salinidad_absoluta_g_kg", "temperatura_conservativa_C", "sigma0_kg_m3"]:
        base[c] = pd.to_numeric(base[c], errors="coerce")

    o2col = find_o2_column(base)
    if o2col:
        base[o2col] = pd.to_numeric(base[o2col], errors="coerce")

    base["periodo"] = base["estacion"].map(periodo_from_station)
    base["numero_estacion"] = base["estacion"].map(station_number)

    profile_col = "id_perfil" if "id_perfil" in base else "estacion"

    # Estado operativo por estación, para decidir si O2 se usa en OMP.
    estado_o2 = {}
    if not decisiones.empty and {"estacion", "estado_O2"}.issubset(decisiones.columns):
        estado_o2 = dict(zip(decisiones["estacion"].astype(str), decisiones["estado_O2"].astype(str)))
    base["estado_O2_perfil"] = base["estacion"].astype(str).map(estado_o2).fillna("NO_EVALUADO")

    return base, decisiones, pcol, o2col, profile_col


def filtro_hidrografico(df: pd.DataFrame) -> pd.DataFrame:
    x = df.copy()
    if CFG["usar_qartod"] and "qartod_resumen_ctd" in x:
        x = x[x["qartod_resumen_ctd"].isin(CFG["qartod_aceptados"])].copy()
    req = ["profundidad_m", "salinidad_absoluta_g_kg", "temperatura_conservativa_C", "sigma0_kg_m3"]
    x = x[x[req].notna().all(axis=1)].copy()
    return x


def filtro_omp(df: pd.DataFrame, o2col: str | None) -> pd.DataFrame:
    if o2col is None:
        return df.iloc[0:0].copy()
    x = filtro_hidrografico(df)
    if CFG["usar_qartod"] and "qartod_oxigeno" in x:
        x = x[x["qartod_oxigeno"].isin(CFG["qartod_aceptados"])].copy()
    if CFG["omp_requiere_O2_AUTO_OK"]:
        x = x[x["estado_O2_perfil"].eq("AUTO_OK")].copy()
    x = x[x[o2col].notna()].copy()
    return x


# =============================================================================
# MLD, ILD, CSW Y NÚCLEOS HIDROGRÁFICOS
# =============================================================================


def analizar_perfil_hidrografico(g: pd.DataFrame, pcol: str, o2col: str | None) -> dict:
    x = filtro_hidrografico(g).sort_values("profundidad_m").copy()
    out = {
        "estacion": str(g["estacion"].iloc[0]),
        "periodo": str(g["periodo"].iloc[0]),
        "numero_estacion": float(g["numero_estacion"].iloc[0]),
        "latitud": float(pd.to_numeric(g.get("latitud", pd.Series([np.nan])).iloc[0], errors="coerce")),
        "longitud": float(pd.to_numeric(g.get("longitud", pd.Series([np.nan])).iloc[0], errors="coerce")),
        "profundidad_max_m": float(pd.to_numeric(x["profundidad_m"], errors="coerce").max()) if len(x) else np.nan,
        "n_bins_hidrografia": int(len(x)),
    }
    if len(x) < 5:
        out["estado_nucleos"] = "DATOS_INSUFICIENTES"
        return out

    z = x["profundidad_m"].to_numpy(float)
    sig = x["sigma0_kg_m3"].to_numpy(float)
    ct = x["temperatura_conservativa_C"].to_numpy(float)
    sa = x["salinidad_absoluta_g_kg"].to_numpy(float)
    rho = x["densidad_kg_m3"].to_numpy(float) if "densidad_kg_m3" in x else (1000.0 + sig)

    zref = CFG["ref_depth_m"]
    sig10 = safe_interp(z, sig, zref)
    ct10 = safe_interp(z, ct, zref)
    sa10 = safe_interp(z, sa, zref)

    if np.isfinite(sig10):
        mm = z >= zref
        mld = first_threshold_depth(z[mm], sig[mm] - sig10, CFG["mld_delta_sigma0"], mode="ge")
        if not np.isfinite(mld):
            mld = zref
    else:
        mld = np.nan

    if np.isfinite(ct10):
        mm = z >= zref
        ild = first_threshold_depth(z[mm], np.abs(ct[mm] - ct10), CFG["ild_delta_CT_C"], mode="ge")
        if not np.isfinite(ild):
            ild = zref
    else:
        ild = np.nan

    out.update(
        {
            "sigma0_10m": sig10,
            "CT_10m_C": ct10,
            "SA_10m_g_kg": sa10,
            "MLD_m": mld,
            "ILD_m": ild,
            "BLT_m": max(0.0, ild - mld) if np.isfinite(ild) and np.isfinite(mld) else np.nan,
            "CSW_CT_MLD_C": safe_interp(z, ct, mld) if np.isfinite(mld) else np.nan,
            "CSW_SA_MLD_g_kg": safe_interp(z, sa, mld) if np.isfinite(mld) else np.nan,
            "CSW_sigma0_MLD": safe_interp(z, sig, mld) if np.isfinite(mld) else np.nan,
        }
    )

    # Núcleo SUW: máximo de SA por debajo de 70 m (Torres et al. 2023).
    suw = x[x["profundidad_m"] >= CFG["suw_min_depth_m"]].copy()
    if len(suw):
        i = suw["salinidad_absoluta_g_kg"].idxmax()
        r = suw.loc[i]
        suw_z = float(r["profundidad_m"])
        out.update(
            {
                "SUW_core_depth_m": suw_z,
                "SUW_core_SA_g_kg": float(r["salinidad_absoluta_g_kg"]),
                "SUW_core_CT_C": float(r["temperatura_conservativa_C"]),
                "SUW_core_sigma0": float(r["sigma0_kg_m3"]),
                "SUW_core_O2_umol_kg": float(r[o2col]) if o2col and pd.notna(r.get(o2col, np.nan)) else np.nan,
            }
        )
    else:
        suw_z = np.nan
        out.update({k: np.nan for k in [
            "SUW_core_depth_m", "SUW_core_SA_g_kg", "SUW_core_CT_C", "SUW_core_sigma0", "SUW_core_O2_umol_kg"
        ]})

    # Estabilidad estática entre MLD y núcleo SUW.
    if np.isfinite(mld) and np.isfinite(suw_z) and suw_z > mld:
        rho_mld = safe_interp(z, rho, mld)
        rho_suw = safe_interp(z, rho, suw_z)
        if np.isfinite(rho_mld) and np.isfinite(rho_suw):
            # z oceanográfico negativo: E = -(1/rho) d rho / dz.
            z1, z2 = -mld, -suw_z
            out["estabilidad_MLD_SUW_m_1"] = -(rho_suw - rho_mld) / (((rho_suw + rho_mld) / 2) * (z2 - z1))
        else:
            out["estabilidad_MLD_SUW_m_1"] = np.nan
    else:
        out["estabilidad_MLD_SUW_m_1"] = np.nan

    # Núcleo AAIW: mínimo de SA por debajo del núcleo SUW, con umbral mínimo
    # operacional de 350 m para evitar confundir Central Water con AAIW.
    if np.isfinite(suw_z):
        zmin = max(CFG["aaiw_min_depth_m"], suw_z + 50.0)
    else:
        zmin = CFG["aaiw_min_depth_m"]
    aaiw = x[x["profundidad_m"] >= zmin].copy()
    if len(aaiw):
        i = aaiw["salinidad_absoluta_g_kg"].idxmin()
        r = aaiw.loc[i]
        aaiw_z = float(r["profundidad_m"])
        bottom = float(x["profundidad_m"].max())
        trunc = (bottom - aaiw_z) <= CFG["aaiw_bottom_margin_m"]
        out.update(
            {
                "AAIW_core_depth_m": aaiw_z,
                "AAIW_core_SA_g_kg": float(r["salinidad_absoluta_g_kg"]),
                "AAIW_core_CT_C": float(r["temperatura_conservativa_C"]),
                "AAIW_core_sigma0": float(r["sigma0_kg_m3"]),
                "AAIW_core_O2_umol_kg": float(r[o2col]) if o2col and pd.notna(r.get(o2col, np.nan)) else np.nan,
                "AAIW_core_truncado_fondo": bool(trunc),
            }
        )
    else:
        out.update({k: np.nan for k in [
            "AAIW_core_depth_m", "AAIW_core_SA_g_kg", "AAIW_core_CT_C", "AAIW_core_sigma0", "AAIW_core_O2_umol_kg"
        ]})
        out["AAIW_core_truncado_fondo"] = True

    out["CSW_rango_climatologico_Torres"] = bool(
        np.isfinite(out.get("CSW_CT_MLD_C", np.nan))
        and np.isfinite(out.get("CSW_SA_MLD_g_kg", np.nan))
        and 25.0 <= out["CSW_CT_MLD_C"] <= 29.0
        and 35.0 <= out["CSW_SA_MLD_g_kg"] <= 36.5
    )
    out["estado_nucleos"] = "OK"
    return out


# =============================================================================
# OMP CLÁSICO: mSUW + mAAIW + mNADW
# =============================================================================


@dataclass
class OMPSolution:
    frac_mSUW_raw: float = np.nan
    frac_mAAIW_raw: float = np.nan
    frac_mNADW_raw: float = np.nan
    sum_frac_raw: float = np.nan
    mass_residual_pct: float = np.nan
    resid_CT_C: float = np.nan
    resid_SA_g_kg: float = np.nan
    resid_O2_umol_kg: float = np.nan
    weighted_rms: float = np.nan
    omp_ok: bool = False
    frac_mSUW: float = np.nan
    frac_mAAIW: float = np.nan
    frac_mNADW: float = np.nan
    masa_dominante: str = "NO_EVALUADO"


def _omp_normalization():
    means = SWT.mean(axis=1)
    stds = SWT.std(axis=1, ddof=1)  # MATLAB std() usado por OMP v2 normaliza con N-1
    if (stds <= 0).any():
        raise ValueError("Alguna propiedad SWT tiene desviación estándar cero; no puede normalizarse.")
    return means, stds


def omp_solve(ct_obs, sa_obs, o2_obs) -> OMPSolution:
    if not all(np.isfinite([ct_obs, sa_obs, o2_obs])):
        return OMPSolution()

    means, stds = _omp_normalization()

    # Matriz SWT normalizada: filas CT, SA, O2; columnas mSUW, mAAIW, mNADW.
    G = np.vstack([
        (SWT.loc["CT_C"].to_numpy(float) - means["CT_C"]) / stds["CT_C"],
        (SWT.loc["SA_g_kg"].to_numpy(float) - means["SA_g_kg"]) / stds["SA_g_kg"],
        (SWT.loc["O2_umol_kg"].to_numpy(float) - means["O2_umol_kg"]) / stds["O2_umol_kg"],
        np.ones(3),
    ])

    d = np.array([
        (ct_obs - means["CT_C"]) / stds["CT_C"],
        (sa_obs - means["SA_g_kg"]) / stds["SA_g_kg"],
        (o2_obs - means["O2_umol_kg"]) / stds["O2_umol_kg"],
        1.0,
    ])
    w = np.array([
        OMP_WEIGHTS["CT_C"], OMP_WEIGHTS["SA_g_kg"],
        OMP_WEIGHTS["O2_umol_kg"], OMP_WEIGHTS["mass"],
    ], float)

    Gw = G * w[:, None]
    dw = d * w

    # No-negatividad física; no se impone x<=1 porque la ecuación de masa ya
    # controla la suma y queremos conservar su residual como diagnóstico.
    sol = lsq_linear(Gw, dw, bounds=(0.0, np.inf), method="trf", lsmr_tol="auto")
    if not sol.success:
        return OMPSolution()

    f = np.asarray(sol.x, float)
    s = float(f.sum())
    mass_resid = 100.0 * abs(s - 1.0)

    # Reconstrucción en unidades físicas para residuales interpretables.
    ct_rec = float(np.dot(SWT.loc["CT_C"].to_numpy(float), f))
    sa_rec = float(np.dot(SWT.loc["SA_g_kg"].to_numpy(float), f))
    o2_rec = float(np.dot(SWT.loc["O2_umol_kg"].to_numpy(float), f))
    res_ct = ct_obs - ct_rec
    res_sa = sa_obs - sa_rec
    res_o2 = o2_obs - o2_rec
    wrms = float(np.sqrt(np.mean((Gw @ f - dw) ** 2)))

    ok = bool(np.isfinite(mass_resid) and mass_resid <= CFG["omp_mass_residual_max_pct"] and s > 0)
    if ok:
        fn = f / s
        names = ["mSUW", "mAAIW", "mNADW"]
        j = int(np.argmax(fn))
        dominant = names[j] if fn[j] >= CFG["fraccion_significativa"] else "MEZCLA"
    else:
        fn = np.full(3, np.nan)
        dominant = "OMP_NO_VALIDO"

    return OMPSolution(
        frac_mSUW_raw=float(f[0]),
        frac_mAAIW_raw=float(f[1]),
        frac_mNADW_raw=float(f[2]),
        sum_frac_raw=s,
        mass_residual_pct=mass_resid,
        resid_CT_C=float(res_ct),
        resid_SA_g_kg=float(res_sa),
        resid_O2_umol_kg=float(res_o2),
        weighted_rms=wrms,
        omp_ok=ok,
        frac_mSUW=float(fn[0]),
        frac_mAAIW=float(fn[1]),
        frac_mNADW=float(fn[2]),
        masa_dominante=dominant,
    )


def ejecutar_omp(base: pd.DataFrame, nucleos: pd.DataFrame, o2col: str | None, profile_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if o2col is None:
        return pd.DataFrame(), pd.DataFrame([
            {"motivo": "No existe columna de oxígeno en µmol kg-1; no se ejecutó OMP de tres end-members."}
        ])

    suw_depth = dict(zip(nucleos["estacion"].astype(str), nucleos["SUW_core_depth_m"]))
    rows = []
    excl = []

    for pid, g in base.groupby(profile_col, sort=False):
        est = str(g["estacion"].iloc[0])
        z_suw = pd.to_numeric(pd.Series([suw_depth.get(est, np.nan)]), errors="coerce").iloc[0]
        if not np.isfinite(z_suw):
            excl.append({"estacion": est, "id_perfil": pid, "motivo": "Sin núcleo SUW válido"})
            continue

        gg = filtro_omp(g, o2col)
        if gg.empty:
            estado = str(g["estado_O2_perfil"].iloc[0]) if "estado_O2_perfil" in g else "desconocido"
            excl.append({
                "estacion": est,
                "id_perfil": pid,
                "motivo": f"Sin bins OMP válidos o estado O2 no aceptado ({estado})",
            })
            continue

        # Torres et al.: OMP por debajo del máximo de salinidad.
        gg = gg[gg["profundidad_m"] >= float(z_suw)].copy()
        if gg.empty:
            excl.append({"estacion": est, "id_perfil": pid, "motivo": "Sin datos bajo el máximo de salinidad"})
            continue

        for _, r in gg.iterrows():
            sol = omp_solve(
                float(r["temperatura_conservativa_C"]),
                float(r["salinidad_absoluta_g_kg"]),
                float(r[o2col]),
            )
            row = {
                "id_perfil": pid,
                "estacion": est,
                "periodo": r["periodo"],
                "numero_estacion": r["numero_estacion"],
                "latitud": r.get("latitud", np.nan),
                "longitud": r.get("longitud", np.nan),
                "presion_dbar": r.get("bin_centro_dbar", r.get("presion_dbar", np.nan)),
                "profundidad_m": r["profundidad_m"],
                "CT_C": r["temperatura_conservativa_C"],
                "SA_g_kg": r["salinidad_absoluta_g_kg"],
                "O2_umol_kg": r[o2col],
                "sigma0_kg_m3": r.get("sigma0_kg_m3", np.nan),
                "estado_O2_perfil": r.get("estado_O2_perfil", "NO_EVALUADO"),
                **asdict(sol),
            }
            rows.append(row)

    return pd.DataFrame(rows), pd.DataFrame(excl)


def resumen_omp_por_perfil(omp: pd.DataFrame) -> pd.DataFrame:
    if omp.empty:
        return pd.DataFrame()
    rows = []
    for (pid, est), g in omp.groupby(["id_perfil", "estacion"], sort=False):
        q = g[g["omp_ok"]].copy()
        row = {
            "id_perfil": pid,
            "estacion": est,
            "periodo": g["periodo"].iloc[0],
            "numero_estacion": g["numero_estacion"].iloc[0],
            "n_bins_omp": int(len(g)),
            "n_bins_omp_ok": int(len(q)),
            "residual_masa_med_pct": float(q["mass_residual_pct"].mean()) if len(q) else np.nan,
            "residual_masa_max_pct": float(q["mass_residual_pct"].max()) if len(q) else np.nan,
        }
        for mass in ("mSUW", "mAAIW", "mNADW"):
            c = f"frac_{mass}"
            if len(q) and q[c].notna().any():
                im = q[c].idxmax()
                fmax = float(q.loc[im, c])
                zmax = float(q.loc[im, "profundidad_m"])
                significant = q[q[c] >= CFG["fraccion_significativa"]]
                row[f"{mass}_max_frac"] = fmax
                row[f"{mass}_depth_max_m"] = zmax
                row[f"{mass}_presente_ge30"] = bool(fmax >= CFG["fraccion_significativa"])
                if len(significant):
                    row[f"{mass}_top_ge30_m"] = float(significant["profundidad_m"].min())
                    row[f"{mass}_bottom_ge30_m"] = float(significant["profundidad_m"].max())
                    # Espesor aproximado muestreado; para bins ~1 dbar es cercano a metros.
                    zz = np.sort(significant["profundidad_m"].dropna().unique())
                    row[f"{mass}_thickness_ge30_m"] = float(zz[-1] - zz[0]) if len(zz) > 1 else 0.0
                else:
                    row[f"{mass}_top_ge30_m"] = np.nan
                    row[f"{mass}_bottom_ge30_m"] = np.nan
                    row[f"{mass}_thickness_ge30_m"] = 0.0
            else:
                for suf in ["max_frac", "depth_max_m", "top_ge30_m", "bottom_ge30_m", "thickness_ge30_m"]:
                    row[f"{mass}_{suf}"] = np.nan
                row[f"{mass}_presente_ge30"] = False
        rows.append(row)
    return pd.DataFrame(rows)


def omp_en_profundidades_estandar(omp: pd.DataFrame) -> pd.DataFrame:
    if omp.empty:
        return pd.DataFrame()
    rows = []
    tol = CFG["tolerancia_prof_estandar_m"]
    for (pid, est), g in omp[omp["omp_ok"]].groupby(["id_perfil", "estacion"], sort=False):
        for z0 in CFG["profundidades_estandar_m"]:
            if g.empty:
                continue
            d = np.abs(g["profundidad_m"].to_numpy(float) - float(z0))
            j = int(np.nanargmin(d))
            if d[j] > tol:
                continue
            r = g.iloc[j]
            rows.append({
                "id_perfil": pid,
                "estacion": est,
                "periodo": r["periodo"],
                "numero_estacion": r["numero_estacion"],
                "prof_objetivo_m": z0,
                "prof_muestreada_m": r["profundidad_m"],
                "frac_mSUW": r["frac_mSUW"],
                "frac_mAAIW": r["frac_mAAIW"],
                "frac_mNADW": r["frac_mNADW"],
                "mass_residual_pct": r["mass_residual_pct"],
            })
    return pd.DataFrame(rows)


# =============================================================================
# SECCIONES
# =============================================================================


def station_order_table(df: pd.DataFrame) -> pd.DataFrame:
    s = (
        df.groupby(["numero_estacion", "estacion"], as_index=False)
        .agg(latitud=("latitud", "median"), longitud=("longitud", "median"))
    )
    if CFG["orden_estaciones"]:
        order = {int(v): i for i, v in enumerate(CFG["orden_estaciones"])}
        s["_ord"] = s["numero_estacion"].map(lambda x: order.get(int(x), 9999) if np.isfinite(x) else 9999)
        s = s.sort_values(["_ord", "numero_estacion"]).drop(columns="_ord")
    else:
        s = s.sort_values("numero_estacion")
    return s


def profile_station_positions(df: pd.DataFrame) -> pd.DataFrame:
    s = (
        df.groupby("numero_estacion", as_index=False)
        .agg(latitud=("latitud", "median"), longitud=("longitud", "median"))
        .sort_values("numero_estacion")
    )
    if CFG["orden_estaciones"]:
        order = {int(v): i for i, v in enumerate(CFG["orden_estaciones"])}
        s["_ord"] = s["numero_estacion"].map(lambda x: order.get(int(x), 9999))
        s = s.sort_values("_ord").drop(columns="_ord")
    dist = [0.0]
    for i in range(1, len(s)):
        dist.append(
            dist[-1] + haversine_km(
                s.iloc[i - 1].latitud, s.iloc[i - 1].longitud,
                s.iloc[i].latitud, s.iloc[i].longitud,
            )
        )
    s["dist_km"] = dist
    return s


def section_matrix(df: pd.DataFrame, value_col: str, periodo: str, zgrid=None):
    x = df[df["periodo"].eq(periodo)].copy()
    if x.empty or value_col not in x:
        return None
    pos = profile_station_positions(x)
    if zgrid is None:
        zmax = float(np.nanmax(x["profundidad_m"]))
        zgrid = np.arange(0.0, zmax + CFG["seccion_dz_m"], CFG["seccion_dz_m"])
    M = np.full((len(zgrid), len(pos)), np.nan)
    labels = []
    for j, rr in pos.reset_index(drop=True).iterrows():
        st = rr.numero_estacion
        g = x[x["numero_estacion"].eq(st)].copy()
        # Si accidentalmente hay más de un perfil del mismo periodo/estación,
        # se promedia sobre bins de profundidad cercanos antes de interpolar.
        gg = (
            g[["profundidad_m", value_col]].dropna()
            .groupby("profundidad_m", as_index=False)[value_col].mean()
            .sort_values("profundidad_m")
        )
        if len(gg) >= 2:
            z = gg["profundidad_m"].to_numpy(float)
            v = gg[value_col].to_numpy(float)
            inside = (zgrid >= z.min()) & (zgrid <= z.max())
            M[inside, j] = np.interp(zgrid[inside], z, v)
        labels.append(f"E{int(st):02d}" if np.isfinite(st) else str(st))
    return pos["dist_km"].to_numpy(float), zgrid, M, labels, pos


# =============================================================================
# FIGURAS
# =============================================================================


def add_sigma0(ax, sa_data, ct_data):
    sa = np.asarray(sa_data, float)
    ct = np.asarray(ct_data, float)
    m = np.isfinite(sa) & np.isfinite(ct)
    if m.sum() < 10:
        return
    dx = max(0.05, np.ptp(sa[m]) * 0.05)
    dy = max(0.1, np.ptp(ct[m]) * 0.05)
    xx = np.linspace(sa[m].min() - dx, sa[m].max() + dx, 130)
    yy = np.linspace(ct[m].min() - dy, ct[m].max() + dy, 130)
    X, Y = np.meshgrid(xx, yy)
    sig = gsw.sigma0(X, Y)
    lev = MaxNLocator(nbins=8).tick_values(np.nanmin(sig), np.nanmax(sig))
    cs = ax.contour(X, Y, sig, levels=lev, linewidths=0.6, alpha=0.6)
    ax.clabel(cs, fmt="%.1f", fontsize=10)


def fig_mapa_estaciones(base, out):
    s = (
        base.groupby("numero_estacion", as_index=False)
        .agg(latitud=("latitud", "median"), longitud=("longitud", "median"), profundidad_max_m=("profundidad_m", "max"))
        .sort_values("numero_estacion")
    )
    fig, ax = plt.subplots(figsize=(8, 6), layout="constrained")
    sc = ax.scatter(s["longitud"], s["latitud"], c=s["profundidad_max_m"], s=70)
    for _, r in s.iterrows():
        ax.annotate(f"E{int(r.numero_estacion):02d}", (r.longitud, r.latitud), xytext=(4, 4), textcoords="offset points")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Profundidad máxima CTD [m]")
    ax.set_xlabel("Longitud [°]")
    ax.set_ylabel("Latitud [°]")
    ax.set_title("Estaciones CTD y profundidad máxima muestreada")
    ax.grid(alpha=0.2)
    savefig(fig, out / "Fig01_mapa_estaciones.png")


def fig_ts(base, o2col, out):
    x = filtro_hidrografico(base)
    # Profundidad como z negativa, según la convención gráfica solicitada.
    zplot = -np.abs(x["profundidad_m"].to_numpy(float))
    fig, ax = plt.subplots(figsize=(8.5, 6.8), layout="constrained")
    sc = ax.scatter(x["salinidad_absoluta_g_kg"], x["temperatura_conservativa_C"], c=zplot, s=7, alpha=0.75)
    add_sigma0(ax, x["salinidad_absoluta_g_kg"], x["temperatura_conservativa_C"])
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Coordenada vertical, z [m]")
    markers = {"mSUW": "o", "mAAIW": "s", "mNADW": "^"}
    for name in SWT.columns:
        ax.scatter(SWT.loc["SA_g_kg", name], SWT.loc["CT_C", name], marker=markers[name], s=90, label=name)
        ax.annotate(name, (SWT.loc["SA_g_kg", name], SWT.loc["CT_C", name]), xytext=(5, 5), textcoords="offset points")
    ax.set_xlabel(r"Salinidad absoluta, $S_A$ [g kg$^{-1}$]")
    ax.set_ylabel(r"Temperatura conservativa, $\Theta$ [°C]")
    ax.set_title("Diagrama Θ–SA de campaña · color por z · SWT de Torres et al. (2023)")
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, out / "Fig02_TS_profundidad_endmembers.png")

    if o2col is not None:
        xo = filtro_omp(base, o2col)
        if not xo.empty:
            fig, ax = plt.subplots(figsize=(8.5, 6.8), layout="constrained")
            sc = ax.scatter(xo["salinidad_absoluta_g_kg"], xo["temperatura_conservativa_C"], c=xo[o2col], s=7, alpha=0.75)
            add_sigma0(ax, xo["salinidad_absoluta_g_kg"], xo["temperatura_conservativa_C"])
            cb = fig.colorbar(sc, ax=ax)
            cb.set_label(r"Oxígeno disuelto [$\mu$mol kg$^{-1}$]")
            for name in SWT.columns:
                ax.scatter(SWT.loc["SA_g_kg", name], SWT.loc["CT_C", name], marker=markers[name], s=90, label=name)
                ax.annotate(name, (SWT.loc["SA_g_kg", name], SWT.loc["CT_C", name]), xytext=(5, 5), textcoords="offset points")
            ax.set_xlabel(r"Salinidad absoluta, $S_A$ [g kg$^{-1}$]")
            ax.set_ylabel(r"Temperatura conservativa, $\Theta$ [°C]")
            ax.set_title("Diagrama Θ–SA de campaña · color por O₂")
            ax.grid(alpha=0.2)
            ax.legend()
            savefig(fig, out / "Fig03_TS_oxigeno_endmembers.png")


def fig_capas_superficiales(nucleos, out):
    if nucleos.empty:
        return
    fig, ax = plt.subplots(figsize=(9, 6), layout="constrained")
    for periodo, g in nucleos.groupby("periodo"):
        g = g.sort_values("numero_estacion")
        ax.plot(g["numero_estacion"], -g["MLD_m"], marker="o", label=f"MLD {periodo}")
        ax.plot(g["numero_estacion"], -g["ILD_m"], marker="s", label=f"ILD {periodo}")
        ax.plot(g["numero_estacion"], -g["SUW_core_depth_m"], marker="^", label=f"Núcleo SUW {periodo}")
    ax.set_xlabel("Número de estación")
    ax.set_ylabel("Coordenada vertical, z [m]")
    ax.set_title("MLD, ILD y profundidad del máximo de salinidad (SUW)")
    ax.grid(alpha=0.2)
    ax.legend(ncol=2, fontsize=8)
    savefig(fig, out / "Fig04_MLD_ILD_SUW.png")


def fig_nucleos_intermedios(nucleos, out):
    if nucleos.empty:
        return
    fig, ax = plt.subplots(figsize=(9, 6), layout="constrained")
    for periodo, g in nucleos.groupby("periodo"):
        g = g.sort_values("numero_estacion")
        ax.plot(g["numero_estacion"], -g["SUW_core_depth_m"], marker="o", label=f"SUW {periodo}")
        good = g[~g["AAIW_core_truncado_fondo"].fillna(True)]
        ax.plot(good["numero_estacion"], -good["AAIW_core_depth_m"], marker="s", label=f"AAIW {periodo}")
    ax.set_xlabel("Número de estación")
    ax.set_ylabel("Coordenada vertical, z [m]")
    ax.set_title("Profundidad de los núcleos hidrográficos SUW y AAIW")
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, out / "Fig05_nucleos_SUW_AAIW.png")


def fig_secciones_hidrograficas(base, o2col, out):
    hydro = filtro_hidrografico(base)
    oxy = filtro_omp(base, o2col) if o2col else hydro.iloc[0:0].copy()
    variables = [
        ("salinidad_absoluta_g_kg", r"$S_A$ [g kg$^{-1}$]", hydro),
        ("temperatura_conservativa_C", r"$\Theta$ [°C]", hydro),
        ("sigma0_kg_m3", r"$\sigma_0$ [kg m$^{-3}$]", hydro),
    ]
    if o2col:
        variables.append((o2col, r"O$_2$ [$\mu$mol kg$^{-1}$]", oxy))

    for periodo in [p for p in ("Diurno", "Nocturno") if (hydro["periodo"] == p).any()]:
        fig, axs = plt.subplots(1, len(variables), figsize=(5.0 * len(variables), 6.5), layout="constrained", sharey=True)
        if len(variables) == 1:
            axs = [axs]
        for ax, (col, label, source) in zip(axs, variables):
            sec = section_matrix(source, col, periodo)
            if sec is None:
                ax.set_visible(False)
                continue
            xdist, zgrid, M, labels, _ = sec
            pc = ax.pcolormesh(xdist, -zgrid, M, shading="auto")
            cb = fig.colorbar(pc, ax=ax)
            cb.set_label(label)
            ax.set_xlabel("Distancia acumulada [km]")
            ax.set_title(label)
            ax.set_xticks(xdist)
            ax.set_xticklabels(labels, rotation=45)
            ax.grid(alpha=0.12)
        axs[0].set_ylabel("Coordenada vertical, z [m]")
        fig.suptitle(f"Secciones hidrográficas · {periodo}", fontsize=14)
        savefig(fig, out / f"Fig06_seccion_hidrografica_{periodo.lower()}.png")


def fig_secciones_omp(omp, out):
    if omp.empty:
        return
    q = omp[omp["omp_ok"]].copy()
    for periodo in [p for p in ("Diurno", "Nocturno") if (q["periodo"] == p).any()]:
        fig, axs = plt.subplots(1, 3, figsize=(16, 6.5), layout="constrained", sharey=True)
        for ax, mass in zip(axs, ["mSUW", "mAAIW", "mNADW"]):
            col = f"frac_{mass}"
            sec = section_matrix(q, col, periodo)
            if sec is None:
                ax.set_visible(False)
                continue
            xdist, zgrid, M, labels, _ = sec
            pc = ax.pcolormesh(xdist, -zgrid, 100 * M, shading="auto", vmin=0, vmax=100)
            cb = fig.colorbar(pc, ax=ax)
            cb.set_label(f"{mass} [%]")
            # Isolínea de 30 % para reproducir el umbral descriptivo del artículo.
            try:
                ax.contour(xdist, -zgrid, 100 * M, levels=[30], linewidths=1.0)
            except Exception:
                pass
            ax.set_xlabel("Distancia acumulada [km]")
            ax.set_title(mass)
            ax.set_xticks(xdist)
            ax.set_xticklabels(labels, rotation=45)
            ax.grid(alpha=0.12)
        axs[0].set_ylabel("Coordenada vertical, z [m]")
        fig.suptitle(f"Fracciones OMP mSUW–mAAIW–mNADW · {periodo}", fontsize=14)
        savefig(fig, out / f"Fig07_OMP_{periodo.lower()}.png")


def fig_omp_perfil_medio(omp, out):
    if omp.empty:
        return
    q = omp[omp["omp_ok"]].copy()
    if q.empty:
        return
    dz = 25.0
    q["zbin"] = np.round(q["profundidad_m"] / dz) * dz
    m = q.groupby("zbin")[["frac_mSUW", "frac_mAAIW", "frac_mNADW"]].mean().reset_index()
    fig, ax = plt.subplots(figsize=(7, 7), layout="constrained")
    ax.plot(100 * m["frac_mSUW"], -m["zbin"], label="mSUW")
    ax.plot(100 * m["frac_mAAIW"], -m["zbin"], label="mAAIW")
    ax.plot(100 * m["frac_mNADW"], -m["zbin"], label="mNADW")
    ax.axvline(30, linestyle="--", linewidth=1)
    ax.set_xlabel("Fracción media en los perfiles muestreados [%]")
    ax.set_ylabel("Coordenada vertical, z [m]")
    ax.set_title("Estructura vertical media de fracciones OMP\n(promedio muestral, NO fracción volumétrica del Caribe)")
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, out / "Fig08_OMP_perfil_medio_muestral.png")


def fig_residual_omp(omp, out):
    if omp.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 6), layout="constrained")
    ax.scatter(omp["mass_residual_pct"], -omp["profundidad_m"], s=8, alpha=0.6)
    ax.axvline(CFG["omp_mass_residual_max_pct"], linestyle="--", linewidth=1)
    ax.set_xlabel("Residual de conservación de masa [%]")
    ax.set_ylabel("Coordenada vertical, z [m]")
    ax.set_title("Diagnóstico del ajuste OMP")
    ax.grid(alpha=0.2)
    savefig(fig, out / "Fig09_residual_masa_OMP.png")


# =============================================================================
# COMPARACIÓN DIURNO / NOCTURNO
# =============================================================================


def comparar_dia_noche(nucleos: pd.DataFrame, omp_std: pd.DataFrame) -> pd.DataFrame:
    if nucleos.empty:
        return pd.DataFrame(), pd.DataFrame()
    rows = []
    for st, g in nucleos.groupby("numero_estacion"):
        d = g[g["periodo"].eq("Diurno")]
        n = g[g["periodo"].eq("Nocturno")]
        if len(d) != 1 or len(n) != 1:
            continue
        rd, rn = d.iloc[0], n.iloc[0]
        row = {
            "numero_estacion": st,
            "MLD_dia_minus_noche_m": rd.get("MLD_m", np.nan) - rn.get("MLD_m", np.nan),
            "ILD_dia_minus_noche_m": rd.get("ILD_m", np.nan) - rn.get("ILD_m", np.nan),
            "BLT_dia_minus_noche_m": rd.get("BLT_m", np.nan) - rn.get("BLT_m", np.nan),
            "SUW_depth_dia_minus_noche_m": rd.get("SUW_core_depth_m", np.nan) - rn.get("SUW_core_depth_m", np.nan),
            "SUW_SA_dia_minus_noche": rd.get("SUW_core_SA_g_kg", np.nan) - rn.get("SUW_core_SA_g_kg", np.nan),
            "AAIW_depth_dia_minus_noche_m": rd.get("AAIW_core_depth_m", np.nan) - rn.get("AAIW_core_depth_m", np.nan),
            "AAIW_SA_dia_minus_noche": rd.get("AAIW_core_SA_g_kg", np.nan) - rn.get("AAIW_core_SA_g_kg", np.nan),
        }
        rows.append(row)
    comp = pd.DataFrame(rows)

    # Añade diferencias OMP en profundidades estándar si hay pares completos.
    if not omp_std.empty and not comp.empty:
        extra = []
        for st, g in omp_std.groupby("numero_estacion"):
            for z0, gz in g.groupby("prof_objetivo_m"):
                d = gz[gz["periodo"].eq("Diurno")]
                n = gz[gz["periodo"].eq("Nocturno")]
                if len(d) == 1 and len(n) == 1:
                    extra.append({
                        "numero_estacion": st,
                        "prof_objetivo_m": z0,
                        "d_mSUW_dia_minus_noche": float(d.iloc[0]["frac_mSUW"] - n.iloc[0]["frac_mSUW"]),
                        "d_mAAIW_dia_minus_noche": float(d.iloc[0]["frac_mAAIW"] - n.iloc[0]["frac_mAAIW"]),
                        "d_mNADW_dia_minus_noche": float(d.iloc[0]["frac_mNADW"] - n.iloc[0]["frac_mNADW"]),
                    })
        # Se guarda aparte por naturaleza long-format; no se fusiona con comp.
        if extra:
            return comp, pd.DataFrame(extra)
    return comp, pd.DataFrame()


# =============================================================================
# EXPORTACIÓN
# =============================================================================


def write_methodology(out: Path, path_input: Path, o2col: str | None, omp: pd.DataFrame):
    if omp.empty:
        omp_note = "No se ejecutó OMP de tres end-members porque no hubo O2 utilizable según la configuración."
    else:
        ok = omp[omp["omp_ok"]]
        omp_note = (
            f"OMP ejecutado en {len(omp)} bins; {len(ok)} cumplieron residual de masa <= "
            f"{CFG['omp_mass_residual_max_pct']:.0f} %. Residual medio de masa en bins aceptados: "
            f"{ok['mass_residual_pct'].mean():.2f} %."
            if len(ok) else "OMP ejecutado, pero ningún bin pasó el criterio de residual de masa."
        )

    text = f"""# Metodología — análisis de masas de agua

Generado: {datetime.now().isoformat(timespec='seconds')}
Script: analisis_masas_agua_caribe.py v{VERSION}
Entrada: {path_input}

## Marco conceptual

El análisis sigue como referencia principal a Torres et al. (2023), Ocean Dynamics 73:39–57,
DOI 10.1007/s10236-022-01529-5.

- MLD: Δsigma0 = 0.03 kg m-3 respecto a 10 m.
- ILD: |ΔΘ| = 0.2 °C respecto a 10 m.
- Núcleo SUW: máximo de SA por debajo de 70 m.
- Núcleo AAIW: mínimo de SA por debajo del máximo de salinidad; en este script se exige además
  profundidad >= {CFG['aaiw_min_depth_m']:.0f} m como control operacional para no confundir el mínimo
  con aguas centrales superiores.
- OMP mSUW–mAAIW–mNADW: aplicado únicamente por debajo del máximo de salinidad.
- Fracción >= {100*CFG['fraccion_significativa']:.0f} %: contribución significativa descriptiva.
- Criterio de residual de masa: <= {CFG['omp_mass_residual_max_pct']:.0f} %.

## Source water types y pesos

{SWT.to_string()}

Pesos: {OMP_WEIGHTS}

El sistema OMP se normaliza por propiedad usando la media y desviación estándar de los SWT,
se pondera con los pesos publicados y se resuelve mediante mínimos cuadrados no negativos.
El residual de masa se conserva como diagnóstico. Las fracciones normalizadas se calculan
solo para soluciones que cumplen el criterio de residual.

Columna O2 utilizada: {o2col or 'NO DISPONIBLE'}

{omp_note}

## QC

Los bins QARTOD 4 (FAIL) y 9 (MISSING) no se usan. Los bins 3 (SUSPECT/HIGH INTEREST) se
conservan porque pueden representar gradientes oceanográficos reales. Para OMP, por defecto
se exige estado dinámico de O2 = AUTO_OK por perfil.

## Interpretación

Las fracciones OMP son locales al dominio y profundidades muestreadas. El promedio muestral
NO debe interpretarse como porcentaje volumétrico del Caribe. Torres et al. (2023) obtuvieron
fracciones volumétricas basin-wide usando WOA18 y una geometría espacial que este crucero no
reproduce.

La comparación Diurno–Nocturno es descriptiva. No debe presentarse como ciclo diurno sin
considerar simultaneidad, fechas de muestreo, separación temporal, advección y variabilidad
meso/submesoescala.
"""
    (out / "METODOLOGIA_MASAS_AGUA.md").write_text(text, encoding="utf-8")


def exportar_resultados(out, nucleos, omp, omp_res, omp_std, excl, comp_dn, comp_omp_dn):
    out.mkdir(parents=True, exist_ok=True)

    nucleos.to_csv(out / "nucleos_masas_agua_por_perfil.csv", index=False, encoding="utf-8-sig")
    if not omp.empty:
        omp.to_csv(out / "OMP_por_bin.csv", index=False, encoding="utf-8-sig")
    if not omp_res.empty:
        omp_res.to_csv(out / "OMP_resumen_por_perfil.csv", index=False, encoding="utf-8-sig")
    if not omp_std.empty:
        omp_std.to_csv(out / "OMP_profundidades_estandar.csv", index=False, encoding="utf-8-sig")
    if not excl.empty:
        excl.to_csv(out / "OMP_perfiles_excluidos.csv", index=False, encoding="utf-8-sig")
    if not comp_dn.empty:
        comp_dn.to_csv(out / "comparacion_diurno_nocturno_nucleos.csv", index=False, encoding="utf-8-sig")
    if not comp_omp_dn.empty:
        comp_omp_dn.to_csv(out / "comparacion_diurno_nocturno_OMP.csv", index=False, encoding="utf-8-sig")

    swt_out = SWT.copy()
    swt_out["OMP_weight"] = [OMP_WEIGHTS["CT_C"], OMP_WEIGHTS["SA_g_kg"], OMP_WEIGHTS["O2_umol_kg"]]

    xlsx = out / "RESULTADOS_MASAS_AGUA.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
        nucleos.to_excel(writer, sheet_name="Nucleos", index=False)
        if not omp.empty:
            omp.to_excel(writer, sheet_name="OMP_bins", index=False)
        if not omp_res.empty:
            omp_res.to_excel(writer, sheet_name="OMP_resumen", index=False)
        if not omp_std.empty:
            omp_std.to_excel(writer, sheet_name="OMP_prof_estandar", index=False)
        if not excl.empty:
            excl.to_excel(writer, sheet_name="OMP_excluidos", index=False)
        if not comp_dn.empty:
            comp_dn.to_excel(writer, sheet_name="Dia_Noche_nucleos", index=False)
        if not comp_omp_dn.empty:
            comp_omp_dn.to_excel(writer, sheet_name="Dia_Noche_OMP", index=False)
        swt_out.to_excel(writer, sheet_name="SWT_Torres2023")
        pd.DataFrame([CFG]).astype(str).to_excel(writer, sheet_name="Configuracion", index=False)

    return xlsx




# =============================================================================
# V1.1 — REVISIÓN PARA MANUSCRITO
# =============================================================================


def savefig(fig, path: Path):
    """Guarda PNG de alta resolución y, opcionalmente, PDF vectorial."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=CFG["dpi"], bbox_inches="tight", facecolor="white")
    if CFG.get("guardar_pdf", True):
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    if CFG["mostrar_figuras"]:
        plt.show()
    else:
        plt.close(fig)


def _filter_qc_class(df: pd.DataFrame, qcol: str, flag: int) -> pd.DataFrame:
    if qcol in df:
        return df[df[qcol].eq(flag)].copy()
    # Si no hay QARTOD, todo se considera PASS a efectos gráficos.
    return df.copy() if flag == CFG["qartod_pass"] else df.iloc[0:0].copy()


def filtro_hidrografico_pass(df: pd.DataFrame) -> pd.DataFrame:
    x = _filter_qc_class(df, "qartod_resumen_ctd", CFG["qartod_pass"])
    req = ["profundidad_m", "salinidad_absoluta_g_kg", "temperatura_conservativa_C", "sigma0_kg_m3"]
    return x[x[req].notna().all(axis=1)].copy()


def filtro_hidrografico_suspect(df: pd.DataFrame) -> pd.DataFrame:
    x = _filter_qc_class(df, "qartod_resumen_ctd", CFG["qartod_suspect"])
    req = ["profundidad_m", "salinidad_absoluta_g_kg", "temperatura_conservativa_C", "sigma0_kg_m3"]
    return x[x[req].notna().all(axis=1)].copy()


def filtro_o2_pass(df: pd.DataFrame, o2col: str | None) -> pd.DataFrame:
    if o2col is None:
        return df.iloc[0:0].copy()
    x = filtro_hidrografico_pass(df)
    if "qartod_oxigeno" in x:
        x = x[x["qartod_oxigeno"].eq(CFG["qartod_pass"])].copy()
    x = x[x[o2col].notna()].copy()
    return x


def filtro_o2_suspect(df: pd.DataFrame, o2col: str | None) -> pd.DataFrame:
    if o2col is None:
        return df.iloc[0:0].copy()
    x = filtro_hidrografico(df)
    if "qartod_oxigeno" in x:
        x = x[x["qartod_oxigeno"].eq(CFG["qartod_suspect"])].copy()
    else:
        return x.iloc[0:0].copy()
    return x[x[o2col].notna()].copy()


def _robust_limits(values, qlo=0.005, qhi=0.995, pad=0.04):
    a = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(float)
    if len(a) == 0:
        return None
    lo, hi = np.quantile(a, [qlo, qhi]) if len(a) >= 20 else (a.min(), a.max())
    span = max(float(hi-lo), 1e-8)
    return float(lo-pad*span), float(hi+pad*span)


def add_sigma0_fixed(ax, xlim, ylim):
    xx = np.linspace(xlim[0], xlim[1], 160)
    yy = np.linspace(ylim[0], ylim[1], 160)
    X, Y = np.meshgrid(xx, yy)
    sig = gsw.sigma0(X, Y)
    lev = MaxNLocator(nbins=8).tick_values(float(np.nanmin(sig)), float(np.nanmax(sig)))
    cs = ax.contour(X, Y, sig, levels=lev, linewidths=0.65, alpha=0.55)
    ax.clabel(cs, fmt="%.1f", fontsize=10)


@dataclass
class OMPTSSolution:
    frac_mSUW: float = np.nan
    frac_mAAIW: float = np.nan
    frac_mNADW: float = np.nan
    mass_residual_pct: float = np.nan
    resid_CT_C: float = np.nan
    resid_SA_g_kg: float = np.nan
    weighted_rms: float = np.nan
    omp_ok: bool = False


def omp_solve_ts(ct_obs, sa_obs) -> OMPTSSolution:
    """OMP de sensibilidad usando solo Θ, SA y conservación de masa.

    NO sustituye al OMP completo. Se usa para cuantificar cuánto dependen las
    fracciones mAAIW/mNADW del O2 absoluto del SBE43.
    """
    if not all(np.isfinite([ct_obs, sa_obs])):
        return OMPTSSolution()
    means = SWT.loc[["CT_C", "SA_g_kg"]].mean(axis=1)
    stds = SWT.loc[["CT_C", "SA_g_kg"]].std(axis=1, ddof=1)
    G = np.vstack([
        (SWT.loc["CT_C"].to_numpy(float) - means["CT_C"]) / stds["CT_C"],
        (SWT.loc["SA_g_kg"].to_numpy(float) - means["SA_g_kg"]) / stds["SA_g_kg"],
        np.ones(3),
    ])
    d = np.array([
        (ct_obs - means["CT_C"]) / stds["CT_C"],
        (sa_obs - means["SA_g_kg"]) / stds["SA_g_kg"],
        1.0,
    ])
    w = np.array([OMP_WEIGHTS["CT_C"], OMP_WEIGHTS["SA_g_kg"], OMP_WEIGHTS["mass"]], float)
    sol = lsq_linear(G*w[:, None], d*w, bounds=(0.0, np.inf), method="trf", lsmr_tol="auto")
    if not sol.success:
        return OMPTSSolution()
    f = np.asarray(sol.x, float)
    s = float(f.sum())
    if s <= 0:
        return OMPTSSolution()
    fn = f/s
    mass = 100.0*abs(s-1.0)
    ct_rec = float(np.dot(SWT.loc["CT_C"].to_numpy(float), fn))
    sa_rec = float(np.dot(SWT.loc["SA_g_kg"].to_numpy(float), fn))
    wrms = float(np.sqrt(np.mean((G@f-d)**2)))
    ok = bool(mass <= CFG["omp_mass_residual_max_pct"])
    return OMPTSSolution(float(fn[0]), float(fn[1]), float(fn[2]), mass,
                         float(ct_obs-ct_rec), float(sa_obs-sa_rec), wrms, ok)


def _nnls_fixed_matrix(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """NNLS vectorizado para una matriz A fija y muchos vectores observados.

    Con tres end-members se evalúan las 7 combinaciones activas posibles. Es
    algebraicamente equivalente a buscar la mejor solución no negativa sobre
    las caras del ortante, pero evita miles de llamadas a un optimizador.
    """
    from itertools import combinations
    n=B.shape[0]
    best=np.zeros((n,A.shape[1]),float)
    best_cost=np.full(n,np.inf)
    for k in range(1,A.shape[1]+1):
        for idx in combinations(range(A.shape[1]),k):
            idx=list(idx)
            P=np.linalg.pinv(A[:,idx])              # k x m
            X=B @ P.T                               # n x k
            valid=np.all(X >= -1e-10,axis=1)
            X=np.where(X<0,0,X)
            pred=X @ A[:,idx].T
            cost=np.sum((pred-B)**2,axis=1)
            upd=valid & (cost < best_cost)
            if np.any(upd):
                best[upd]=0.0
                best[np.ix_(upd,idx)]=X[upd]
                best_cost[upd]=cost[upd]
    return best


def ejecutar_omp_ts(base: pd.DataFrame, nucleos: pd.DataFrame, profile_col: str) -> pd.DataFrame:
    """Sensibilidad Θ-SA vectorizada para todos los perfiles hidrográficos."""
    suw_depth = dict(zip(nucleos["estacion"].astype(str), nucleos["SUW_core_depth_m"]))
    pieces=[]
    for pid,g in base.groupby(profile_col,sort=False):
        est=str(g["estacion"].iloc[0])
        z0=pd.to_numeric(pd.Series([suw_depth.get(est,np.nan)]),errors="coerce").iloc[0]
        if not np.isfinite(z0): continue
        gg=filtro_hidrografico(g)
        gg=gg[gg["profundidad_m"] >= float(z0)].copy()
        if gg.empty: continue
        gg["id_perfil"]=pid
        pieces.append(gg)
    if not pieces: return pd.DataFrame()
    q=pd.concat(pieces,ignore_index=True,sort=False)

    means=SWT.loc[["CT_C","SA_g_kg"]].mean(axis=1)
    stds=SWT.loc[["CT_C","SA_g_kg"]].std(axis=1,ddof=1)
    G=np.vstack([
        (SWT.loc["CT_C"].to_numpy(float)-means["CT_C"])/stds["CT_C"],
        (SWT.loc["SA_g_kg"].to_numpy(float)-means["SA_g_kg"])/stds["SA_g_kg"],
        np.ones(3),
    ])
    w=np.array([OMP_WEIGHTS["CT_C"],OMP_WEIGHTS["SA_g_kg"],OMP_WEIGHTS["mass"]],float)
    A=G*w[:,None]
    B=np.column_stack([
        (q["temperatura_conservativa_C"].to_numpy(float)-means["CT_C"])/stds["CT_C"],
        (q["salinidad_absoluta_g_kg"].to_numpy(float)-means["SA_g_kg"])/stds["SA_g_kg"],
        np.ones(len(q)),
    ])*w
    raw=_nnls_fixed_matrix(A,B)
    ss=raw.sum(axis=1)
    good=ss>0
    fn=np.full_like(raw,np.nan);fn[good]=raw[good]/ss[good,None]
    mass=100*np.abs(ss-1.0)
    pred=raw @ A.T
    wrms=np.sqrt(np.mean((pred-B)**2,axis=1))
    ct_rec=np.nansum(fn*SWT.loc["CT_C"].to_numpy(float),axis=1)
    sa_rec=np.nansum(fn*SWT.loc["SA_g_kg"].to_numpy(float),axis=1)
    ok=good & (mass <= CFG["omp_mass_residual_max_pct"])

    out=pd.DataFrame({
        "id_perfil":q["id_perfil"].to_numpy(),"estacion":q["estacion"].astype(str).to_numpy(),
        "periodo":q["periodo"].to_numpy(),"numero_estacion":q["numero_estacion"].to_numpy(),
        "latitud":q.get("latitud",pd.Series(np.nan,index=q.index)).to_numpy(),
        "longitud":q.get("longitud",pd.Series(np.nan,index=q.index)).to_numpy(),
        "profundidad_m":q["profundidad_m"].to_numpy(float),
        "CT_C":q["temperatura_conservativa_C"].to_numpy(float),
        "SA_g_kg":q["salinidad_absoluta_g_kg"].to_numpy(float),
        "sigma0_kg_m3":q.get("sigma0_kg_m3",pd.Series(np.nan,index=q.index)).to_numpy(),
        "frac_mSUW":fn[:,0],"frac_mAAIW":fn[:,1],"frac_mNADW":fn[:,2],
        "mass_residual_pct":mass,"resid_CT_C":q["temperatura_conservativa_C"].to_numpy(float)-ct_rec,
        "resid_SA_g_kg":q["salinidad_absoluta_g_kg"].to_numpy(float)-sa_rec,
        "weighted_rms":wrms,"omp_ok":ok,
    })
    return out


def augmentar_diagnostico_omp(omp: pd.DataFrame, omp_ts: pd.DataFrame) -> tuple[pd.DataFrame,pd.DataFrame]:
    """Añade residuales con fracciones normalizadas y sensibilidad al O2."""
    if omp.empty:
        return omp, pd.DataFrame()
    x=omp.copy()
    F=x[["frac_mSUW","frac_mAAIW","frac_mNADW"]].to_numpy(float)
    rec=F @ SWT.loc[["CT_C","SA_g_kg","O2_umol_kg"]].to_numpy(float).T
    x["resid_CT_norm_C"] = x["CT_C"].to_numpy(float)-rec[:,0]
    x["resid_SA_norm_g_kg"] = x["SA_g_kg"].to_numpy(float)-rec[:,1]
    x["resid_O2_norm_umol_kg"] = x["O2_umol_kg"].to_numpy(float)-rec[:,2]

    # Cruce con OMP Θ-SA por perfil y profundidad (nearest merge con tolerancia 1.5 m).
    merged=[]
    if not omp_ts.empty:
        for pid, g in x.groupby("id_perfil", sort=False):
            h=omp_ts[omp_ts["id_perfil"].eq(pid)].sort_values("profundidad_m")
            gg=g.sort_values("profundidad_m")
            if h.empty:
                merged.append(gg)
                continue
            cols=["profundidad_m","frac_mSUW","frac_mAAIW","frac_mNADW"]
            hh=h[cols].rename(columns={
                "frac_mSUW":"fracTS_mSUW","frac_mAAIW":"fracTS_mAAIW","frac_mNADW":"fracTS_mNADW"})
            mm=pd.merge_asof(gg,hh,on="profundidad_m",direction="nearest",tolerance=1.5)
            merged.append(mm)
        x=pd.concat(merged,ignore_index=True,sort=False)
        if {"fracTS_mSUW","fracTS_mAAIW","fracTS_mNADW"}.issubset(x):
            FT=x[["fracTS_mSUW","fracTS_mAAIW","fracTS_mNADW"]].to_numpy(float)
            pred_o2=FT @ SWT.loc["O2_umol_kg"].to_numpy(float)
            x["O2_predicho_desde_TS_umol_kg"]=pred_o2
            x["resid_O2_vs_TS_umol_kg"]=x["O2_umol_kg"].to_numpy(float)-pred_o2
            for m in ("mSUW","mAAIW","mNADW"):
                x[f"delta_{m}_full_minus_TS"] = x[f"frac_{m}"]-x[f"fracTS_{m}"]

    rows=[]
    for c in ["resid_CT_norm_C","resid_SA_norm_g_kg","resid_O2_norm_umol_kg","resid_O2_vs_TS_umol_kg"]:
        if c not in x:
            continue
        v=pd.to_numeric(x[c],errors="coerce").dropna().to_numpy(float)
        if len(v):
            med=float(np.median(v)); mad=float(np.median(np.abs(v-med)))
            same=max(float(np.mean(v>0)),float(np.mean(v<0)))
            rows.append({"variable":c,"n":len(v),"media":float(np.mean(v)),"mediana":med,
                         "MAD":mad,"p05":float(np.quantile(v,.05)),"p95":float(np.quantile(v,.95)),
                         "fraccion_mismo_signo":same})
    return x,pd.DataFrame(rows)


def _interp_omp_standard(omp: pd.DataFrame, cols: list[str], zlist=None, tol=None) -> pd.DataFrame:
    if omp.empty:
        return pd.DataFrame()
    zlist = CFG["profundidades_estandar_m"] if zlist is None else zlist
    tol = CFG["tolerancia_prof_estandar_m"] if tol is None else tol
    rows=[]
    for (pid,est),g in omp.groupby(["id_perfil","estacion"],sort=False):
        for z0 in zlist:
            d=np.abs(g["profundidad_m"].to_numpy(float)-float(z0))
            if len(d)==0 or not np.isfinite(d).any(): continue
            j=int(np.nanargmin(d))
            if d[j] > tol: continue
            r=g.iloc[j]
            row={"id_perfil":pid,"estacion":est,"periodo":r.get("periodo",""),
                 "numero_estacion":r.get("numero_estacion",np.nan),"prof_objetivo_m":z0,
                 "prof_muestreada_m":r["profundidad_m"]}
            for c in cols:
                row[c]=r.get(c,np.nan)
            rows.append(row)
    return pd.DataFrame(rows)


def resumen_hidrografico_masas(nucleos: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    # SUW todos los perfiles.
    for mass,prefix,mask in [
        ("SUW","SUW_core",np.ones(len(nucleos),bool)),
        ("AAIW","AAIW_core",~nucleos["AAIW_core_truncado_fondo"].fillna(True).to_numpy(bool)),
    ]:
        g=nucleos.loc[mask]
        row={"masa":mass,"n":len(g)}
        mapping={"profundidad_m":f"{prefix}_depth_m","SA_g_kg":f"{prefix}_SA_g_kg",
                 "CT_C":f"{prefix}_CT_C","sigma0_kg_m3":f"{prefix}_sigma0",
                 "O2_umol_kg":f"{prefix}_O2_umol_kg"}
        for outc,inc in mapping.items():
            v=pd.to_numeric(g.get(inc),errors="coerce").dropna()
            row[outc+"_media"]=float(v.mean()) if len(v) else np.nan
            row[outc+"_sd"]=float(v.std(ddof=1)) if len(v)>1 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _station_summary(base):
    return (base.groupby("numero_estacion",as_index=False)
            .agg(latitud=("latitud","median"),longitud=("longitud","median"),
                 profundidad_max_m=("profundidad_m","max"))
            .sort_values("numero_estacion"))




def fig_mapa_estaciones_vf(base, out):
    """
    Mapa regional + detalle con cartografía REAL.

    Prioridad:
    1) Cartopy + Natural Earth.
    2) Basemap.
    3) Si ninguna librería cartográfica está disponible, se genera un
       gráfico de estaciones SIN inventar una línea de costa.

    No se usa ninguna costa esquemática artificial.
    """
    s = _station_summary(base)
    if s.empty:
        return

    lon = pd.to_numeric(s["longitud"], errors="coerce")
    lat = pd.to_numeric(s["latitud"], errors="coerce")

    lon_min, lon_max = float(lon.min()), float(lon.max())
    lat_min, lat_max = float(lat.min()), float(lat.max())
    lon_c = 0.5 * (lon_min + lon_max)
    lat_c = 0.5 * (lat_min + lat_max)

    # ------------------------------------------------------------------
    # Extensiones
    # ------------------------------------------------------------------
    # Detalle: mantiene la malla visible, pero deja bastante espacio
    # geográfico para que aparezca la costa real si las estaciones son
    # oceánicas.
    half_lon_detail = max(1.35, 0.5 * (lon_max - lon_min) + 0.55)
    half_lat_detail = max(0.95, 0.5 * (lat_max - lat_min) + 0.45)

    detail_extent = [
        lon_c - half_lon_detail,
        lon_c + half_lon_detail,
        lat_c - half_lat_detail,
        lat_c + half_lat_detail,
    ]

    # Contexto regional más amplio para que el lector se ubique en el Caribe
    # colombiano y no sólo vea la malla CTD.
    half_lon_reg = max(2.40, half_lon_detail * 1.55)
    half_lat_reg = max(1.60, half_lat_detail * 1.45)

    regional_extent = [
        lon_c - half_lon_reg,
        lon_c + half_lon_reg,
        lat_c - half_lat_reg,
        lat_c + half_lat_reg,
    ]

    # ------------------------------------------------------------------
    # Opción 1: CARTOPY (preferida)
    # ------------------------------------------------------------------
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        proj = ccrs.PlateCarree()

        fig = plt.figure(figsize=(14.5, 6.6), layout="constrained")
        ax0 = fig.add_subplot(1, 2, 1, projection=proj)
        ax1 = fig.add_subplot(1, 2, 2, projection=proj)
        axs = [ax0, ax1]

        def _decorate_cartopy(ax, extent):
            ax.set_extent(extent, crs=proj)
            ax.add_feature(
                cfeature.NaturalEarthFeature(
                    category="physical",
                    name="land",
                    scale="10m",
                    facecolor="0.86",
                    edgecolor="none"
                ),
                zorder=0
            )
            ax.add_feature(
                cfeature.NaturalEarthFeature(
                    category="physical",
                    name="ocean",
                    scale="10m",
                    facecolor="#eaf2fb",
                    edgecolor="none"
                ),
                zorder=-1
            )
            ax.coastlines(resolution="10m", linewidth=0.9, color="0.25", zorder=2)
            ax.add_feature(cfeature.BORDERS.with_scale("10m"),
                           linewidth=0.45, edgecolor="0.45", zorder=2)
            gl = ax.gridlines(
                crs=proj, draw_labels=True,
                linewidth=0.45, color="0.65",
                alpha=0.55, linestyle="--"
            )
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {"size": 10}
            gl.ylabel_style = {"size": 10}

        _decorate_cartopy(ax0, regional_extent)
        _decorate_cartopy(ax1, detail_extent)

        sc0 = ax0.scatter(
            s.longitud, s.latitud,
            c=s.profundidad_max_m,
            cmap="viridis",
            s=72,
            edgecolor="k",
            linewidth=0.45,
            transform=proj,
            zorder=5
        )

        for _, r in s.iterrows():
            ax0.annotate(
                f"E{int(r.numero_estacion):02d}",
                (r.longitud, r.latitud),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=10,
                transform=proj,
                zorder=6
            )

        sc = ax1.scatter(
            s.longitud, s.latitud,
            c=s.profundidad_max_m,
            cmap="viridis",
            s=100,
            edgecolor="k",
            linewidth=0.50,
            transform=proj,
            zorder=5
        )

        for _, r in s.iterrows():
            ax1.annotate(
                f"E{int(r.numero_estacion):02d}",
                (r.longitud, r.latitud),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=11,
                transform=proj,
                zorder=6
            )

        for name, stations in CFG["transectos"].items():
            gg = s[s.numero_estacion.isin(stations)].copy()
            order = {v: i for i, v in enumerate(stations)}
            gg["_o"] = gg.numero_estacion.map(order)
            gg = gg.sort_values("_o")

            if len(gg) >= 2:
                ax1.plot(
                    gg.longitud, gg.latitud,
                    "--",
                    lw=1.1,
                    alpha=0.75,
                    transform=proj,
                    zorder=4,
                    label=f"Grupo {name}"
                )

        ax0.set_title("Contexto regional")
        ax1.set_title("Detalle de estaciones")
        ax1.legend(fontsize=10, loc="best")

        cb = fig.colorbar(sc, ax=axs, pad=0.025, shrink=0.88)
        cb.set_label("Profundidad máxima muestreada [m]")

        fig.suptitle(
            "Estaciones CTD en el Caribe colombiano",
            fontsize=CFG["font_suptitle"],
            weight="bold"
        )

        savefig(
            fig,
            out / "Fig01_V3_4_mapa_estaciones_cartografia_real.png"
        )
        return

    except Exception as exc_cartopy:
        warnings.warn(
            f"No se pudo usar Cartopy/Natural Earth ({exc_cartopy}). "
            "Se intentará Basemap."
        )

    # ------------------------------------------------------------------
    # Opción 2: BASEMAP
    # ------------------------------------------------------------------
    try:
        from mpl_toolkits.basemap import Basemap

        fig, axs = plt.subplots(
            1, 2,
            figsize=(14.5, 6.6),
            layout="constrained"
        )
        ax0, ax1 = axs

        def _decorate_basemap(ax, extent):
            xmin, xmax, ymin, ymax = extent
            m = Basemap(
                projection="cyl",
                llcrnrlon=xmin,
                llcrnrlat=ymin,
                urcrnrlon=xmax,
                urcrnrlat=ymax,
                resolution="i",
                ax=ax
            )
            m.drawmapboundary(
                fill_color="#eaf2fb",
                linewidth=0.8
            )
            m.fillcontinents(
                color="0.86",
                lake_color="#eaf2fb",
                zorder=0
            )
            m.drawcoastlines(
                linewidth=0.9,
                color="0.25",
                zorder=2
            )
            m.drawcountries(
                linewidth=0.45,
                color="0.45",
                zorder=2
            )
            m.drawparallels(
                np.arange(np.floor(ymin), np.ceil(ymax) + 0.1, 0.5),
                labels=[1, 0, 0, 0],
                linewidth=0.35,
                color="0.7",
                fontsize=9
            )
            m.drawmeridians(
                np.arange(np.floor(xmin), np.ceil(xmax) + 0.1, 0.5),
                labels=[0, 0, 0, 1],
                linewidth=0.35,
                color="0.7",
                fontsize=9
            )
            return m

        _decorate_basemap(ax0, regional_extent)
        _decorate_basemap(ax1, detail_extent)

        sc0 = ax0.scatter(
            s.longitud, s.latitud,
            c=s.profundidad_max_m,
            cmap="viridis",
            s=72,
            edgecolor="k",
            linewidth=0.45,
            zorder=5
        )

        for _, r in s.iterrows():
            ax0.annotate(
                f"E{int(r.numero_estacion):02d}",
                (r.longitud, r.latitud),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=10,
                zorder=6
            )

        sc = ax1.scatter(
            s.longitud, s.latitud,
            c=s.profundidad_max_m,
            cmap="viridis",
            s=100,
            edgecolor="k",
            linewidth=0.50,
            zorder=5
        )

        for _, r in s.iterrows():
            ax1.annotate(
                f"E{int(r.numero_estacion):02d}",
                (r.longitud, r.latitud),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=11,
                zorder=6
            )

        for name, stations in CFG["transectos"].items():
            gg = s[s.numero_estacion.isin(stations)].copy()
            order = {v: i for i, v in enumerate(stations)}
            gg["_o"] = gg.numero_estacion.map(order)
            gg = gg.sort_values("_o")

            if len(gg) >= 2:
                ax1.plot(
                    gg.longitud, gg.latitud,
                    "--",
                    lw=1.1,
                    alpha=0.75,
                    zorder=4,
                    label=f"Grupo {name}"
                )

        ax0.set_title("Contexto regional")
        ax1.set_title("Detalle de estaciones")
        ax1.legend(fontsize=10, loc="best")

        cb = fig.colorbar(sc, ax=axs, pad=0.025, shrink=0.88)
        cb.set_label("Profundidad máxima muestreada [m]")

        fig.suptitle(
            "Estaciones CTD en el Caribe colombiano",
            fontsize=CFG["font_suptitle"],
            weight="bold"
        )

        savefig(
            fig,
            out / "Fig01_V3_4_mapa_estaciones_cartografia_real.png"
        )
        return

    except Exception as exc_basemap:
        warnings.warn(
            f"No se pudo usar Basemap ({exc_basemap}). "
            "La figura se generará sin línea de costa artificial."
        )

    # ------------------------------------------------------------------
    # Opción 3: sólo estaciones, sin inventar geografía
    # ------------------------------------------------------------------
    fig, axs = plt.subplots(
        1, 2,
        figsize=(14.5, 6.6),
        layout="constrained"
    )
    ax0, ax1 = axs

    sc0 = ax0.scatter(
        s.longitud, s.latitud,
        c=s.profundidad_max_m,
        cmap="viridis",
        s=72,
        edgecolor="k",
        linewidth=0.45
    )
    sc = ax1.scatter(
        s.longitud, s.latitud,
        c=s.profundidad_max_m,
        cmap="viridis",
        s=100,
        edgecolor="k",
        linewidth=0.50
    )

    for ax, extent, title in [
        (ax0, regional_extent, "Contexto regional (sin cartografía)"),
        (ax1, detail_extent, "Detalle de estaciones (sin cartografía)")
    ]:
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_xlabel("Longitud [°]")
        ax.set_ylabel("Latitud [°]")
        ax.set_title(title)
        ax.grid(alpha=0.18)

    for _, r in s.iterrows():
        ax0.annotate(
            f"E{int(r.numero_estacion):02d}",
            (r.longitud, r.latitud),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=10
        )
        ax1.annotate(
            f"E{int(r.numero_estacion):02d}",
            (r.longitud, r.latitud),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=11
        )

    cb = fig.colorbar(sc, ax=axs, pad=0.025, shrink=0.88)
    cb.set_label("Profundidad máxima muestreada [m]")

    fig.suptitle(
        "Estaciones CTD — instalar Cartopy o Basemap para mostrar costa real",
        fontsize=CFG["font_suptitle"],
        weight="bold"
    )

    savefig(
        fig,
        out / "Fig01_V3_4_mapa_estaciones_sin_cartografia.png"
    )

def fig_ts_vf(base,nucleos,o2col,out):
    p=filtro_hidrografico_pass(base); s=filtro_hidrografico_suspect(base)
    if p.empty: p=filtro_hidrografico(base)
    xlim=_robust_limits(p["salinidad_absoluta_g_kg"]); ylim=_robust_limits(p["temperatura_conservativa_C"])
    if xlim is None or ylim is None:return
    # Asegura que SWT queden visibles, sin dejar que outliers determinen la escala.
    xlim=(min(xlim[0],float(SWT.loc["SA_g_kg"].min())-.15),max(xlim[1],float(SWT.loc["SA_g_kg"].max())+.15))
    ylim=(min(ylim[0],float(SWT.loc["CT_C"].min())-.4),max(ylim[1],float(SWT.loc["CT_C"].max())+.4))
    fig,ax=plt.subplots(figsize=(8.8,7.0),layout="constrained")
    z=-np.abs(p["profundidad_m"].to_numpy(float))
    sc=ax.scatter(p["salinidad_absoluta_g_kg"],p["temperatura_conservativa_C"],c=z,s=8,alpha=.82,rasterized=True)
    if not s.empty:
        ax.scatter(s["salinidad_absoluta_g_kg"],s["temperatura_conservativa_C"],s=17,marker="x",alpha=.45,label="QARTOD 3")
    add_sigma0_fixed(ax,xlim,ylim)
    cb=fig.colorbar(sc,ax=ax,pad=.02);cb.set_label("Coordenada vertical, z [m]")
    markers={"mSUW":"o","mAAIW":"s","mNADW":"^"}
    for name in SWT.columns:
        ax.scatter(SWT.loc["SA_g_kg",name],SWT.loc["CT_C",name],marker=markers[name],s=105,
                   edgecolor="k",linewidth=.5,label=f"SWT {name}")
    # Núcleos observados: media ± 1 SD.
    for label,prefix,mask,marker in [
        ("SUW observado","SUW_core",np.ones(len(nucleos),bool),"D"),
        ("AAIW observado","AAIW_core",~nucleos["AAIW_core_truncado_fondo"].fillna(True).to_numpy(bool),"P")]:
        g=nucleos.loc[mask]
        sa=pd.to_numeric(g[f"{prefix}_SA_g_kg"],errors="coerce");ct=pd.to_numeric(g[f"{prefix}_CT_C"],errors="coerce")
        ax.errorbar(sa.mean(),ct.mean(),xerr=sa.std(ddof=1),yerr=ct.std(ddof=1),fmt=marker,ms=8,
                    capsize=3,elinewidth=1.0,label=label)
    ax.set(xlim=xlim,ylim=ylim,xlabel=r"Salinidad absoluta, $S_A$ [g kg$^{-1}$]",
           ylabel=r"Temperatura conservativa, $\Theta$ [°C]",
           title="Diagrama Θ–SA de campaña: QC y núcleos observados")
    ax.grid(alpha=.18);ax.legend(fontsize=8,ncol=2,loc="best")
    savefig(fig,out/"Fig02_VF_TS_profundidad_QC.png")

    if o2col is not None:
        op=filtro_o2_pass(base,o2col); os=filtro_o2_suspect(base,o2col)
        if not op.empty:
            fig,ax=plt.subplots(figsize=(8.8,7.0),layout="constrained")
            clim=_robust_limits(op[o2col],.01,.99,.02)
            kw={"vmin":clim[0],"vmax":clim[1]} if clim else {}
            sc=ax.scatter(op["salinidad_absoluta_g_kg"],op["temperatura_conservativa_C"],c=op[o2col],
                          s=8,alpha=.82,rasterized=True,**kw)
            if not os.empty:
                ax.scatter(os["salinidad_absoluta_g_kg"],os["temperatura_conservativa_C"],s=17,marker="x",alpha=.45,label="O₂ QARTOD 3")
            add_sigma0_fixed(ax,xlim,ylim)
            cb=fig.colorbar(sc,ax=ax,pad=.02);cb.set_label(r"O$_2$ [$\mu$mol kg$^{-1}$]")
            for name in SWT.columns:
                ax.scatter(SWT.loc["SA_g_kg",name],SWT.loc["CT_C",name],marker=markers[name],s=105,
                           edgecolor="k",linewidth=.5,label=f"SWT {name}")
            ax.set(xlim=xlim,ylim=ylim,xlabel=r"Salinidad absoluta, $S_A$ [g kg$^{-1}$]",
                   ylabel=r"Temperatura conservativa, $\Theta$ [°C]",
                   title="Diagrama Θ–SA: oxígeno como trazador auxiliar")
            ax.grid(alpha=.18);ax.legend(fontsize=8,ncol=2,loc="best")
            savefig(fig,out/"Fig03_VF_TS_oxigeno_QC.png")


def fig_capas_superficiales_vf(nucleos,out):
    if nucleos.empty:return
    fig,axs=plt.subplots(1,2,figsize=(12,5.3),layout="constrained")
    for periodo,g in nucleos.groupby("periodo"):
        g=g.sort_values("numero_estacion")
        axs[0].plot(g.numero_estacion,-g.MLD_m,marker="o",label=f"MLD {periodo}")
        axs[0].plot(g.numero_estacion,-g.ILD_m,marker="s",label=f"ILD {periodo}")
        axs[1].plot(g.numero_estacion,g.BLT_m,marker="o",label=periodo)
    axs[0].set(xlabel="Número de estación",ylabel="Coordenada vertical, z [m]",title="Profundidad de MLD e ILD")
    axs[1].set(xlabel="Número de estación",ylabel="BLT = ILD − MLD [m]",title="Espesor de capa barrera")
    for ax in axs: ax.grid(alpha=.2);ax.legend(fontsize=10)
    fig.suptitle("Estructura de la capa superficial",fontsize=14)
    savefig(fig,out/"Fig04_VF_MLD_ILD_BLT.png")


def fig_nucleos_vf(nucleos,out):
    if nucleos.empty:return
    fig,axs=plt.subplots(1,2,figsize=(12,5.5),layout="constrained",sharex=True)
    for periodo,g in nucleos.groupby("periodo"):
        g=g.sort_values("numero_estacion")
        axs[0].plot(g.numero_estacion,-g.SUW_core_depth_m,marker="o",label=periodo)
        good=g[~g.AAIW_core_truncado_fondo.fillna(True)]
        bad=g[g.AAIW_core_truncado_fondo.fillna(True)]
        axs[1].plot(good.numero_estacion,-good.AAIW_core_depth_m,marker="s",label=periodo)
        if len(bad): axs[1].scatter(bad.numero_estacion,-bad.AAIW_core_depth_m,marker="x",s=55,label=f"{periodo} truncado")
    axs[0].set(title="Núcleo SUW (máximo de $S_A$)",ylabel="Coordenada vertical, z [m]",xlabel="Número de estación")
    axs[1].set(title="Núcleo AAIW (mínimo de $S_A$)",xlabel="Número de estación")
    for ax in axs:ax.grid(alpha=.2);ax.legend(fontsize=10)
    fig.suptitle("Profundidad de núcleos hidrográficos",fontsize=14)
    savefig(fig,out/"Fig05_VF_nucleos_SUW_AAIW.png")


def _section_matrix_stations(df,value_col,periodo,stations,dz=None):
    dz=CFG["seccion_dz_m"] if dz is None else dz
    x=df[df.periodo.eq(periodo) & df.numero_estacion.isin(stations)].copy()
    if x.empty or value_col not in x:return None
    pos=(x.groupby("numero_estacion",as_index=False).agg(latitud=("latitud","median"),longitud=("longitud","median")))
    order={st:i for i,st in enumerate(stations)};pos["_o"]=pos.numero_estacion.map(order);pos=pos.sort_values("_o")
    dist=[0.0]
    for i in range(1,len(pos)):
        dist.append(dist[-1]+haversine_km(pos.iloc[i-1].latitud,pos.iloc[i-1].longitud,pos.iloc[i].latitud,pos.iloc[i].longitud))
    pos["dist_km"]=dist
    zmax=float(np.nanmax(x.profundidad_m));zgrid=np.arange(0,zmax+dz,dz)
    M=np.full((len(zgrid),len(pos)),np.nan)
    labels=[]
    for j,r in pos.reset_index(drop=True).iterrows():
        g=x[x.numero_estacion.eq(r.numero_estacion)][["profundidad_m",value_col]].dropna().groupby("profundidad_m",as_index=False)[value_col].mean().sort_values("profundidad_m")
        if len(g)>=2:
            z=g.profundidad_m.to_numpy(float);v=g[value_col].to_numpy(float);inside=(zgrid>=z.min())&(zgrid<=z.max());M[inside,j]=np.interp(zgrid[inside],z,v)
        labels.append(f"E{int(r.numero_estacion):02d}")
    return pos.dist_km.to_numpy(float),zgrid,M,labels,pos


def fig_secciones_hidrograficas_vf(base,o2col,out):
    hydro=filtro_hidrografico_pass(base)
    oxy=filtro_o2_pass(base,o2col) if o2col else hydro.iloc[0:0]
    variables=[("salinidad_absoluta_g_kg",r"$S_A$ [g kg$^{-1}$]",hydro),
               ("temperatura_conservativa_C",r"$\Theta$ [°C]",hydro),
               ("sigma0_kg_m3",r"$\sigma_0$ [kg m$^{-3}$]",hydro)]
    if o2col: variables.append((o2col,r"O$_2$ [$\mu$mol kg$^{-1}$]",oxy))
    periods=[p for p in ("Diurno","Nocturno") if (base.periodo==p).any()]
    for tname,stations in CFG["transectos"].items():
        # omite agrupaciones con menos de dos estaciones presentes
        present=set(pd.to_numeric(base.numero_estacion,errors="coerce").dropna().astype(int))
        use=[s for s in stations if s in present]
        if len(use)<2:continue
        fig,axs=plt.subplots(len(periods),len(variables),figsize=(4.4*len(variables),4.4*len(periods)),
                             layout="constrained",sharey="row",squeeze=False)
        for i,periodo in enumerate(periods):
            for j,(col,label,source) in enumerate(variables):
                ax=axs[i,j];sec=_section_matrix_stations(source,col,periodo,use)
                if sec is None:ax.set_visible(False);continue
                xd,z,M,labs,_=sec
                pc=ax.pcolormesh(xd,-z,M,shading="nearest")
                cb=fig.colorbar(pc,ax=ax,pad=.02);cb.set_label(label)
                ax.set_xticks(xd);ax.set_xticklabels(labs);ax.grid(alpha=.12)
                if i==len(periods)-1:ax.set_xlabel("Distancia dentro del grupo [km]")
                if j==0:ax.set_ylabel(f"{periodo}\nz [m]")
                if i==0:ax.set_title(label)
        fig.suptitle(f"Sección hidrográfica · Grupo {tname} · estaciones {', '.join(f'E{s:02d}' for s in use)}\nAgrupación configurable: confirmar con la bitácora",fontsize=13)
        savefig(fig,out/f"Fig06_VF_seccion_{slug(tname)}.png")


def fig_omp_sensibilidad_vf(omp_full,omp_ts,out):
    if omp_full.empty or omp_ts.empty:return
    # Promedios por bins de 25 m.
    dz=25.0
    a=omp_full[omp_full.omp_ok].copy();b=omp_ts[omp_ts.omp_ok].copy()
    a["zbin"]=np.round(a.profundidad_m/dz)*dz;b["zbin"]=np.round(b.profundidad_m/dz)*dz
    ma=a.groupby("zbin")[["frac_mSUW","frac_mAAIW","frac_mNADW"]].mean()
    mb=b.groupby("zbin")[["frac_mSUW","frac_mAAIW","frac_mNADW"]].mean()
    fig,axs=plt.subplots(1,3,figsize=(13,6.2),layout="constrained",sharey=True)
    for ax,mass in zip(axs,["mSUW","mAAIW","mNADW"]):
        c=f"frac_{mass}"
        if c in ma.columns: ax.plot(100*ma[c],-ma.index,label="Θ+SA+O₂")
        if c in mb.columns: ax.plot(100*mb[c],-mb.index,"--",label="Θ+SA (sensibilidad)")
        ax.axvline(30,ls=":",lw=1);ax.set(xlabel="Fracción media [%]",title=mass);ax.grid(alpha=.2);ax.legend(fontsize=8)
    axs[0].set_ylabel("Coordenada vertical, z [m]")
    fig.suptitle("Sensibilidad de las fracciones OMP al uso de O₂\nPromedio muestral; no equivale a fracción volumétrica del Caribe",fontsize=13)
    savefig(fig,out/"Fig07_VF_OMP_sensibilidad_vertical.png")


def fig_omp_residuales_vf(omp_diag,out):
    if omp_diag.empty:return
    cols=[("resid_CT_norm_C",r"Residual $\Theta$ [°C]"),("resid_SA_norm_g_kg",r"Residual $S_A$ [g kg$^{-1}$]"),("resid_O2_norm_umol_kg",r"Residual O$_2$ [$\mu$mol kg$^{-1}$]")]
    fig,axs=plt.subplots(1,3,figsize=(13,6),layout="constrained",sharey=True)
    dz=25.0
    q=omp_diag.copy();q["zbin"]=np.round(q.profundidad_m/dz)*dz
    for ax,(c,label) in zip(axs,cols):
        g=q.groupby("zbin")[c]
        med=g.median();q25=g.quantile(.25);q75=g.quantile(.75)
        ax.fill_betweenx(-med.index,q25,q75,alpha=.2,label="IQR")
        ax.plot(med,-med.index,lw=1.5,label="Mediana")
        ax.axvline(0,ls="--",lw=1);ax.set(xlabel=label,title=label);ax.grid(alpha=.2);ax.legend(fontsize=8)
    axs[0].set_ylabel("Coordenada vertical, z [m]")
    fig.suptitle("Diagnóstico de residuales del OMP completo (fracciones normalizadas)",fontsize=13)
    savefig(fig,out/"Fig08_VF_OMP_residuales.png")


def fig_omp_estandar_vf(omp_diag,omp_ts,out):
    if omp_diag.empty or omp_ts.empty:return
    full=_interp_omp_standard(omp_diag,["frac_mSUW","frac_mAAIW","frac_mNADW"])
    ts=_interp_omp_standard(omp_ts,["frac_mSUW","frac_mAAIW","frac_mNADW"])
    if full.empty or ts.empty:return
    rows=[]
    for z in CFG["profundidades_estandar_m"]:
        rf=full[full.prof_objetivo_m.eq(z)];rt=ts[ts.prof_objetivo_m.eq(z)]
        for mass in ["mSUW","mAAIW","mNADW"]:
            rows.append({"profundidad_m":z,"masa":mass,
                         "full_media":100*rf[f"frac_{mass}"].mean() if len(rf) else np.nan,
                         "ts_media":100*rt[f"frac_{mass}"].mean() if len(rt) else np.nan})
    d=pd.DataFrame(rows)
    fig,axs=plt.subplots(1,3,figsize=(13,5.2),layout="constrained",sharey=True)
    for ax,mass in zip(axs,["mSUW","mAAIW","mNADW"]):
        g=d[d.masa.eq(mass)].sort_values("profundidad_m")
        ax.plot(g.full_media,-g.profundidad_m,marker="o",label="Θ+SA+O₂")
        ax.plot(g.ts_media,-g.profundidad_m,marker="s",ls="--",label="Θ+SA")
        ax.set(xlabel="Fracción media [%]",title=mass);ax.grid(alpha=.2);ax.legend(fontsize=8)
    axs[0].set_ylabel("Profundidad estándar, z [m]")
    fig.suptitle("Comparación OMP en profundidades estándar",fontsize=13)
    savefig(fig,out/"Fig09_VF_OMP_profundidades_estandar.png")


def write_captions_vf(out,nucleos,omp_stats,omp_status):
    summary=resumen_hidrografico_masas(nucleos)
    vals={r.masa:r for _,r in summary.iterrows()}
    suw=vals.get("SUW");aa=vals.get("AAIW")
    txt="# Leyendas y notas para figuras VF\n\n"
    txt += "**Fig. 01.** Localización de las estaciones CTD y profundidad máxima muestreada. Las líneas discontinuas representan agrupaciones espaciales configurables y deben verificarse con la bitácora del crucero antes de denominarlas transectos.\n\n"
    txt += "**Fig. 02.** Diagrama Θ–SA de campaña. QARTOD PASS controla la nube principal y los límites de ejes; QARTOD 3 se conserva como diagnóstico. Se incluyen los SWT de Torres et al. (2023) y las medias ±1 SD de los núcleos hidrográficos observados.\n\n"
    txt += "**Fig. 03.** Diagrama Θ–SA coloreado por O₂. El oxígeno se interpreta como trazador auxiliar mientras no exista validación absoluta independiente con Winkler.\n\n"
    txt += "**Fig. 04.** MLD, ILD y espesor de capa barrera por estación para muestreos diurnos y nocturnos.\n\n"
    txt += "**Fig. 05.** Profundidad de los núcleos SUW y AAIW. Los mínimos de AAIW truncados por el fondo se señalan y no se usan en la estadística del núcleo.\n\n"
    txt += "**Fig. 06.** Secciones hidrográficas por agrupaciones espaciales; evitar conectar estaciones de grupos distintos.\n\n"
    txt += "**Fig. 07.** Sensibilidad de las fracciones OMP a la inclusión de O₂. La solución Θ+SA se presenta como diagnóstico de sensibilidad, no como sustituto automático del OMP completo.\n\n"
    txt += "**Fig. 08.** Residuales físicos del OMP completo reconstruidos con fracciones normalizadas. Un sesgo de O₂ persistente indica que el trazador no es reproducido por los SWT/mediciones bajo los supuestos del modelo.\n\n"
    txt += "**Fig. 09.** Comparación de fracciones OMP completas y Θ+SA en profundidades estándar.\n\n"
    if suw is not None:
        txt += f"SUW observado: z = {suw['profundidad_m_media']:.1f} ± {suw['profundidad_m_sd']:.1f} m; SA = {suw['SA_g_kg_media']:.3f} ± {suw['SA_g_kg_sd']:.3f} g kg⁻¹; Θ = {suw['CT_C_media']:.2f} ± {suw['CT_C_sd']:.2f} °C.\n\n"
    if aa is not None:
        txt += f"AAIW observado (núcleos no truncados): z = {aa['profundidad_m_media']:.1f} ± {aa['profundidad_m_sd']:.1f} m; SA = {aa['SA_g_kg_media']:.3f} ± {aa['SA_g_kg_sd']:.3f} g kg⁻¹; Θ = {aa['CT_C_media']:.2f} ± {aa['CT_C_sd']:.2f} °C.\n\n"
    txt += f"Estado OMP: **{omp_status}**.\n"
    if not omp_stats.empty:
        txt += "\nResumen de residuales (texto):\n\n```text\n" + omp_stats.to_string(index=False) + "\n```\n"
    (out/"FIGURAS_VF_LEYENDAS.md").write_text(txt,encoding="utf-8")



# =============================================================================
# V2.0 — SÍNTESIS HIDROGRÁFICA DE CINCO NIVELES + FIGURAS DE PUBLICACIÓN
# =============================================================================

# La V2.0 conserva la arquitectura numérica de V1.1 y añade una separación
# explícita entre:
#   (a) cinco niveles hidrográficos OBSERVABLES: CSW, SUW, aguas centrales,
#       AAIW y NADW; y
#   (b) tres source water types INDEPENDIENTES para OMP: mSUW, mAAIW, mNADW.
# Las aguas centrales no se fuerzan como end-member independiente.

VERSION = "2.0"
CFG.update({
    "carpeta_salida": DEFAULT_SESSION / "ANALISIS_MASAS_AGUA_V2_0",
    "font_base": 13,
    "font_axes": 14,
    "font_title": 16,
    "font_suptitle": 18,
    "font_legend": 11.5,
    "ts_marker_size": 30,
    "ts_marker_alpha": 0.82,
    "ts_marker_edge_lw": 0.45,
    "csw_max_depth_operational_m": 75.0,
    "suw_core_halfwindow_m": 45.0,
    "aaiw_core_halfwindow_m": 110.0,
    "central_suw_offset_m": 45.0,
    "central_aaiw_offset_m": 110.0,
    # Los siguientes umbrales son controles OPERACIONALES de confianza,
    # no definiciones universales de masas de agua.
    "suw_prominence_high_gkg": 0.05,
    "suw_prominence_medium_gkg": 0.02,
    "aaiw_prominence_high_gkg": 0.03,
    "aaiw_prominence_medium_gkg": 0.015,
    "nadw_min_depth_support_m": 1000.0,
    "nadw_resolved_depth_m": 1200.0,
    "nadw_sa_recovery_gkg": 0.05,
    "nadw_o2_recovery_umolkg": 10.0,
    "bootstrap_n": 2000,
    "bootstrap_seed": 2026,
})


def apply_publication_style():
    """Tipografía legible para figuras insertadas en Word y manuscrito."""
    plt.rcParams.update({
        "font.size": CFG["font_base"],
        "axes.labelsize": CFG["font_axes"],
        "axes.titlesize": CFG["font_title"],
        "figure.titlesize": CFG["font_suptitle"],
        "legend.fontsize": CFG["font_legend"],
        "xtick.labelsize": CFG["font_base"],
        "ytick.labelsize": CFG["font_base"],
        "axes.linewidth": 1.05,
        "lines.linewidth": 1.6,
        "savefig.dpi": CFG["dpi"],
    })


def analizar_perfil_hidrografico_v2(g: pd.DataFrame, pcol: str, o2col: str | None,
                                     profile_id=None) -> dict:
    """Núcleos con PASS como fuente primaria y fallback al filtro hidrográfico aceptado."""
    gp = filtro_hidrografico_pass(g)
    fuente = "QARTOD_PASS"
    if len(gp) < 5:
        gp = filtro_hidrografico(g)
        fuente = "PASS_SUSPECT_FALLBACK"
    # conservar metadatos que pueden haberse perdido por filtrado
    if gp.empty:
        gp = g.copy()
    res = analizar_perfil_hidrografico(gp, pcol, o2col)
    res["id_perfil"] = profile_id if profile_id is not None else g.get("id_perfil", pd.Series([res.get("estacion", "")])).iloc[0]
    res["fuente_QC_nucleos"] = fuente
    return res


def _local_prominence(profile: pd.DataFrame, z0: float, value_col: str, kind="max",
                      inner=25.0, outer=100.0):
    """Prominencia robusta respecto a hombros por encima y debajo del núcleo."""
    if not np.isfinite(z0) or value_col not in profile:
        return np.nan, np.nan, np.nan
    q = profile[["profundidad_m", value_col]].dropna().sort_values("profundidad_m")
    if q.empty:
        return np.nan, np.nan, np.nan
    z = q["profundidad_m"]
    above = q[(z >= z0-outer) & (z <= z0-inner)][value_col]
    below = q[(z >= z0+inner) & (z <= z0+outer)][value_col]
    core = safe_interp(q["profundidad_m"], q[value_col], z0)
    va = float(above.median()) if len(above) else np.nan
    vb = float(below.median()) if len(below) else np.nan
    diffs=[]
    if np.isfinite(va): diffs.append(core-va if kind=="max" else va-core)
    if np.isfinite(vb): diffs.append(core-vb if kind=="max" else vb-core)
    prom=float(min(diffs)) if diffs else np.nan
    return prom, va, vb


def _confidence_from_prominence(prom, high, medium, truncated=False, fallback=False):
    if truncated:
        return "BAJA_TRUNCADO"
    if not np.isfinite(prom):
        return "BAJA_SIN_HOMBROS" if not fallback else "BAJA_QC_FALLBACK"
    if prom >= high and not fallback:
        return "ALTA"
    if prom >= medium:
        return "MEDIA" if not fallback else "MEDIA_QC_FALLBACK"
    return "BAJA"


def enriquecer_nucleos_v2(base: pd.DataFrame, nucleos: pd.DataFrame, omp_diag: pd.DataFrame,
                           o2col: str | None, profile_col: str) -> pd.DataFrame:
    """Añade confianza de núcleos, aguas centrales y evidencia de NADW por perfil."""
    if nucleos.empty:
        return nucleos
    out=nucleos.copy()
    rows=[]
    omp_groups={str(k):v.copy() for k,v in omp_diag.groupby("id_perfil",sort=False)} if not omp_diag.empty else {}
    for _,r in out.iterrows():
        pid=str(r.get("id_perfil",r.get("estacion","")))
        g=base[base[profile_col].astype(str).eq(pid)].copy()
        gp=filtro_hidrografico_pass(g)
        fallback=False
        if len(gp)<5:
            gp=filtro_hidrografico(g); fallback=True
        suw_z=float(pd.to_numeric(pd.Series([r.get("SUW_core_depth_m",np.nan)]),errors="coerce").iloc[0])
        aa_z=float(pd.to_numeric(pd.Series([r.get("AAIW_core_depth_m",np.nan)]),errors="coerce").iloc[0])
        trunc=bool(r.get("AAIW_core_truncado_fondo",True))
        psuw,_,_=_local_prominence(gp,suw_z,"salinidad_absoluta_g_kg","max")
        paa,_,_=_local_prominence(gp,aa_z,"salinidad_absoluta_g_kg","min",inner=35,outer=150)
        suw_conf=_confidence_from_prominence(psuw,CFG["suw_prominence_high_gkg"],CFG["suw_prominence_medium_gkg"],False,fallback)
        aa_conf=_confidence_from_prominence(paa,CFG["aaiw_prominence_high_gkg"],CFG["aaiw_prominence_medium_gkg"],trunc,fallback)

        central_top=np.nan; central_bottom=np.nan
        central_sa_min=np.nan; central_sa_max=np.nan; central_ct_min=np.nan; central_ct_max=np.nan; central_o2_min=np.nan
        if np.isfinite(suw_z) and np.isfinite(aa_z) and aa_z>suw_z:
            central_top=suw_z+CFG["central_suw_offset_m"]
            central_bottom=aa_z-CFG["central_aaiw_offset_m"]
            if central_bottom<=central_top:
                central_top=suw_z
                central_bottom=aa_z
            cg=gp[(gp.profundidad_m>=central_top)&(gp.profundidad_m<=central_bottom)]
            if len(cg):
                central_sa_min=float(cg.salinidad_absoluta_g_kg.min()); central_sa_max=float(cg.salinidad_absoluta_g_kg.max())
                central_ct_min=float(cg.temperatura_conservativa_C.min()); central_ct_max=float(cg.temperatura_conservativa_C.max())
                if o2col and o2col in cg and cg[o2col].notna().any(): central_o2_min=float(cg[o2col].min())

        zmax=float(pd.to_numeric(gp.profundidad_m,errors="coerce").max()) if len(gp) else np.nan
        nadw_state="NADW_NO_EVALUADA"
        nadw_maxfrac=np.nan; nadw_top=np.nan; sa_rec=np.nan; o2_rec=np.nan
        if not np.isfinite(aa_z):
            nadw_state="NADW_NO_EVALUADA_SIN_AAIW"
        elif not np.isfinite(zmax) or zmax < max(CFG["nadw_min_depth_support_m"],aa_z+150):
            nadw_state="NADW_NO_ALCANZADA"
        else:
            deep=gp[gp.profundidad_m>=aa_z+100].copy()
            aa_sa=float(r.get("AAIW_core_SA_g_kg",np.nan)); aa_o2=float(r.get("AAIW_core_O2_umol_kg",np.nan))
            if len(deep) and np.isfinite(aa_sa): sa_rec=float(deep.salinidad_absoluta_g_kg.tail(max(1,min(20,len(deep)))).median()-aa_sa)
            if o2col and o2col in deep and deep[o2col].notna().any() and np.isfinite(aa_o2):
                o2_rec=float(deep[o2col].dropna().tail(max(1,min(20,deep[o2col].notna().sum()))).median()-aa_o2)
            og=omp_groups.get(pid,pd.DataFrame())
            if not og.empty:
                oq=og[og.get("omp_ok",False).astype(bool) & (og.profundidad_m>=aa_z+100)]
                if len(oq):
                    nadw_maxfrac=float(oq.frac_mNADW.max())
                    sig=oq[oq.frac_mNADW>=CFG["fraccion_significativa"]]
                    if len(sig): nadw_top=float(sig.profundidad_m.min())
            sa_support=np.isfinite(sa_rec) and sa_rec>=CFG["nadw_sa_recovery_gkg"]
            o2_support=(not np.isfinite(o2_rec)) or (o2_rec>=CFG["nadw_o2_recovery_umolkg"])
            omp_support=np.isfinite(nadw_maxfrac) and nadw_maxfrac>=CFG["fraccion_significativa"]
            if zmax>=CFG["nadw_resolved_depth_m"] and sa_support and o2_support and omp_support:
                nadw_state="NADW_RESUELTA"
            elif sa_support or omp_support:
                nadw_state="NADW_INFLUENCIA"
            else:
                nadw_state="NADW_NO_RESUELTA"

        rows.append({
            "id_perfil":r.get("id_perfil",pid),
            "SUW_prominencia_SA_g_kg":psuw,"SUW_confianza":suw_conf,
            "AAIW_prominencia_SA_g_kg":paa,"AAIW_confianza":aa_conf,
            "Central_top_m":central_top,"Central_bottom_m":central_bottom,
            "Central_SA_min_g_kg":central_sa_min,"Central_SA_max_g_kg":central_sa_max,
            "Central_CT_min_C":central_ct_min,"Central_CT_max_C":central_ct_max,
            "Central_O2_min_umol_kg":central_o2_min,
            "NADW_estado":nadw_state,"NADW_top_aprox_m":nadw_top,"NADW_max_frac_OMP":nadw_maxfrac,
            "NADW_recuperacion_SA_g_kg":sa_rec,"NADW_recuperacion_O2_umol_kg":o2_rec,
        })
    add=pd.DataFrame(rows)
    return out.merge(add,on="id_perfil",how="left")


def _bootstrap_ci(values, nboot=None, seed=None):
    v=pd.to_numeric(pd.Series(values),errors="coerce").dropna().to_numpy(float)
    if len(v)<2: return (np.nan,np.nan)
    nboot=CFG["bootstrap_n"] if nboot is None else nboot
    rng=np.random.default_rng(CFG["bootstrap_seed"] if seed is None else seed)
    idx=rng.integers(0,len(v),size=(nboot,len(v)))
    med=np.median(v[idx],axis=1)
    return float(np.quantile(med,.025)),float(np.quantile(med,.975))


def resumen_hidrografico_robusto_v2(nucleos: pd.DataFrame) -> pd.DataFrame:
    """Síntesis para manuscrito: mediana, IQR e IC95 bootstrap de los núcleos."""
    rows=[]
    specs=[("SUW","SUW_core",np.ones(len(nucleos),bool)),
           ("AAIW","AAIW_core",~nucleos["AAIW_core_truncado_fondo"].fillna(True).to_numpy(bool))]
    for mass,prefix,mask in specs:
        g=nucleos.loc[mask].copy()
        row={"masa":mass,"n_perfiles":len(g)}
        for outc,inc in {"profundidad_m":f"{prefix}_depth_m","SA_g_kg":f"{prefix}_SA_g_kg","CT_C":f"{prefix}_CT_C","sigma0":f"{prefix}_sigma0","O2_umol_kg":f"{prefix}_O2_umol_kg"}.items():
            v=pd.to_numeric(g.get(inc),errors="coerce").dropna()
            row[f"{outc}_media"]=float(v.mean()) if len(v) else np.nan
            row[f"{outc}_sd"]=float(v.std(ddof=1)) if len(v)>1 else np.nan
            row[f"{outc}_mediana"]=float(v.median()) if len(v) else np.nan
            row[f"{outc}_p25"]=float(v.quantile(.25)) if len(v) else np.nan
            row[f"{outc}_p75"]=float(v.quantile(.75)) if len(v) else np.nan
            lo,hi=_bootstrap_ci(v)
            row[f"{outc}_IC95_mediana_lo"]=lo; row[f"{outc}_IC95_mediana_hi"]=hi
        if f"{mass}_confianza" in g:
            vc=g[f"{mass}_confianza"].value_counts(dropna=False)
            row["confianza_ALTA_n"]=int(vc.get("ALTA",0));row["confianza_MEDIA_n"]=int(sum(v for k,v in vc.items() if str(k).startswith("MEDIA")))
            row["confianza_BAJA_n"]=int(sum(v for k,v in vc.items() if str(k).startswith("BAJA")))
        rows.append(row)
    return pd.DataFrame(rows)


def preparar_dataset_masas_v2(base: pd.DataFrame, nucleos: pd.DataFrame,
                               omp_diag: pd.DataFrame, profile_col: str, o2col: str | None):
    """Asigna a cada bin un nivel hidrográfico para VISUALIZACIÓN, no un nuevo OMP."""
    x=filtro_hidrografico_pass(base)
    if x.empty: x=filtro_hidrografico(base)
    x=x.copy();x["_pid_v2"]=x[profile_col].astype(str)
    nuc=nucleos.copy();nuc["_pid_v2"]=nuc["id_perfil"].astype(str)
    keep=["_pid_v2","MLD_m","ILD_m","SUW_core_depth_m","AAIW_core_depth_m","AAIW_core_truncado_fondo","Central_top_m","Central_bottom_m","NADW_top_aprox_m","NADW_estado"]
    keep=[c for c in keep if c in nuc]
    x=x.merge(nuc[keep],on="_pid_v2",how="left")

    # Incorporar fracción mNADW por proximidad vertical sólo como apoyo diagnóstico.
    x["frac_mNADW_v2"]=np.nan;x["omp_ok_v2"]=False
    if not omp_diag.empty:
        pieces=[]
        for pid,g in x.groupby("_pid_v2",sort=False):
            h=omp_diag[omp_diag.id_perfil.astype(str).eq(pid)].sort_values("profundidad_m")
            gg=g.sort_values("profundidad_m")
            if h.empty:
                pieces.append(gg);continue
            hh=h[["profundidad_m","frac_mNADW","omp_ok"]].rename(columns={"frac_mNADW":"frac_mNADW_v2_r","omp_ok":"omp_ok_v2_r"})
            mm=pd.merge_asof(gg,hh,on="profundidad_m",direction="nearest",tolerance=1.5)
            mm["frac_mNADW_v2"]=mm["frac_mNADW_v2_r"].combine_first(mm["frac_mNADW_v2"])
            mm["omp_ok_v2"]=mm["omp_ok_v2_r"].fillna(False).astype(bool)
            mm=mm.drop(columns=["frac_mNADW_v2_r","omp_ok_v2_r"])
            pieces.append(mm)
        x=pd.concat(pieces,ignore_index=True,sort=False)

    def classify(r):
        z=float(r.get("profundidad_m",np.nan)); suw=float(r.get("SUW_core_depth_m",np.nan)); aa=float(r.get("AAIW_core_depth_m",np.nan))
        mld=float(r.get("MLD_m",np.nan));ild=float(r.get("ILD_m",np.nan))
        csw=max(CFG["csw_max_depth_operational_m"], *(v for v in [mld,ild] if np.isfinite(v))) if any(np.isfinite(v) for v in [mld,ild]) else CFG["csw_max_depth_operational_m"]
        csw=min(csw,120.0)
        if z<=csw: return "CSW"
        ctop=float(r.get("Central_top_m",np.nan));cbot=float(r.get("Central_bottom_m",np.nan))
        if np.isfinite(ctop) and z<ctop: return "SUW"
        if np.isfinite(ctop) and np.isfinite(cbot) and ctop<=z<=cbot: return "Aguas centrales"
        nadw_top=float(r.get("NADW_top_aprox_m",np.nan));nadw_state=str(r.get("NADW_estado",""))
        if (np.isfinite(nadw_top) and z>=nadw_top) or ("NADW_" in nadw_state and nadw_state in {"NADW_RESUELTA","NADW_INFLUENCIA"} and bool(r.get("omp_ok_v2",False)) and float(r.get("frac_mNADW_v2",np.nan))>=CFG["fraccion_significativa"]):
            return "NADW"
        if np.isfinite(cbot) and z>cbot: return "AAIW"
        if np.isfinite(suw) and z<=suw+CFG["suw_core_halfwindow_m"]: return "SUW"
        if np.isfinite(aa) and z<aa: return "Aguas centrales"
        if np.isfinite(aa): return "AAIW"
        return "Sin clasificar"
    x["masa_hidrografica_v2"]=x.apply(classify,axis=1)
    return x


def fig_ts_bipanel_v2(base,nucleos,omp_diag,o2col,profile_col,out):
    """Figura principal: marcador=masa; color=profundidad/O2; barras independientes."""
    from matplotlib.lines import Line2D
    from matplotlib.colors import Normalize
    x=preparar_dataset_masas_v2(base,nucleos,omp_diag,profile_col,o2col)
    if x.empty:return
    xp=filtro_hidrografico_pass(base)
    if xp.empty:xp=filtro_hidrografico(base)
    xlim=_robust_limits(xp["salinidad_absoluta_g_kg"]);ylim=_robust_limits(xp["temperatura_conservativa_C"])
    if xlim is None or ylim is None:return
    xlim=(min(xlim[0],float(SWT.loc["SA_g_kg"].min())-.15),max(xlim[1],float(SWT.loc["SA_g_kg"].max())+.15))
    ylim=(min(ylim[0],float(SWT.loc["CT_C"].min())-.4),max(ylim[1],float(SWT.loc["CT_C"].max())+.4))
    styles={"CSW":"o","SUW":"^","Aguas centrales":"s","AAIW":"D","NADW":"P"}
    order=list(styles)
    zvals=pd.to_numeric(x.profundidad_m,errors="coerce")
    znorm=Normalize(vmin=float(np.nanmin(zvals)),vmax=float(np.nanmax(zvals)))
    cmap_z=plt.cm.viridis
    xo=x[x[o2col].notna()].copy() if o2col and o2col in x else x.iloc[0:0].copy()
    if len(xo):
        ov=pd.to_numeric(xo[o2col],errors="coerce").dropna();lo,hi=np.quantile(ov,[.01,.99]) if len(ov)>=20 else (ov.min(),ov.max())
        onorm=Normalize(vmin=float(lo),vmax=float(hi));cmap_o=plt.cm.plasma
    else:onorm=None;cmap_o=plt.cm.plasma

    fig,axs=plt.subplots(1,2,figsize=(17.5,7.8),sharex=True,sharey=True,layout="constrained")
    for ax in axs:
        add_sigma0_fixed(ax,xlim,ylim);ax.set_xlim(xlim);ax.set_ylim(ylim);ax.grid(alpha=.16)
        ax.set_xlabel("Salinidad absoluta, SA [g kg⁻¹]")
    axs[0].set_ylabel("Temperatura conservativa, Θ [°C]")

    for mass in order:
        g=x[x.masa_hidrografica_v2.eq(mass)]
        if len(g):
            axs[0].scatter(g.salinidad_absoluta_g_kg,g.temperatura_conservativa_C,c=g.profundidad_m,cmap=cmap_z,norm=znorm,
                           marker=styles[mass],s=CFG["ts_marker_size"],alpha=CFG["ts_marker_alpha"],edgecolors="k",linewidths=CFG["ts_marker_edge_lw"],rasterized=True,zorder=3)
    unc=x[x.masa_hidrografica_v2.eq("Sin clasificar")]
    if len(unc):axs[0].scatter(unc.salinidad_absoluta_g_kg,unc.temperatura_conservativa_C,c=unc.profundidad_m,cmap=cmap_z,norm=znorm,s=11,marker=".",alpha=.25,rasterized=True)
    sm=plt.cm.ScalarMappable(norm=znorm,cmap=cmap_z);sm.set_array([]);cb=fig.colorbar(sm,ax=axs[0],pad=.018);cb.set_label("Profundidad [m]",fontsize=CFG["font_axes"]);cb.ax.tick_params(labelsize=CFG["font_base"])
    axs[0].set_title("a) Estructura Θ–SA y profundidad")

    if len(xo):
        for mass in order:
            g=xo[xo.masa_hidrografica_v2.eq(mass)]
            if len(g):
                axs[1].scatter(g.salinidad_absoluta_g_kg,g.temperatura_conservativa_C,c=g[o2col],cmap=cmap_o,norm=onorm,
                               marker=styles[mass],s=CFG["ts_marker_size"],alpha=CFG["ts_marker_alpha"],edgecolors="k",linewidths=CFG["ts_marker_edge_lw"],rasterized=True,zorder=3)
        sm2=plt.cm.ScalarMappable(norm=onorm,cmap=cmap_o);sm2.set_array([]);cb2=fig.colorbar(sm2,ax=axs[1],pad=.018);cb2.set_label("Oxígeno disuelto [µmol kg⁻¹]",fontsize=CFG["font_axes"]);cb2.ax.tick_params(labelsize=CFG["font_base"])
        axs[1].set_title("b) Estructura Θ–SA y oxígeno disuelto")
    else:
        axs[1].text(.5,.5,"Oxígeno no disponible o no válido",transform=axs[1].transAxes,ha="center",va="center",fontsize=CFG["font_title"])
        axs[1].set_title("b) Oxígeno disuelto")

    # SWT siempre con estrella blanca/borde negro: no compiten con las barras de color.
    for ax in axs:
        for name in SWT.columns:
            ax.scatter(SWT.loc["SA_g_kg",name],SWT.loc["CT_C",name],marker="*",s=210,facecolor="white",edgecolor="k",linewidth=1.1,zorder=7)
            ax.annotate(name,(SWT.loc["SA_g_kg",name],SWT.loc["CT_C",name]),xytext=(6,6),textcoords="offset points",fontsize=11.5,weight="bold")

    handles=[Line2D([0],[0],marker=styles[m],linestyle="None",color="k",markerfacecolor="white",markeredgecolor="k",markersize=9,label=m) for m in order]
    handles.append(Line2D([0],[0],marker="*",linestyle="None",color="k",markerfacecolor="white",markeredgecolor="k",markersize=12,label="SWT OMP (Torres et al., 2023)"))
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(0.5,-0.045),ncol=6,frameon=True,title="Interpretación hidrográfica (forma del marcador)",fontsize=CFG["font_legend"],title_fontsize=CFG["font_legend"]+0.5)
    fig.suptitle("Diagramas Θ–SA conjuntos: cinco niveles hidrográficos y trazadores",fontsize=CFG["font_suptitle"],weight="bold")
    savefig(fig,out/"Fig02_V2_TS_bipanel_profundidad_oxigeno_masas.png")


def fig_perfiles_medianos_v2(base,nucleos,o2col,out):
    """Perfil vertical mediano e IQR para SA, Θ, sigma0 y O2."""
    q=filtro_hidrografico_pass(base)
    if q.empty:q=filtro_hidrografico(base)
    if q.empty:return
    dz=10.0;q=q.copy();q["zbin"]=np.round(q.profundidad_m/dz)*dz
    vars_=[("salinidad_absoluta_g_kg","SA [g kg⁻¹]"),("temperatura_conservativa_C","Θ [°C]"),("sigma0_kg_m3","sigma0 [kg m⁻³]")]
    if o2col and o2col in q:vars_.append((o2col,"O₂ [µmol kg⁻¹]"))
    fig,axs=plt.subplots(1,len(vars_),figsize=(4.5*len(vars_),7.2),sharey=True,layout="constrained",squeeze=False);axs=axs[0]
    suw_med=float(pd.to_numeric(nucleos.SUW_core_depth_m,errors="coerce").median()) if "SUW_core_depth_m" in nucleos else np.nan
    aa_good=nucleos[~nucleos.AAIW_core_truncado_fondo.fillna(True)] if "AAIW_core_truncado_fondo" in nucleos else nucleos.iloc[0:0]
    aa_med=float(pd.to_numeric(aa_good.AAIW_core_depth_m,errors="coerce").median()) if len(aa_good) else np.nan
    for ax,(c,lab) in zip(axs,vars_):
        g=q.groupby("zbin")[c];med=g.median();p25=g.quantile(.25);p75=g.quantile(.75)
        ax.fill_betweenx(med.index,p25,p75,alpha=.20,label="IQR");ax.plot(med,med.index,lw=2,label="Mediana")
        if np.isfinite(suw_med):ax.axhline(suw_med,ls="--",lw=1.2,label="SUW mediana")
        if np.isfinite(aa_med):ax.axhline(aa_med,ls=":",lw=1.4,label="AAIW mediana")
        ax.set_xlabel(lab);ax.set_title(lab);ax.grid(alpha=.18);ax.invert_yaxis()
    axs[0].set_ylabel("Profundidad [m]");axs[-1].legend(loc="best",fontsize=CFG["font_legend"])
    fig.suptitle("Estructura vertical mediana de la campaña (banda = IQR)",weight="bold")
    savefig(fig,out/"Fig03_V2_perfiles_medianos_IQR.png")


def tabla_maestra_v2(nucleos,omp_res):
    x=nucleos.copy()
    if not omp_res.empty:
        cols=[c for c in omp_res.columns if c not in {"id_perfil","estacion","periodo","numero_estacion"}]
        if "id_perfil" in omp_res and "id_perfil" in x:
            x=x.merge(omp_res[["id_perfil"]+cols],on="id_perfil",how="left",suffixes=("","_OMP"))
    return x


def write_methodology_v2(out: Path,path_input: Path,o2col: str|None,omp_status: str):
    text=f"""# Metodología — masas de agua del Caribe colombiano — V2.0

Generado: {datetime.now().isoformat(timespec='seconds')}
Entrada: {path_input}

## Síntesis conceptual

La interpretación hidrográfica se organiza en cinco niveles observables: CSW, SUW, aguas centrales, AAIW y NADW. Esta clasificación se utiliza para describir la columna de agua y las figuras. No implica que los cinco niveles sean end-members independientes.

Para OMP se emplean únicamente mSUW, mAAIW y mNADW de Torres et al. (2023). Las aguas centrales se interpretan como región de transición/mezcla entre SUW y AAIW y no se fuerzan como cuarto end-member.

## Criterios principales

- MLD: delta sigma0 = 0.03 kg m-3 respecto a 10 m.
- ILD: |delta temperatura conservativa| = 0.2 °C respecto a 10 m.
- SUW: máximo de SA por debajo de {CFG['suw_min_depth_m']:.0f} m.
- AAIW: mínimo de SA profundo, con control operacional >= {CFG['aaiw_min_depth_m']:.0f} m y evaluación de truncamiento por fondo.
- Aguas centrales: intervalo hidrográfico entre los núcleos SUW y AAIW; no se asigna un end-member independiente.
- NADW: se reporta como no alcanzada, influencia, no resuelta o resuelta combinando profundidad alcanzada, recuperación de SA/O2 y soporte OMP.
- Los índices de prominencia y categorías ALTA/MEDIA/BAJA son controles operacionales de robustez del núcleo, no umbrales universales de definición de masa.
- OMP: mSUW + mAAIW + mNADW debajo del máximo de salinidad; residual de masa aceptable <= {CFG['omp_mass_residual_max_pct']:.0f} %.
- Fracción >= {100*CFG['fraccion_significativa']:.0f} % se usa sólo como umbral descriptivo de contribución importante.

## SWT Torres et al. (2023)

{SWT.to_string()}

Pesos OMP: {OMP_WEIGHTS}

Columna de O2: {o2col or 'NO DISPONIBLE'}
Estado diagnóstico OMP: {omp_status}

## QC e interpretación

QARTOD PASS controla la detección de núcleos y las figuras principales. SUSPECT se conserva sólo como apoyo/fallback cuando la cobertura PASS es insuficiente. Las fracciones OMP son locales al dominio muestreado y no representan porcentajes volumétricos de todo el mar Caribe.
"""
    (out/"METODOLOGIA_MASAS_AGUA_V2.md").write_text(text,encoding="utf-8")


def write_captions_v2(out:Path):
    txt="""# Leyendas sugeridas — figuras V2.0

**Figura 1.** Localización de las estaciones CTD y profundidad máxima muestreada. Las agrupaciones espaciales configuradas deben verificarse con la bitácora antes de denominarse transectos definitivos.

**Figura 2.** Diagramas temperatura conservativa–salinidad absoluta (Θ–SA) de la campaña. (a) Observaciones coloreadas por profundidad; (b) observaciones coloreadas por oxígeno disuelto. La forma de los marcadores representa la interpretación hidrográfica en cinco niveles: CSW, SUW, aguas centrales, AAIW y NADW. Las estrellas corresponden a los source water types mSUW, mAAIW y mNADW empleados en el OMP según Torres et al. (2023). Las formas representan niveles hidrográficos observables y no implican cinco end-members independientes.

**Figura 3.** Perfiles verticales medianos de SA, Θ, sigma0 y O2 para los datos QARTOD PASS; las bandas corresponden al rango intercuartílico. Se muestran como referencia las profundidades medianas de los núcleos SUW y AAIW cuando están resueltos.

**Figuras posteriores.** Las secciones y diagnósticos OMP conservan la interpretación local del crucero. El OMP Θ+SA se usa como sensibilidad al O2 y no sustituye automáticamente al OMP completo.
"""
    (out/"FIGURAS_V2_LEYENDAS.md").write_text(txt,encoding="utf-8")


def main():
    apply_publication_style()
    print("="*82)
    print("ANÁLISIS DE MASAS DE AGUA — V2.0 PARA PUBLICACIÓN")
    print("5 niveles hidrográficos observables + OMP mSUW/mAAIW/mNADW + sensibilidad")
    print("="*82)
    path=choose_input_file(Path(CFG["archivo_entrada"]))
    out=Path(CFG["carpeta_salida"]);out.mkdir(parents=True,exist_ok=True)
    figout=out/"figuras_publicacion";figout.mkdir(exist_ok=True)
    base,decisiones,pcol,o2col,profile_col=cargar_dataset(path)
    print(f"Perfiles: {base[profile_col].nunique()} | bins: {len(base):,} | O2: {o2col}")

    nuc=[]
    for pid,g in base.groupby(profile_col,sort=False):
        nuc.append(analizar_perfil_hidrografico_v2(g,pcol,o2col,profile_id=pid))
    nucleos=pd.DataFrame(nuc)

    # OMP principal (Θ+SA+O2) y análisis de sensibilidad (Θ+SA).
    omp,excl=ejecutar_omp(base,nucleos,o2col,profile_col)
    omp_ts=ejecutar_omp_ts(base,nucleos,profile_col)
    omp_diag,omp_stats=augmentar_diagnostico_omp(omp,omp_ts)

    omp_status="NO_EVALUADO"
    if not omp_stats.empty and "resid_O2_norm_umol_kg" in set(omp_stats.variable):
        rr=omp_stats[omp_stats.variable.eq("resid_O2_norm_umol_kg")].iloc[0]
        if abs(rr.mediana)>=CFG["omp_o2_bias_warn_umol_kg"] and rr.fraccion_mismo_signo>=CFG["omp_same_sign_warn_fraction"]:
            omp_status="EXPLORATORIO_DESACOPLE_O2_SISTEMATICO"
        else:
            omp_status="SIN_ALARMA_SISTEMATICA_O2"

    # Enriquecimiento hidrográfico posterior al OMP: confianza, aguas centrales, NADW.
    nucleos=enriquecer_nucleos_v2(base,nucleos,omp_diag,o2col,profile_col)
    omp_res=resumen_omp_por_perfil(omp_diag) if not omp_diag.empty else pd.DataFrame()
    omp_std=omp_en_profundidades_estandar(omp_diag) if not omp_diag.empty else pd.DataFrame()
    ts_std=_interp_omp_standard(omp_ts,["frac_mSUW","frac_mAAIW","frac_mNADW"])
    hydro_summary=resumen_hidrografico_robusto_v2(nucleos)
    master=tabla_maestra_v2(nucleos,omp_res)

    sens=[]
    if not omp_std.empty or not ts_std.empty:
        for z in CFG["profundidades_estandar_m"]:
            a=omp_std[omp_std.prof_objetivo_m.eq(z)] if not omp_std.empty else pd.DataFrame()
            b=ts_std[ts_std.prof_objetivo_m.eq(z)] if not ts_std.empty else pd.DataFrame()
            for mass in ["mSUW","mAAIW","mNADW"]:
                sens.append({"profundidad_m":z,"masa":mass,
                             "OMP_full_media":float(a[f"frac_{mass}"].mean()) if len(a) else np.nan,
                             "OMP_TS_media":float(b[f"frac_{mass}"].mean()) if len(b) else np.nan,
                             "diferencia_full_menos_TS":float(a[f"frac_{mass}"].mean()-b[f"frac_{mass}"].mean()) if len(a) and len(b) else np.nan})
    sens=pd.DataFrame(sens)

    # Dataset de clasificación usado en la figura Θ-SA; se exporta para trazabilidad.
    clasif=preparar_dataset_masas_v2(base,nucleos,omp_diag,profile_col,o2col)

    # CSV reproducibles.
    master.to_csv(out/"MASAS_AGUA_POR_PERFIL_V2.csv",index=False,encoding="utf-8-sig")
    hydro_summary.to_csv(out/"RESUMEN_HIDROGRAFICO_ROBUSTO_V2.csv",index=False,encoding="utf-8-sig")
    clasif.to_csv(out/"CLASIFICACION_HIDROGRAFICA_POR_BIN_V2.csv",index=False,encoding="utf-8-sig")
    if not omp_diag.empty:omp_diag.to_csv(out/"OMP_FULL_DIAGNOSTICO_V2.csv",index=False,encoding="utf-8-sig")
    if not omp_ts.empty:omp_ts.to_csv(out/"OMP_TS_SENSIBILIDAD_V2.csv",index=False,encoding="utf-8-sig")
    if not omp_stats.empty:omp_stats.to_csv(out/"OMP_RESIDUALES_RESUMEN_V2.csv",index=False,encoding="utf-8-sig")
    if not sens.empty:sens.to_csv(out/"OMP_SENSIBILIDAD_PROF_ESTANDAR_V2.csv",index=False,encoding="utf-8-sig")
    if not excl.empty:excl.to_csv(out/"OMP_PERFILES_EXCLUIDOS_V2.csv",index=False,encoding="utf-8-sig")

    # Workbook para análisis y manuscrito.
    with pd.ExcelWriter(out/"RESULTADOS_MASAS_AGUA_V2.xlsx",engine="openpyxl") as wr:
        master.to_excel(wr,sheet_name="Masas_por_perfil",index=False)
        hydro_summary.to_excel(wr,sheet_name="Resumen_robusto",index=False)
        clasif.to_excel(wr,sheet_name="Clasificacion_bins",index=False)
        if not omp_diag.empty:omp_diag.to_excel(wr,sheet_name="OMP_full",index=False)
        if not omp_ts.empty:omp_ts.to_excel(wr,sheet_name="OMP_TS_sens",index=False)
        if not omp_stats.empty:omp_stats.to_excel(wr,sheet_name="OMP_residuales",index=False)
        if not sens.empty:sens.to_excel(wr,sheet_name="OMP_prof_estandar",index=False)
        if not excl.empty:excl.to_excel(wr,sheet_name="OMP_excluidos",index=False)
        SWT.to_excel(wr,sheet_name="SWT_Torres2023")
        pd.DataFrame([CFG]).astype(str).to_excel(wr,sheet_name="Configuracion",index=False)

    # Figuras de publicación. PDF vectorial se guarda en paralelo por savefig().
    fig_mapa_estaciones_vf(base,figout)
    fig_ts_bipanel_v2(base,nucleos,omp_diag,o2col,profile_col,figout)
    fig_perfiles_medianos_v2(base,nucleos,o2col,figout)
    fig_capas_superficiales_vf(nucleos,figout)
    fig_nucleos_vf(nucleos,figout)
    fig_secciones_hidrograficas_vf(base,o2col,figout)
    fig_omp_sensibilidad_vf(omp_diag,omp_ts,figout)
    fig_omp_residuales_vf(omp_diag,figout)
    fig_omp_estandar_vf(omp_diag,omp_ts,figout)

    write_methodology_v2(out,path,o2col,omp_status)
    write_captions_v2(out)

    meta={"version":VERSION,"input":str(path),"output":str(out),"source_reference":SOURCE_REFERENCE,
          "SWT":SWT.to_dict(),"OMP_weights":OMP_WEIGHTS,"OMP_status":omp_status,
          "sintesis_hidrografica":["CSW","SUW","Aguas centrales","AAIW","NADW"],
          "endmembers_OMP":["mSUW","mAAIW","mNADW"],
          "nota":"La clasificación de cinco niveles se usa para interpretación hidrográfica/visual. Aguas centrales no se fuerza como end-member. Umbrales de prominencia/confianza son controles operacionales, no definiciones universales.",
          "transectos_configurados":CFG["transectos"]}
    (out/"metadata_masas_agua_v2.json").write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding="utf-8")

    print("\nRESUMEN HIDROGRÁFICO ROBUSTO")
    print(hydro_summary.to_string(index=False))
    print("\nDIAGNÓSTICO OMP:",omp_status)
    if not omp_stats.empty:print(omp_stats.to_string(index=False))
    print("\nResultados:",out)
    print("Figuras publicación:",figout)
    return master,omp_diag,omp_ts,clasif



# =============================================================================
# V3.0 — CONVENCIÓN VERTICAL ÚNICA + DIAGRAMA Θ–SA TRIPANEL
# =============================================================================
# Cambios principales respecto de V2.0:
#   1) Todas las figuras con dimensión vertical usan z = -|profundidad| [m].
#      Superficie = 0 m; profundidades mayores = valores más negativos y se
#      representan más abajo en la figura.
#   2) El diagrama Θ–SA principal se presenta en tres paneles:
#         a) color continuo = profundidad positiva [m];
#         b) color continuo = O2 disuelto [µmol kg-1];
#         c) color categórico = nivel hidrográfico identificado.
#      Se usa un único marcador circular para los datos; la masa de agua se
#      diferencia exclusivamente por color en el panel (c).
#   3) Los SWT mSUW, mAAIW y mNADW se muestran como estrellas de referencia,
#      sin confundirse con los cinco niveles hidrográficos observables.

VERSION = "3.0"
CFG.update({
    "carpeta_salida": DEFAULT_SESSION / "ANALISIS_MASAS_AGUA_V3_0",
    "font_base": 13,
    "font_axes": 14,
    "font_title": 16,
    "font_suptitle": 18,
    "font_legend": 11.5,
    "ts_marker_size": 28,
    "ts_marker_alpha": 0.82,
    "ts_marker_edge_lw": 0.35,
})


def to_negative_z(values):
    """Convención gráfica única: z=0 en superficie y z<0 hacia el fondo."""
    s = pd.to_numeric(pd.Series(values), errors="coerce")
    return -np.abs(s.to_numpy(float))


def _set_vertical_axis_negative(ax, depths_positive, pad_fraction=0.03):
    """Fija 0 m arriba y profundidades mayores (más negativas) abajo."""
    z = to_negative_z(depths_positive)
    z = z[np.isfinite(z)]
    if len(z) == 0:
        return
    zmin = float(np.nanmin(z))
    pad = max(1.0, abs(zmin) * pad_fraction)
    ax.set_ylim(zmin - pad, 0.0)


def fig_ts_tripanel_v3(base, nucleos, omp_diag, o2col, profile_col, out):
    """Figura principal Θ–SA: profundidad, O2 y clasificación hidrográfica."""
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D

    x = preparar_dataset_masas_v2(base, nucleos, omp_diag, profile_col, o2col)
    if x.empty:
        return

    # Límites robustos controlados por hidrografía PASS.
    xp = filtro_hidrografico_pass(base)
    if xp.empty:
        xp = filtro_hidrografico(base)
    xlim = _robust_limits(xp["salinidad_absoluta_g_kg"])
    ylim = _robust_limits(xp["temperatura_conservativa_C"])
    if xlim is None or ylim is None:
        return
    xlim = (
        min(xlim[0], float(SWT.loc["SA_g_kg"].min()) - 0.15),
        max(xlim[1], float(SWT.loc["SA_g_kg"].max()) + 0.15),
    )
    ylim = (
        min(ylim[0], float(SWT.loc["CT_C"].min()) - 0.4),
        max(ylim[1], float(SWT.loc["CT_C"].max()) + 0.4),
    )

    # Panel (a): profundidad positiva como trazador continuo.
    depth = pd.to_numeric(x["profundidad_m"], errors="coerce")
    dvalid = depth[np.isfinite(depth)]
    if len(dvalid) == 0:
        return
    dnorm = Normalize(vmin=float(dvalid.min()), vmax=float(dvalid.max()))
    cmap_depth = plt.cm.viridis

    # Panel (b): O2 robusto (P01–P99) para evitar que pocos outliers dominen.
    if o2col and o2col in x:
        xo = x[x[o2col].notna()].copy()
    else:
        xo = x.iloc[0:0].copy()
    onorm = None
    cmap_o2 = plt.cm.plasma
    if not xo.empty:
        ov = pd.to_numeric(xo[o2col], errors="coerce").dropna().to_numpy(float)
        if len(ov):
            if len(ov) >= 20:
                olo, ohi = np.quantile(ov, [0.01, 0.99])
            else:
                olo, ohi = np.nanmin(ov), np.nanmax(ov)
            if np.isfinite(olo) and np.isfinite(ohi) and ohi > olo:
                onorm = Normalize(vmin=float(olo), vmax=float(ohi))

    # Panel (c): colores categóricos fijos para interpretación hidrográfica.
    mass_colors = {
        "CSW": "#377eb8",
        "SUW": "#e41a1c",
        "Aguas centrales": "#ff7f00",
        "AAIW": "#4daf4a",
        "NADW": "#984ea3",
        "Sin clasificar": "#bdbdbd",
    }
    mass_order = ["CSW", "SUW", "Aguas centrales", "AAIW", "NADW", "Sin clasificar"]

    fig, axs = plt.subplots(
        1, 3,
        figsize=(22.0, 7.8),
        sharex=True,
        sharey=True,
        layout="constrained",
    )

    # Base común.
    for ax in axs:
        add_sigma0_fixed(ax, xlim, ylim)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("Salinidad absoluta, SA [g kg⁻¹]")
        ax.grid(alpha=0.14)
    axs[0].set_ylabel("Temperatura conservativa, Θ [°C]")

    # ------------------------------------------------------------------
    # (a) Profundidad
    # ------------------------------------------------------------------
    sc1 = axs[0].scatter(
        x["salinidad_absoluta_g_kg"],
        x["temperatura_conservativa_C"],
        c=x["profundidad_m"],
        cmap=cmap_depth,
        norm=dnorm,
        marker="o",
        s=CFG["ts_marker_size"],
        alpha=CFG["ts_marker_alpha"],
        edgecolors="k",
        linewidths=CFG["ts_marker_edge_lw"],
        rasterized=True,
        zorder=3,
    )
    cb1 = fig.colorbar(sc1, ax=axs[0], pad=0.018)
    cb1.set_label("Profundidad [m]", fontsize=CFG["font_axes"])
    cb1.ax.tick_params(labelsize=CFG["font_base"])
    axs[0].set_title("a) Estructura Θ–SA y profundidad")

    # ------------------------------------------------------------------
    # (b) Oxígeno disuelto
    # ------------------------------------------------------------------
    if not xo.empty and onorm is not None:
        sc2 = axs[1].scatter(
            xo["salinidad_absoluta_g_kg"],
            xo["temperatura_conservativa_C"],
            c=xo[o2col],
            cmap=cmap_o2,
            norm=onorm,
            marker="o",
            s=CFG["ts_marker_size"],
            alpha=CFG["ts_marker_alpha"],
            edgecolors="k",
            linewidths=CFG["ts_marker_edge_lw"],
            rasterized=True,
            zorder=3,
        )
        cb2 = fig.colorbar(sc2, ax=axs[1], pad=0.018)
        cb2.set_label("Oxígeno disuelto [µmol kg⁻¹]", fontsize=CFG["font_axes"])
        cb2.ax.tick_params(labelsize=CFG["font_base"])
    else:
        axs[1].text(
            0.5, 0.5,
            "Oxígeno no disponible o no válido",
            transform=axs[1].transAxes,
            ha="center", va="center",
            fontsize=CFG["font_title"],
        )
    axs[1].set_title("b) Estructura Θ–SA y oxígeno disuelto")

    # ------------------------------------------------------------------
    # (c) Masas de agua: un mismo marcador, color categórico.
    # ------------------------------------------------------------------
    for mass in mass_order:
        g = x[x["masa_hidrografica_v2"].eq(mass)]
        if g.empty:
            continue
        alpha = 0.30 if mass == "Sin clasificar" else 0.88
        size = 16 if mass == "Sin clasificar" else CFG["ts_marker_size"]
        axs[2].scatter(
            g["salinidad_absoluta_g_kg"],
            g["temperatura_conservativa_C"],
            color=mass_colors[mass],
            marker="o",
            s=size,
            alpha=alpha,
            edgecolors="k" if mass != "Sin clasificar" else "none",
            linewidths=CFG["ts_marker_edge_lw"] if mass != "Sin clasificar" else 0,
            rasterized=True,
            zorder=4 if mass != "Sin clasificar" else 2,
        )
    axs[2].set_title("c) Estructura Θ–SA y masas de agua")

    # SWT: referencia matemática del OMP. Se mantienen como estrellas blancas
    # y se etiquetan para distinguirlos de la clasificación de cinco niveles.
    for ax in axs:
        for name in SWT.columns:
            ax.scatter(
                SWT.loc["SA_g_kg", name],
                SWT.loc["CT_C", name],
                marker="*", s=200,
                facecolor="white", edgecolor="k", linewidth=1.1,
                zorder=8,
            )
            ax.annotate(
                name,
                (SWT.loc["SA_g_kg", name], SWT.loc["CT_C", name]),
                xytext=(6, 6), textcoords="offset points",
                fontsize=11.5, weight="bold",
            )

    # Leyenda categórica únicamente en el tercer panel.
    handles = []
    for mass in mass_order:
        if mass == "Sin clasificar" and not x["masa_hidrografica_v2"].eq(mass).any():
            continue
        handles.append(Line2D(
            [0], [0], marker="o", linestyle="None",
            markerfacecolor=mass_colors[mass], markeredgecolor="k" if mass != "Sin clasificar" else mass_colors[mass],
            markersize=8.5, label=mass,
        ))
    handles.append(Line2D(
        [0], [0], marker="*", linestyle="None",
        markerfacecolor="white", markeredgecolor="k", color="k",
        markersize=11, label="SWT OMP",
    ))
    axs[2].legend(
        handles=handles,
        loc="best",
        frameon=True,
        title="Clasificación hidrográfica",
        fontsize=CFG["font_legend"],
        title_fontsize=CFG["font_legend"] + 0.5,
    )

    fig.suptitle(
        "Diagramas Θ–SA conjuntos: profundidad, oxígeno y clasificación hidrográfica",
        fontsize=CFG["font_suptitle"], weight="bold",
    )
    savefig(fig, out / "Fig02_V3_TS_tripanel_profundidad_oxigeno_masas.png")


def fig_perfiles_medianos_v3(base, nucleos, o2col, out):
    """Perfiles medianos con IQR y convención z negativa uniforme."""
    q = filtro_hidrografico_pass(base)
    if q.empty:
        q = filtro_hidrografico(base)
    if q.empty:
        return

    dz = 10.0
    q = q.copy()
    q["zbin_m"] = np.round(pd.to_numeric(q["profundidad_m"], errors="coerce") / dz) * dz
    q["z_plot_m"] = -np.abs(q["zbin_m"])

    vars_ = [
        ("salinidad_absoluta_g_kg", "SA [g kg⁻¹]"),
        ("temperatura_conservativa_C", "Θ [°C]"),
        ("sigma0_kg_m3", "sigma0 [kg m⁻³]"),
    ]
    if o2col and o2col in q:
        vars_.append((o2col, "O₂ [µmol kg⁻¹]"))

    fig, axs = plt.subplots(
        1, len(vars_),
        figsize=(4.7 * len(vars_), 7.4),
        sharey=True,
        layout="constrained",
        squeeze=False,
    )
    axs = axs.ravel()

    stats = {}
    for c, _ in vars_:
        g = q.groupby("zbin_m")[c]
        stats[c] = pd.DataFrame({
            "med": g.median(),
            "q25": g.quantile(0.25),
            "q75": g.quantile(0.75),
        }).dropna(how="all")

    suw_med = -float(pd.to_numeric(nucleos.get("SUW_core_depth_m"), errors="coerce").median()) if "SUW_core_depth_m" in nucleos else np.nan
    if "AAIW_core_truncado_fondo" in nucleos:
        aa_src = nucleos[~nucleos["AAIW_core_truncado_fondo"].fillna(True)]
    else:
        aa_src = nucleos
    aaiw_med = -float(pd.to_numeric(aa_src.get("AAIW_core_depth_m"), errors="coerce").median()) if "AAIW_core_depth_m" in aa_src else np.nan

    zall = []
    for ax, (c, label) in zip(axs, vars_):
        st = stats[c]
        if st.empty:
            ax.set_visible(False)
            continue
        zplot = -np.abs(st.index.to_numpy(float))
        zall.extend(zplot[np.isfinite(zplot)].tolist())
        ax.fill_betweenx(zplot, st["q25"], st["q75"], alpha=0.18, label="IQR")
        ax.plot(st["med"], zplot, lw=2.2, label="Mediana")
        if np.isfinite(suw_med):
            ax.axhline(suw_med, ls="--", lw=1.3, label="SUW mediana")
        if np.isfinite(aaiw_med):
            ax.axhline(aaiw_med, ls=":", lw=1.5, label="AAIW mediana")
        ax.set_xlabel(label)
        ax.set_title(label)
        ax.grid(alpha=0.15)

    axs[0].set_ylabel("Coordenada vertical, z [m]")
    if zall:
        zmin = float(np.nanmin(zall))
        axs[0].set_ylim(zmin - max(10.0, 0.02 * abs(zmin)), 0.0)
    # Leyenda sólo en el último panel visible para reducir ruido.
    visible = [ax for ax in axs if ax.get_visible()]
    if visible:
        visible[-1].legend(loc="best", fontsize=CFG["font_legend"])
    fig.suptitle("Estructura vertical mediana de la campaña (banda = IQR)", weight="bold")
    savefig(fig, out / "Fig03_V3_perfiles_medianos_IQR_zneg.png")


def fig_capas_superficiales_v3(nucleos, out):
    """MLD/ILD en z negativa; BLT permanece como espesor positivo."""
    if nucleos.empty:
        return
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.8), layout="constrained")
    depths=[]
    for periodo, g in nucleos.groupby("periodo"):
        g = g.sort_values("numero_estacion")
        zmld = -np.abs(pd.to_numeric(g["MLD_m"], errors="coerce"))
        zild = -np.abs(pd.to_numeric(g["ILD_m"], errors="coerce"))
        depths.extend((-zmld[np.isfinite(zmld)]).tolist())
        depths.extend((-zild[np.isfinite(zild)]).tolist())
        axs[0].plot(g["numero_estacion"], zmld, marker="o", label=f"MLD {periodo}")
        axs[0].plot(g["numero_estacion"], zild, marker="s", label=f"ILD {periodo}")
        axs[1].plot(g["numero_estacion"], g["BLT_m"], marker="o", label=periodo)
    axs[0].set(xlabel="Número de estación", ylabel="Coordenada vertical, z [m]", title="Profundidad de MLD e ILD")
    if depths:
        _set_vertical_axis_negative(axs[0], depths, pad_fraction=0.08)
    axs[1].set(xlabel="Número de estación", ylabel="BLT = ILD − MLD [m]", title="Espesor de capa barrera")
    for ax in axs:
        ax.grid(alpha=0.18)
        ax.legend(fontsize=CFG["font_legend"])
    fig.suptitle("Estructura de la capa superficial", weight="bold")
    savefig(fig, out / "Fig04_V3_MLD_ILD_BLT_zneg.png")


def fig_nucleos_v3(nucleos, out):
    """Núcleos SUW/AAIW con 0 m en superficie y z más negativa hacia abajo."""
    if nucleos.empty:
        return
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.8), layout="constrained", sharex=True)
    suw_depth=[]; aa_depth=[]
    for periodo, g in nucleos.groupby("periodo"):
        g = g.sort_values("numero_estacion")
        zsuw = -np.abs(pd.to_numeric(g["SUW_core_depth_m"], errors="coerce"))
        suw_depth.extend((-zsuw[np.isfinite(zsuw)]).tolist())
        axs[0].plot(g["numero_estacion"], zsuw, marker="o", label=periodo)

        good = g[~g["AAIW_core_truncado_fondo"].fillna(True)]
        bad = g[g["AAIW_core_truncado_fondo"].fillna(True)]
        zgood = -np.abs(pd.to_numeric(good["AAIW_core_depth_m"], errors="coerce"))
        zbad = -np.abs(pd.to_numeric(bad["AAIW_core_depth_m"], errors="coerce"))
        aa_depth.extend((-zgood[np.isfinite(zgood)]).tolist())
        aa_depth.extend((-zbad[np.isfinite(zbad)]).tolist())
        if len(good):
            axs[1].plot(good["numero_estacion"], zgood, marker="s", label=periodo)
        if len(bad):
            axs[1].scatter(bad["numero_estacion"], zbad, marker="x", s=70, label=f"{periodo} truncado")

    axs[0].set(title="Núcleo SUW (máximo de SA)", ylabel="Coordenada vertical, z [m]", xlabel="Número de estación")
    axs[1].set(title="Núcleo AAIW (mínimo de SA)", ylabel="Coordenada vertical, z [m]", xlabel="Número de estación")
    if suw_depth:
        _set_vertical_axis_negative(axs[0], suw_depth, pad_fraction=0.08)
    if aa_depth:
        _set_vertical_axis_negative(axs[1], aa_depth, pad_fraction=0.08)
    for ax in axs:
        ax.grid(alpha=0.18)
        ax.legend(fontsize=CFG["font_legend"])
    fig.suptitle("Profundidad de núcleos hidrográficos", weight="bold")
    savefig(fig, out / "Fig05_V3_nucleos_SUW_AAIW_zneg.png")


def fig_secciones_hidrograficas_v3(base, o2col, out):
    """Secciones con z negativa: superficie 0 m, fondo hacia valores negativos."""
    hydro = filtro_hidrografico_pass(base)
    oxy = filtro_o2_pass(base, o2col) if o2col else hydro.iloc[0:0]
    variables = [
        ("salinidad_absoluta_g_kg", "SA [g kg⁻¹]", hydro),
        ("temperatura_conservativa_C", "Θ [°C]", hydro),
        ("sigma0_kg_m3", "sigma0 [kg m⁻³]", hydro),
    ]
    if o2col:
        variables.append((o2col, "O₂ [µmol kg⁻¹]", oxy))
    periods = [p for p in ("Diurno", "Nocturno") if (base["periodo"] == p).any()]

    for tname, stations in CFG["transectos"].items():
        present = set(pd.to_numeric(base["numero_estacion"], errors="coerce").dropna().astype(int))
        use = [s for s in stations if s in present]
        if len(use) < 2:
            continue
        fig, axs = plt.subplots(
            len(periods), len(variables),
            figsize=(4.6 * len(variables), 4.6 * len(periods)),
            layout="constrained", sharey="row", squeeze=False,
        )
        for i, periodo in enumerate(periods):
            row_zmins=[]
            for j, (col, label, source) in enumerate(variables):
                ax = axs[i, j]
                sec = _section_matrix_stations(source, col, periodo, use)
                if sec is None:
                    ax.set_visible(False)
                    continue
                xd, zpos, M, labs, _ = sec
                zplot = -np.abs(zpos)
                row_zmins.append(float(np.nanmin(zplot)))
                pc = ax.pcolormesh(xd, zplot, M, shading="nearest")
                cb = fig.colorbar(pc, ax=ax, pad=0.02)
                cb.set_label(label, fontsize=CFG["font_base"])
                cb.ax.tick_params(labelsize=CFG["font_base"]-1)
                ax.set_xticks(xd)
                ax.set_xticklabels(labs)
                ax.grid(alpha=0.10)
                if i == len(periods)-1:
                    ax.set_xlabel("Distancia dentro del grupo [km]")
                if j == 0:
                    ax.set_ylabel(f"{periodo}\nCoordenada vertical, z [m]")
                if i == 0:
                    ax.set_title(label)
            if row_zmins:
                zmin=min(row_zmins)
                for j in range(len(variables)):
                    if axs[i,j].get_visible():
                        axs[i,j].set_ylim(zmin, 0.0)
        fig.suptitle(
            f"Sección hidrográfica · Grupo {tname} · estaciones {', '.join(f'E{s:02d}' for s in use)}\n"
            "Agrupación configurable: confirmar con la bitácora",
            fontsize=CFG["font_title"],
        )
        savefig(fig, out / f"Fig06_V3_seccion_{slug(tname)}_zneg.png")


def fig_omp_sensibilidad_v3(omp_full, omp_ts, out):
    if omp_full.empty or omp_ts.empty:
        return
    dz=25.0
    a=omp_full[omp_full.omp_ok].copy(); b=omp_ts[omp_ts.omp_ok].copy()
    a["zbin"]=np.round(a.profundidad_m/dz)*dz; b["zbin"]=np.round(b.profundidad_m/dz)*dz
    ma=a.groupby("zbin")[["frac_mSUW","frac_mAAIW","frac_mNADW"]].mean()
    mb=b.groupby("zbin")[["frac_mSUW","frac_mAAIW","frac_mNADW"]].mean()
    fig,axs=plt.subplots(1,3,figsize=(14.5,6.8),layout="constrained",sharey=True)
    zall=[]
    for ax,mass in zip(axs,["mSUW","mAAIW","mNADW"]):
        c=f"frac_{mass}"
        if c in ma.columns:
            z=-np.abs(ma.index.to_numpy(float)); zall.extend(z.tolist())
            ax.plot(100*ma[c],z,label="Θ+SA+O₂")
        if c in mb.columns:
            z=-np.abs(mb.index.to_numpy(float)); zall.extend(z.tolist())
            ax.plot(100*mb[c],z,"--",label="Θ+SA (sensibilidad)")
        ax.axvline(30,ls=":",lw=1.2)
        ax.set(xlabel="Fracción media [%]",title=mass)
        ax.grid(alpha=.18); ax.legend(fontsize=CFG["font_legend"])
    axs[0].set_ylabel("Coordenada vertical, z [m]")
    if zall:
        axs[0].set_ylim(float(np.nanmin(zall))-25,0)
    fig.suptitle("Sensibilidad de las fracciones OMP al uso de O₂\nPromedio muestral; no equivale a fracción volumétrica del Caribe",fontsize=CFG["font_title"])
    savefig(fig,out/"Fig07_V3_OMP_sensibilidad_vertical_zneg.png")


def fig_omp_residuales_v3(omp_diag,out):
    if omp_diag.empty:
        return
    cols=[("resid_CT_norm_C","Residual Θ [°C]"),("resid_SA_norm_g_kg","Residual SA [g kg⁻¹]"),("resid_O2_norm_umol_kg","Residual O₂ [µmol kg⁻¹]")]
    fig,axs=plt.subplots(1,3,figsize=(14.5,6.8),layout="constrained",sharey=True)
    dz=25.0
    q=omp_diag.copy(); q["zbin"]=np.round(q.profundidad_m/dz)*dz
    zall=[]
    for ax,(c,label) in zip(axs,cols):
        if c not in q:
            ax.set_visible(False); continue
        g=q.groupby("zbin")[c]
        med=g.median();q25=g.quantile(.25);q75=g.quantile(.75)
        z=-np.abs(med.index.to_numpy(float));zall.extend(z.tolist())
        ax.fill_betweenx(z,q25,q75,alpha=.2,label="IQR")
        ax.plot(med,z,lw=1.8,label="Mediana")
        ax.axvline(0,ls="--",lw=1.1)
        ax.set(xlabel=label,title=label);ax.grid(alpha=.18);ax.legend(fontsize=CFG["font_legend"])
    axs[0].set_ylabel("Coordenada vertical, z [m]")
    if zall:
        axs[0].set_ylim(float(np.nanmin(zall))-25,0)
    fig.suptitle("Diagnóstico de residuales del OMP completo (fracciones normalizadas)",fontsize=CFG["font_title"])
    savefig(fig,out/"Fig08_V3_OMP_residuales_zneg.png")


def fig_omp_estandar_v3(omp_diag,omp_ts,out):
    if omp_diag.empty or omp_ts.empty:
        return
    full=_interp_omp_standard(omp_diag,["frac_mSUW","frac_mAAIW","frac_mNADW"])
    ts=_interp_omp_standard(omp_ts,["frac_mSUW","frac_mAAIW","frac_mNADW"])
    if full.empty or ts.empty:
        return
    rows=[]
    for z in CFG["profundidades_estandar_m"]:
        rf=full[full.prof_objetivo_m.eq(z)];rt=ts[ts.prof_objetivo_m.eq(z)]
        for mass in ["mSUW","mAAIW","mNADW"]:
            rows.append({"profundidad_m":z,"masa":mass,
                         "full_media":100*rf[f"frac_{mass}"].mean() if len(rf) else np.nan,
                         "ts_media":100*rt[f"frac_{mass}"].mean() if len(rt) else np.nan})
    d=pd.DataFrame(rows)
    fig,axs=plt.subplots(1,3,figsize=(14.5,5.8),layout="constrained",sharey=True)
    zall=[]
    for ax,mass in zip(axs,["mSUW","mAAIW","mNADW"]):
        g=d[d.masa.eq(mass)].sort_values("profundidad_m")
        z=-np.abs(g.profundidad_m.to_numpy(float));zall.extend(z.tolist())
        ax.plot(g.full_media,z,marker="o",label="Θ+SA+O₂")
        ax.plot(g.ts_media,z,marker="s",ls="--",label="Θ+SA")
        ax.set(xlabel="Fracción media [%]",title=mass);ax.grid(alpha=.18);ax.legend(fontsize=CFG["font_legend"])
    axs[0].set_ylabel("Coordenada vertical, z [m]")
    if zall:
        axs[0].set_ylim(float(np.nanmin(zall))-50,0)
    fig.suptitle("Comparación OMP en profundidades estándar",fontsize=CFG["font_title"])
    savefig(fig,out/"Fig09_V3_OMP_profundidades_estandar_zneg.png")


def write_methodology_v3(out:Path,path_input:Path,o2col: str|None,omp_status: str):
    text=f"""# Metodología — masas de agua del Caribe colombiano — V3.0

Generado: {datetime.now().isoformat(timespec='seconds')}
Entrada: {path_input}

## Síntesis conceptual

La interpretación hidrográfica se organiza en cinco niveles observables: CSW, SUW, aguas centrales, AAIW y NADW. Para OMP se emplean únicamente mSUW, mAAIW y mNADW de Torres et al. (2023). Las aguas centrales se interpretan como región de transición/mezcla entre SUW y AAIW y no se fuerzan como cuarto end-member.

## Convención vertical de figuras

Todas las figuras que contienen un eje vertical utilizan z = -|profundidad| [m]. La superficie se representa en 0 m y las profundidades mayores se muestran como valores más negativos hacia la parte inferior de la figura. La profundidad usada como barra de color en el diagrama Θ–SA permanece positiva [m], porque funciona allí como trazador continuo y no como eje vertical.

## Figura Θ–SA principal

Se utilizan tres paneles: (a) color por profundidad, (b) color por oxígeno disuelto y (c) color categórico por masa de agua identificada. Los datos usan un marcador circular uniforme; en el tercer panel el color distingue CSW, SUW, aguas centrales, AAIW y NADW. Los SWT mSUW, mAAIW y mNADW se muestran como estrellas de referencia del OMP.

## Criterios principales

- MLD: delta sigma0 = 0.03 kg m-3 respecto a 10 m.
- ILD: |delta temperatura conservativa| = 0.2 °C respecto a 10 m.
- SUW: máximo de SA por debajo de {CFG['suw_min_depth_m']:.0f} m.
- AAIW: mínimo de SA profundo, con control operacional >= {CFG['aaiw_min_depth_m']:.0f} m y evaluación de truncamiento por fondo.
- Aguas centrales: intervalo hidrográfico entre los núcleos SUW y AAIW; no se asigna un end-member independiente.
- NADW: se reporta como no alcanzada, influencia, no resuelta o resuelta combinando profundidad alcanzada, recuperación de SA/O2 y soporte OMP.
- OMP: mSUW + mAAIW + mNADW debajo del máximo de salinidad; residual de masa aceptable <= {CFG['omp_mass_residual_max_pct']:.0f} %.
- Fracción >= {100*CFG['fraccion_significativa']:.0f} % se usa como umbral descriptivo de contribución importante.

## SWT Torres et al. (2023)

{SWT.to_string()}

Pesos OMP: {OMP_WEIGHTS}

Columna de O2: {o2col or 'NO DISPONIBLE'}
Estado diagnóstico OMP: {omp_status}
"""
    (out/"METODOLOGIA_MASAS_AGUA_V3.md").write_text(text,encoding="utf-8")


def write_captions_v3(out:Path):
    txt="""# Leyendas sugeridas — figuras V3.0

**Figura 1.** Localización de las estaciones CTD y profundidad máxima muestreada. Las agrupaciones espaciales configuradas deben verificarse con la bitácora antes de denominarse transectos definitivos.

**Figura 2.** Diagramas conjuntos de temperatura conservativa (Θ) y salinidad absoluta (SA). (a) Observaciones coloreadas por profundidad; (b) observaciones coloreadas por oxígeno disuelto; (c) observaciones coloreadas según la interpretación hidrográfica en CSW, SUW, aguas centrales, AAIW y NADW. Los puntos utilizan un marcador uniforme para evitar superposición de códigos visuales. Las estrellas corresponden a los source water types mSUW, mAAIW y mNADW empleados en el OMP según Torres et al. (2023). Las isolíneas corresponden a sigma0.

**Figura 3.** Perfiles verticales medianos de SA, Θ, sigma0 y O2 para datos QARTOD PASS; la banda corresponde al rango intercuartílico. La coordenada vertical se expresa como z negativa: 0 m en la superficie y valores progresivamente más negativos hacia el fondo.

**Figura 4.** Profundidad de MLD e ILD y espesor de capa barrera. MLD e ILD se representan en coordenada vertical z negativa; BLT se mantiene como espesor positivo.

**Figura 5.** Profundidad de los núcleos SUW y AAIW. La coordenada vertical es z negativa, con la superficie en 0 m. Los núcleos AAIW truncados por el fondo se identifican separadamente.

**Figura 6.** Secciones hidrográficas por agrupación espacial y periodo. Todas las secciones utilizan z negativa, con 0 m en superficie y profundidades mayores hacia valores más negativos.

**Figura 7.** Sensibilidad de las fracciones OMP al uso de O2. El eje vertical usa z negativa. La solución Θ+SA se interpreta como análisis de sensibilidad y no como sustituto automático del OMP completo.

**Figura 8.** Residuales del OMP completo en función de z negativa. La línea representa la mediana y la banda el rango intercuartílico.

**Figura 9.** Comparación de fracciones OMP en profundidades estándar representadas como z negativa.
"""
    (out/"FIGURAS_V3_LEYENDAS.md").write_text(txt,encoding="utf-8")


def main():
    apply_publication_style()
    print("="*84)
    print("ANÁLISIS DE MASAS DE AGUA — V3.0 PARA PUBLICACIÓN")
    print("5 niveles hidrográficos + OMP + Θ–SA tripanel + convención vertical z negativa")
    print("="*84)

    path=choose_input_file(Path(CFG["archivo_entrada"]))
    out=Path(CFG["carpeta_salida"]);out.mkdir(parents=True,exist_ok=True)
    figout=out/"figuras_publicacion";figout.mkdir(exist_ok=True)
    base,decisiones,pcol,o2col,profile_col=cargar_dataset(path)
    base["z_m"]=-np.abs(pd.to_numeric(base["profundidad_m"],errors="coerce"))
    print(f"Perfiles: {base[profile_col].nunique()} | bins: {len(base):,} | O2: {o2col}")

    nuc=[]
    for pid,g in base.groupby(profile_col,sort=False):
        nuc.append(analizar_perfil_hidrografico_v2(g,pcol,o2col,profile_id=pid))
    nucleos=pd.DataFrame(nuc)

    omp,excl=ejecutar_omp(base,nucleos,o2col,profile_col)
    omp_ts=ejecutar_omp_ts(base,nucleos,profile_col)
    omp_diag,omp_stats=augmentar_diagnostico_omp(omp,omp_ts)
    if not omp_diag.empty:
        omp_diag["z_m"]=-np.abs(pd.to_numeric(omp_diag["profundidad_m"],errors="coerce"))
    if not omp_ts.empty:
        omp_ts["z_m"]=-np.abs(pd.to_numeric(omp_ts["profundidad_m"],errors="coerce"))

    omp_status="NO_EVALUADO"
    if not omp_stats.empty and "resid_O2_norm_umol_kg" in set(omp_stats.variable):
        rr=omp_stats[omp_stats.variable.eq("resid_O2_norm_umol_kg")].iloc[0]
        if abs(rr.mediana)>=CFG["omp_o2_bias_warn_umol_kg"] and rr.fraccion_mismo_signo>=CFG["omp_same_sign_warn_fraction"]:
            omp_status="EXPLORATORIO_DESACOPLE_O2_SISTEMATICO"
        else:
            omp_status="SIN_ALARMA_SISTEMATICA_O2"

    nucleos=enriquecer_nucleos_v2(base,nucleos,omp_diag,o2col,profile_col)
    omp_res=resumen_omp_por_perfil(omp_diag) if not omp_diag.empty else pd.DataFrame()
    omp_std=omp_en_profundidades_estandar(omp_diag) if not omp_diag.empty else pd.DataFrame()
    ts_std=_interp_omp_standard(omp_ts,["frac_mSUW","frac_mAAIW","frac_mNADW"])
    hydro_summary=resumen_hidrografico_robusto_v2(nucleos)
    master=tabla_maestra_v2(nucleos,omp_res)

    sens=[]
    if not omp_std.empty or not ts_std.empty:
        for z in CFG["profundidades_estandar_m"]:
            a=omp_std[omp_std.prof_objetivo_m.eq(z)] if not omp_std.empty else pd.DataFrame()
            b=ts_std[ts_std.prof_objetivo_m.eq(z)] if not ts_std.empty else pd.DataFrame()
            for mass in ["mSUW","mAAIW","mNADW"]:
                sens.append({"profundidad_m":z,"z_m":-abs(float(z)),"masa":mass,
                             "OMP_full_media":float(a[f"frac_{mass}"].mean()) if len(a) else np.nan,
                             "OMP_TS_media":float(b[f"frac_{mass}"].mean()) if len(b) else np.nan,
                             "diferencia_full_menos_TS":float(a[f"frac_{mass}"].mean()-b[f"frac_{mass}"].mean()) if len(a) and len(b) else np.nan})
    sens=pd.DataFrame(sens)

    clasif=preparar_dataset_masas_v2(base,nucleos,omp_diag,profile_col,o2col)
    clasif["z_m"]=-np.abs(pd.to_numeric(clasif["profundidad_m"],errors="coerce"))

    # CSV reproducibles V3.
    master.to_csv(out/"MASAS_AGUA_POR_PERFIL_V3.csv",index=False,encoding="utf-8-sig")
    hydro_summary.to_csv(out/"RESUMEN_HIDROGRAFICO_ROBUSTO_V3.csv",index=False,encoding="utf-8-sig")
    clasif.to_csv(out/"CLASIFICACION_HIDROGRAFICA_POR_BIN_V3.csv",index=False,encoding="utf-8-sig")
    if not omp_diag.empty:omp_diag.to_csv(out/"OMP_FULL_DIAGNOSTICO_V3.csv",index=False,encoding="utf-8-sig")
    if not omp_ts.empty:omp_ts.to_csv(out/"OMP_TS_SENSIBILIDAD_V3.csv",index=False,encoding="utf-8-sig")
    if not omp_stats.empty:omp_stats.to_csv(out/"OMP_RESIDUALES_RESUMEN_V3.csv",index=False,encoding="utf-8-sig")
    if not sens.empty:sens.to_csv(out/"OMP_SENSIBILIDAD_PROF_ESTANDAR_V3.csv",index=False,encoding="utf-8-sig")
    if not excl.empty:excl.to_csv(out/"OMP_PERFILES_EXCLUIDOS_V3.csv",index=False,encoding="utf-8-sig")

    with pd.ExcelWriter(out/"RESULTADOS_MASAS_AGUA_V3.xlsx",engine="openpyxl") as wr:
        master.to_excel(wr,sheet_name="Masas_por_perfil",index=False)
        hydro_summary.to_excel(wr,sheet_name="Resumen_robusto",index=False)
        clasif.to_excel(wr,sheet_name="Clasificacion_bins",index=False)
        if not omp_diag.empty:omp_diag.to_excel(wr,sheet_name="OMP_full",index=False)
        if not omp_ts.empty:omp_ts.to_excel(wr,sheet_name="OMP_TS_sens",index=False)
        if not omp_stats.empty:omp_stats.to_excel(wr,sheet_name="OMP_residuales",index=False)
        if not sens.empty:sens.to_excel(wr,sheet_name="OMP_prof_estandar",index=False)
        if not excl.empty:excl.to_excel(wr,sheet_name="OMP_excluidos",index=False)
        SWT.to_excel(wr,sheet_name="SWT_Torres2023")
        pd.DataFrame([CFG]).astype(str).to_excel(wr,sheet_name="Configuracion",index=False)

    # Figuras V3: convención vertical uniforme y Θ–SA tripanel.
    fig_mapa_estaciones_vf(base,figout)
    fig_ts_tripanel_v3(base,nucleos,omp_diag,o2col,profile_col,figout)
    fig_perfiles_medianos_v3(base,nucleos,o2col,figout)
    fig_capas_superficiales_v3(nucleos,figout)
    fig_nucleos_v3(nucleos,figout)
    fig_secciones_hidrograficas_v3(base,o2col,figout)
    fig_omp_sensibilidad_v3(omp_diag,omp_ts,figout)
    fig_omp_residuales_v3(omp_diag,figout)
    fig_omp_estandar_v3(omp_diag,omp_ts,figout)

    write_methodology_v3(out,path,o2col,omp_status)
    write_captions_v3(out)

    meta={"version":VERSION,"input":str(path),"output":str(out),"source_reference":SOURCE_REFERENCE,
          "SWT":SWT.to_dict(),"OMP_weights":OMP_WEIGHTS,"OMP_status":omp_status,
          "sintesis_hidrografica":["CSW","SUW","Aguas centrales","AAIW","NADW"],
          "endmembers_OMP":["mSUW","mAAIW","mNADW"],
          "convencion_vertical":"z=-|profundidad| m; superficie=0; fondo=hacia valores más negativos",
          "figura_TS":"tripanel: profundidad / O2 / masa identificada por color",
          "nota":"Cinco niveles para interpretación hidrográfica; tres SWT independientes para OMP. Aguas centrales no se fuerza como end-member.",
          "transectos_configurados":CFG["transectos"]}
    (out/"metadata_masas_agua_v3.json").write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding="utf-8")

    print("\nRESUMEN HIDROGRÁFICO ROBUSTO")
    print(hydro_summary.to_string(index=False))
    print("\nDIAGNÓSTICO OMP:",omp_status)
    if not omp_stats.empty:print(omp_stats.to_string(index=False))
    print("\nResultados:",out)
    print("Figuras publicación:",figout)
    return master,omp_diag,omp_ts,clasif


# __main__ block moved to end by V3.4 patch.


# =====================================================================
# V3.4 PATCH — O2 en mg/L para figuras, TS tripanel más limpio,
# secciones hidrográficas continuas y secciones OMP por fracciones.
# =====================================================================
VERSION = "3.4"
CFG["carpeta_salida"] = "ANALISIS_MASAS_AGUA_V3_4"
O2_UMOLKG_TO_MGL = 0.0328  # aproximación agua de mar ~1.025 kg L-1


def _add_o2_mgl_column(df: pd.DataFrame, o2col: str | None, outcol: str = "O2_mg_L") -> pd.DataFrame:
    x = df.copy()
    if o2col and o2col in x.columns:
        x[outcol] = pd.to_numeric(x[o2col], errors="coerce") * O2_UMOLKG_TO_MGL
    return x


def _robust_vrange(series, qlo=0.02, qhi=0.98, pad_frac=0.02):
    s = pd.to_numeric(series, errors="coerce").to_numpy(float)
    s = s[np.isfinite(s)]
    if len(s) == 0:
        return (0.0, 1.0)
    if len(s) >= 20:
        lo, hi = np.quantile(s, [qlo, qhi])
    else:
        lo, hi = np.nanmin(s), np.nanmax(s)
    if not np.isfinite(lo) or not np.isfinite(hi):
        return (0.0, 1.0)
    if hi <= lo:
        hi = lo + 1e-6
    pad = pad_frac * (hi - lo)
    return (float(lo - pad), float(hi + pad))


def fig_ts_tripanel_v31(base, nucleos, omp_diag, o2col, profile_col, out):
    """Tripanel Θ–SA: z negativa, O2 en mg/L y clasificación hidrográfica."""
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D

    x = preparar_dataset_masas_v2(base, nucleos, omp_diag, profile_col, o2col)
    if x.empty:
        return
    x = _add_o2_mgl_column(x, o2col)
    x['z_m'] = -np.abs(pd.to_numeric(x['profundidad_m'], errors='coerce'))

    xp = filtro_hidrografico_pass(base)
    if xp.empty:
        xp = filtro_hidrografico(base)
    xlim = _robust_limits(xp['salinidad_absoluta_g_kg'])
    ylim = _robust_limits(xp['temperatura_conservativa_C'])
    if xlim is None or ylim is None:
        return
    xlim = (min(xlim[0], float(SWT.loc['SA_g_kg'].min()) - 0.18),
            max(xlim[1], float(SWT.loc['SA_g_kg'].max()) + 0.22))
    ylim = (min(ylim[0], float(SWT.loc['CT_C'].min()) - 0.45),
            max(ylim[1], float(SWT.loc['CT_C'].max()) + 0.45))

    zvals = pd.to_numeric(x['z_m'], errors='coerce').to_numpy(float)
    zvals = zvals[np.isfinite(zvals)]
    if len(zvals) == 0:
        return
    zmin = float(np.nanmin(zvals))
    znorm = Normalize(vmin=zmin, vmax=0.0)

    xo = x[x['O2_mg_L'].notna()].copy() if 'O2_mg_L' in x else x.iloc[0:0].copy()
    onorm = None
    if not xo.empty:
        ov = pd.to_numeric(xo['O2_mg_L'], errors='coerce').dropna().to_numpy(float)
        if len(ov):
            olo, ohi = np.quantile(ov, [0.02, 0.98]) if len(ov) >= 20 else (float(np.nanmin(ov)), float(np.nanmax(ov)))
            if np.isfinite(olo) and np.isfinite(ohi) and ohi > olo:
                onorm = Normalize(vmin=float(olo), vmax=float(ohi))

    mass_colors = {
        'CSW': '#1f78b4',
        'SUW': '#e31a1c',
        'Aguas centrales': '#ff7f00',
        'AAIW': '#33a02c',
        'NADW': '#6a3d9a',
        'Sin clasificar': '#bdbdbd',
    }
    mass_order = ['CSW', 'SUW', 'Aguas centrales', 'AAIW', 'NADW', 'Sin clasificar']

    fig, axs = plt.subplots(1, 3, figsize=(22.2, 7.6), sharex=True, sharey=True, layout='constrained')
    for ax in axs:
        add_sigma0_fixed(ax, xlim, ylim)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel('Salinidad absoluta, SA [g kg⁻¹]')
        ax.grid(alpha=0.12)
    axs[0].set_ylabel('Temperatura conservativa, Θ [°C]')

    sc1 = axs[0].scatter(
        x['salinidad_absoluta_g_kg'], x['temperatura_conservativa_C'],
        c=x['z_m'], cmap=plt.cm.viridis, norm=znorm,
        marker='o', s=18, alpha=0.72, edgecolors='none', linewidths=0.0,
        rasterized=True, zorder=3)
    cb1 = fig.colorbar(sc1, ax=axs[0], pad=0.018)
    cb1.set_label('Coordenada vertical, z [m]', fontsize=CFG['font_axes'])
    # Superficie arriba (0 m), fondo abajo (valores negativos más intensos)
    ticks = np.linspace(zmin, 0.0, 6)
    cb1.set_ticks(ticks)
    cb1.ax.tick_params(labelsize=CFG['font_base'])
    axs[0].set_title('a) Estructura Θ–SA y profundidad')

    if not xo.empty and onorm is not None:
        sc2 = axs[1].scatter(
            xo['salinidad_absoluta_g_kg'], xo['temperatura_conservativa_C'],
            c=xo['O2_mg_L'], cmap=plt.cm.plasma, norm=onorm,
            marker='o', s=18, alpha=0.72, edgecolors='none', linewidths=0.0,
            rasterized=True, zorder=3)
        cb2 = fig.colorbar(sc2, ax=axs[1], pad=0.018)
        cb2.set_label('Oxígeno disuelto [mg L⁻¹]', fontsize=CFG['font_axes'])
        cb2.ax.tick_params(labelsize=CFG['font_base'])
    else:
        axs[1].text(0.5, 0.5, 'Oxígeno no disponible o no válido',
                    transform=axs[1].transAxes, ha='center', va='center',
                    fontsize=CFG['font_title'])
    axs[1].set_title('b) Estructura Θ–SA y oxígeno disuelto')

    for mass in mass_order:
        g = x[x['masa_hidrografica_v2'].eq(mass)]
        if g.empty:
            continue
        alpha = 0.35 if mass == 'Sin clasificar' else 0.88
        size = 14 if mass == 'Sin clasificar' else 19
        axs[2].scatter(
            g['salinidad_absoluta_g_kg'], g['temperatura_conservativa_C'],
            color=mass_colors[mass], marker='o', s=size, alpha=alpha,
            edgecolors='none', linewidths=0.0, rasterized=True,
            zorder=4 if mass != 'Sin clasificar' else 2)
    axs[2].set_title('c) Estructura Θ–SA y masas de agua')

    swt_labels = {'mSUW': (10, 8), 'mAAIW': (10, 8), 'mNADW': (10, 8)}
    for ax in axs:
        for name in SWT.columns:
            xx = float(SWT.loc['SA_g_kg', name])
            yy = float(SWT.loc['CT_C', name])
            ax.scatter(xx, yy, marker='*', s=230, facecolor='white', edgecolor='k', linewidth=1.1, zorder=8)
            ax.annotate(name, (xx, yy), xytext=swt_labels.get(name, (8, 8)), textcoords='offset points',
                        fontsize=11.3, weight='bold',
                        bbox=dict(boxstyle='round,pad=0.18', fc='white', ec='none', alpha=0.86),
                        zorder=9)

    handles = [
        Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=mass_colors[m],
               markeredgecolor='none', markersize=8.5, label=m)
        for m in mass_order if not (m == 'Sin clasificar' and not x['masa_hidrografica_v2'].eq(m).any())
    ]
    handles.append(Line2D([0], [0], marker='*', linestyle='None', markerfacecolor='white',
                          markeredgecolor='k', color='k', markersize=11, label='SWT OMP'))
    axs[2].legend(handles=handles, loc='lower right', frameon=True, title='Clasificación hidrográfica',
                  fontsize=CFG['font_legend'], title_fontsize=CFG['font_legend'] + 0.5)

    fig.suptitle('Diagramas Θ–SA conjuntos: profundidad, oxígeno y clasificación hidrográfica',
                 fontsize=CFG['font_suptitle'], weight='bold')
    savefig(fig, out / 'Fig02_V3_4_TS_tripanel_profundidad_oxigeno_masas.png')

def fig_perfiles_medianos_v31(base, nucleos, o2col, out):
    q = filtro_hidrografico_pass(base)
    if q.empty:
        q = filtro_hidrografico(base)
    if q.empty:
        return
    q = _add_o2_mgl_column(q, o2col)
    dz = 10.0
    q = q.copy()
    q["zbin_m"] = np.round(pd.to_numeric(q["profundidad_m"], errors="coerce") / dz) * dz
    vars_ = [
        ("salinidad_absoluta_g_kg", "Salinidad absoluta, SA [g kg⁻¹]"),
        ("temperatura_conservativa_C", "Temperatura conservativa, Θ [°C]"),
        ("sigma0_kg_m3", "Anomalía de densidad potencial, σ₀ [kg m⁻³]"),
    ]
    if "O2_mg_L" in q.columns:
        vars_.append(("O2_mg_L", "Oxígeno disuelto, O₂ [mg L⁻¹]"))

    fig, axs = plt.subplots(1, len(vars_), figsize=(4.9 * len(vars_), 7.4), sharey=True, layout="constrained", squeeze=False)
    axs = axs.ravel()

    stats = {}
    for c, _ in vars_:
        g = q.groupby("zbin_m")[c]
        stats[c] = pd.DataFrame({"med": g.median(), "q25": g.quantile(0.25), "q75": g.quantile(0.75)}).dropna(how="all")

    suw_med = -float(pd.to_numeric(nucleos.get("SUW_core_depth_m"), errors="coerce").median()) if "SUW_core_depth_m" in nucleos else np.nan
    aa_src = nucleos[~nucleos.get("AAIW_core_truncado_fondo", pd.Series(False, index=nucleos.index)).fillna(True)] if len(nucleos) else nucleos
    aaiw_med = -float(pd.to_numeric(aa_src.get("AAIW_core_depth_m"), errors="coerce").median()) if "AAIW_core_depth_m" in aa_src else np.nan

    zall = []
    for ax, (c, label) in zip(axs, vars_):
        st = stats[c]
        if st.empty:
            ax.set_visible(False)
            continue
        zplot = -np.abs(st.index.to_numpy(float))
        zall.extend(zplot[np.isfinite(zplot)].tolist())
        ax.fill_betweenx(zplot, st["q25"], st["q75"], alpha=0.18, label="IQR")
        ax.plot(st["med"], zplot, lw=2.2, label="Mediana")
        if np.isfinite(suw_med):
            ax.axhline(suw_med, ls="--", lw=1.3, label="SUW mediana")
        if np.isfinite(aaiw_med):
            ax.axhline(aaiw_med, ls=":", lw=1.5, label="AAIW mediana")
        ax.set_xlabel(label)
        ax.set_title(label)
        ax.grid(alpha=0.15)
    axs[0].set_ylabel("Coordenada vertical, z [m]")
    if zall:
        axs[0].set_ylim(float(np.nanmin(zall)) - 20, 0.0)
    visible = [ax for ax in axs if ax.get_visible()]
    if visible:
        visible[-1].legend(loc="best", fontsize=CFG["font_legend"])
    fig.suptitle("Estructura vertical mediana de la campaña (banda = IQR)", fontsize=CFG["font_suptitle"], weight="bold")
    savefig(fig, out / "Fig03_V3_4_perfiles_medianos_IQR_zneg.png")


def fig_secciones_hidrograficas_v31(base, o2col, out):
    """Secciones continuas de toda la campaña con una paleta distinta por variable."""
    hydro = filtro_hidrografico_pass(base)
    if hydro.empty:
        hydro = filtro_hidrografico(base)
    hydro = _add_o2_mgl_column(hydro, o2col)
    oxy = hydro[hydro['O2_mg_L'].notna()].copy() if 'O2_mg_L' in hydro else hydro.iloc[0:0].copy()
    present = sorted(set(pd.to_numeric(hydro['numero_estacion'], errors='coerce').dropna().astype(int).tolist()))
    if len(present) < 2:
        return
    variables = [
        ('salinidad_absoluta_g_kg', 'Salinidad absoluta, SA [g kg⁻¹]', hydro, 'viridis'),
        ('temperatura_conservativa_C', 'Temperatura conservativa, Θ [°C]', hydro, 'turbo'),
        ('sigma0_kg_m3', 'Anomalía de densidad potencial, σ₀ [kg m⁻³]', hydro, 'cividis'),
    ]
    if not oxy.empty:
        variables.append(('O2_mg_L', 'Oxígeno disuelto, O₂ [mg L⁻¹]', oxy, 'plasma'))
    periods = [p for p in ('Diurno', 'Nocturno') if (hydro['periodo'] == p).any()]
    if not periods:
        periods = sorted(hydro['periodo'].dropna().astype(str).unique().tolist())

    ranges = {col: _robust_vrange(src[col]) for col, _, src, _ in variables}

    fig, axs = plt.subplots(len(periods), len(variables), figsize=(5.4 * len(variables), 5.6 * len(periods)),
                             layout='constrained', sharex=True, sharey=True, squeeze=False)
    zmins = []
    for i, periodo in enumerate(periods):
        for j, (col, label, src, cmap) in enumerate(variables):
            ax = axs[i, j]
            sec = _section_matrix_stations(src, col, periodo, present)
            if sec is None:
                ax.text(0.5, 0.5, 'Sin datos válidos', transform=ax.transAxes, ha='center', va='center',
                        fontsize=CFG['font_title'], bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.9))
                ax.set_title(label)
                continue
            xd, z, M, labs, _ = sec
            zplot = -np.abs(z)
            zmins.append(float(np.nanmin(zplot)))
            lo, hi = ranges[col]
            pc = ax.pcolormesh(xd, zplot, M, shading='nearest', cmap=cmap, vmin=lo, vmax=hi)
            cb = fig.colorbar(pc, ax=ax, pad=0.02)
            cb.set_label(label, fontsize=CFG['font_base'])
            cb.ax.tick_params(labelsize=CFG['font_base'] - 1)
            if np.all(~np.isfinite(M)):
                ax.text(0.5, 0.5, 'Sin datos válidos', transform=ax.transAxes, ha='center', va='center',
                        fontsize=CFG['font_title'], bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.9))
            ax.set_title(label)
            ax.set_xticks(xd)
            ax.set_xticklabels(labs, rotation=45)
            ax.grid(alpha=0.10)
            if i == len(periods) - 1:
                ax.set_xlabel('Distancia acumulada [km]')
            if j == 0:
                ax.set_ylabel(f'{periodo}\nCoordenada vertical, z [m]')
    if zmins:
        zmin = min(zmins)
        for ax in axs.ravel():
            ax.set_ylim(zmin, 0.0)
    fig.suptitle('Secciones hidrográficas continuas de la campaña', fontsize=CFG['font_suptitle'], weight='bold')
    savefig(fig, out / 'Fig06_V3_4_secciones_hidrograficas_continuas.png')

def fig_omp_fracciones_secciones_v31(omp_diag, out):
    """Secciones continuas de fracciones OMP por periodo, publicables."""
    if omp_diag.empty:
        return
    present = sorted(set(pd.to_numeric(omp_diag["numero_estacion"], errors="coerce").dropna().astype(int).tolist()))
    if len(present) < 2:
        return
    periods = [p for p in ("Diurno", "Nocturno") if (omp_diag["periodo"] == p).any()]
    for periodo in periods:
        fig, axs = plt.subplots(1, 3, figsize=(6.2 * 3, 6.2), layout="constrained", sharey=True)
        zmins = []
        for ax, mass in zip(axs, ["mSUW", "mAAIW", "mNADW"]):
            col = f"frac_{mass}"
            sec = _section_matrix_stations(omp_diag, col, periodo, present)
            if sec is None:
                ax.text(0.5, 0.5, "Sin datos OMP válidos", transform=ax.transAxes, ha="center", va="center",
                        fontsize=CFG["font_title"], bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
                ax.set_title(mass)
                continue
            xd, z, M, labs, _ = sec
            M = 100.0 * M
            zplot = -np.abs(z)
            zmins.append(float(np.nanmin(zplot)))
            pc = ax.pcolormesh(xd, zplot, M, shading="nearest", cmap="viridis", vmin=0.0, vmax=100.0)
            try:
                if np.isfinite(M).any():
                    ax.contour(xd, zplot, M, levels=[30.0], colors='k', linewidths=1.0, alpha=0.85)
            except Exception:
                pass
            cb = fig.colorbar(pc, ax=ax, pad=0.02)
            cb.set_label(f"{mass} [%]", fontsize=CFG["font_axes"])
            cb.ax.tick_params(labelsize=CFG["font_base"])
            ax.set_title(mass)
            ax.set_xticks(xd)
            ax.set_xticklabels(labs, rotation=45)
            ax.set_xlabel("Distancia acumulada [km]")
            ax.grid(alpha=0.12)
        axs[0].set_ylabel("Coordenada vertical, z [m]")
        if zmins:
            for ax in axs:
                ax.set_ylim(min(zmins), 0.0)
        fig.suptitle(f"Fracciones OMP mSUW–mAAIW–mNADW · {periodo}", fontsize=CFG["font_suptitle"], weight="bold")
        savefig(fig, out / f"Fig07_V3_4_OMP_fracciones_seccion_{slug(periodo)}.png")


def fig_omp_residuales_v31(omp_diag, out):
    if omp_diag.empty:
        return
    q = _add_o2_mgl_column(omp_diag.rename(columns={"resid_O2_norm_umol_kg": "resid_O2_norm_tmp"}), "resid_O2_norm_tmp", outcol="resid_O2_norm_mg_L")
    cols = [
        ("resid_CT_norm_C", "Residual Θ [°C]"),
        ("resid_SA_norm_g_kg", "Residual SA [g kg⁻¹]"),
        ("resid_O2_norm_mg_L", "Residual O₂ [mg L⁻¹]"),
    ]
    fig, axs = plt.subplots(1, 3, figsize=(14.8, 6.8), layout="constrained", sharey=True)
    dz = 25.0
    q = q.copy()
    q["zbin"] = np.round(q.profundidad_m / dz) * dz
    zall = []
    for ax, (c, label) in zip(axs, cols):
        if c not in q:
            ax.set_visible(False)
            continue
        g = q.groupby("zbin")[c]
        med = g.median(); q25 = g.quantile(.25); q75 = g.quantile(.75)
        z = -np.abs(med.index.to_numpy(float)); zall.extend(z.tolist())
        ax.fill_betweenx(z, q25, q75, alpha=.2, label="IQR")
        ax.plot(med, z, lw=1.8, label="Mediana")
        ax.axvline(0, ls="--", lw=1.1)
        ax.set(xlabel=label, title=label)
        ax.grid(alpha=.18); ax.legend(fontsize=CFG["font_legend"])
    axs[0].set_ylabel("Coordenada vertical, z [m]")
    if zall:
        axs[0].set_ylim(float(np.nanmin(zall)) - 25, 0)
    fig.suptitle("Diagnóstico de residuales del OMP completo (fracciones normalizadas)", fontsize=CFG["font_title"])
    savefig(fig, out / "Fig08_V3_4_OMP_residuales_zneg.png")


def write_methodology_v31(out: Path, path_input: Path, o2col: str | None, omp_status: str):
    text = f"""# Metodología — masas de agua del Caribe colombiano — V3.4

Entrada: `{path_input}`

## Novedades de la versión 3.1
- El O₂ se expresa en **mg L⁻¹** en todas las figuras donde aparece como trazador o propiedad hidrográfica. En el panel de profundidad del Θ–SA, la coordenada vertical se muestra como z negativa (superficie arriba, fondo abajo).
- La figura principal Θ–SA se mantiene como **tripanel**: profundidad, oxígeno y clasificación hidrográfica.
- Las **secciones hidrográficas** se presentan de forma **continua a lo largo de todas las estaciones** disponibles, con una barra de color propia y una **paleta diferente para cada variable**. El mapa utiliza una extensión automática mínima para conservar referencia costera incluso cuando las estaciones se encuentran alejadas de tierra.
- Se incorporan **secciones OMP** continuas por periodo para las fracciones mSUW, mAAIW y mNADW.

## Convención vertical
Todas las figuras verticales usan `z = -|profundidad|` [m]. La superficie se representa en 0 m y el fondo hacia valores más negativos.

## Variables hidrográficas
- SA: salinidad absoluta [g kg⁻¹]
- Θ: temperatura conservativa [°C]
- σ₀: anomalía de densidad potencial [kg m⁻³]
- O₂: oxígeno disuelto [mg L⁻¹] en figuras

## O₂
Columna base detectada para O₂: `{o2col or 'NO DISPONIBLE'}`
Conversión usada para figuras: `O₂ [mg L⁻¹] ≈ O₂ [µmol kg⁻¹] × 0.0328`.

## OMP
- SWT utilizados: mSUW, mAAIW y mNADW (Torres et al., 2023).
- Pesos del OMP: {json.dumps(OMP_WEIGHTS, ensure_ascii=False)}
- Estado diagnóstico del OMP: **{omp_status}**
"""
    (out / "METODOLOGIA_V3_4.md").write_text(text, encoding="utf-8")


def write_captions_v31(out: Path):
    text = """# Pies de figura — V3.4

**Figura 1.** Localización de las estaciones CTD y profundidad máxima muestreada.

**Figura 2.** Diagramas conjuntos de temperatura conservativa (Θ) y salinidad absoluta (SA) para la campaña: (a) observaciones coloreadas por profundidad, (b) observaciones coloreadas por oxígeno disuelto en mg L⁻¹ y (c) observaciones coloreadas según la clasificación hidrográfica (CSW, SUW, aguas centrales, AAIW y NADW). Las estrellas indican los SWT del OMP (mSUW, mAAIW y mNADW).

**Figura 3.** Perfiles verticales medianos de salinidad absoluta (SA), temperatura conservativa (Θ), anomalía de densidad potencial (σ₀) y oxígeno disuelto (O₂, mg L⁻¹). La banda muestra el rango intercuartílico (IQR).

**Figura 4.** Estructura de la capa superficial: profundidad de la capa de mezcla (MLD), profundidad de la capa isotérmica (ILD) y espesor de la capa barrera (BLT).

**Figura 5.** Profundidad de los núcleos hidrográficos: máximo subsuperficial de salinidad (SUW) y mínimo intermedio de salinidad (AAIW).

**Figura 6.** Secciones hidrográficas continuas de la campaña para SA, Θ, σ₀ y O₂ (mg L⁻¹), con una fila por periodo de muestreo y una barra de color propia para cada variable.

**Figura 7.** Secciones continuas de las fracciones OMP para mSUW, mAAIW y mNADW, expresadas en porcentaje, para cada periodo de muestreo.

**Figura 8.** Residuales del OMP completo para Θ, SA y O₂ (mg L⁻¹) en función de la coordenada vertical z.

**Figura 9.** Sensibilidad de las fracciones OMP al uso de O₂, expresada como fracción media con la profundidad.

**Figura 10.** Comparación de las fracciones OMP en profundidades estándar para mSUW, mAAIW y mNADW.
"""
    (out / "PIES_DE_FIGURA_V3_4.md").write_text(text, encoding="utf-8")


def main():
    path = Path(CFG["archivo_entrada"]) if CFG.get("archivo_entrada") else Path(input("Ruta del archivo Excel/CSV: ").strip())
    out = Path(CFG["carpeta_salida"]); out.mkdir(parents=True, exist_ok=True)
    figout = out / "figuras_publicacion"; figout.mkdir(exist_ok=True)
    base, decisiones, pcol, o2col, profile_col = cargar_dataset(path)
    base["z_m"] = -np.abs(pd.to_numeric(base["profundidad_m"], errors="coerce"))
    base = _add_o2_mgl_column(base, o2col)
    print(f"Perfiles: {base[profile_col].nunique()} | bins: {len(base):,} | O2: {o2col}")

    nuc = []
    for pid, g in base.groupby(profile_col, sort=False):
        nuc.append(analizar_perfil_hidrografico_v2(g, pcol, o2col, profile_id=pid))
    nucleos = pd.DataFrame(nuc)

    omp, excl = ejecutar_omp(base, nucleos, o2col, profile_col)
    omp_ts = ejecutar_omp_ts(base, nucleos, profile_col)
    omp_diag, omp_stats = augmentar_diagnostico_omp(omp, omp_ts)
    if not omp_diag.empty:
        omp_diag["z_m"] = -np.abs(pd.to_numeric(omp_diag["profundidad_m"], errors="coerce"))
        omp_diag = _add_o2_mgl_column(omp_diag, "O2_umol_kg")
    if not omp_ts.empty:
        omp_ts["z_m"] = -np.abs(pd.to_numeric(omp_ts["profundidad_m"], errors="coerce"))

    omp_status = "NO_EVALUADO"
    if not omp_stats.empty and "resid_O2_norm_umol_kg" in set(omp_stats.variable):
        rr = omp_stats[omp_stats.variable.eq("resid_O2_norm_umol_kg")].iloc[0]
        if abs(rr.mediana) >= CFG["omp_o2_bias_warn_umol_kg"] and rr.fraccion_mismo_signo >= CFG["omp_same_sign_warn_fraction"]:
            omp_status = "EXPLORATORIO_DESACOPLE_O2_SISTEMATICO"
        else:
            omp_status = "SIN_ALARMA_SISTEMATICA_O2"

    nucleos = enriquecer_nucleos_v2(base, nucleos, omp_diag, o2col, profile_col)
    omp_res = resumen_omp_por_perfil(omp_diag) if not omp_diag.empty else pd.DataFrame()
    omp_std = omp_en_profundidades_estandar(omp_diag) if not omp_diag.empty else pd.DataFrame()
    ts_std = _interp_omp_standard(omp_ts, ["frac_mSUW", "frac_mAAIW", "frac_mNADW"])
    hydro_summary = resumen_hidrografico_robusto_v2(nucleos)
    master = tabla_maestra_v2(nucleos, omp_res)

    sens = []
    if not omp_std.empty or not ts_std.empty:
        for z in CFG["profundidades_estandar_m"]:
            a = omp_std[omp_std.prof_objetivo_m.eq(z)] if not omp_std.empty else pd.DataFrame()
            b = ts_std[ts_std.prof_objetivo_m.eq(z)] if not ts_std.empty else pd.DataFrame()
            for mass in ["mSUW", "mAAIW", "mNADW"]:
                sens.append({"profundidad_m": z, "z_m": -abs(float(z)), "masa": mass,
                             "OMP_full_media": float(a[f"frac_{mass}"].mean()) if len(a) else np.nan,
                             "OMP_TS_media": float(b[f"frac_{mass}"].mean()) if len(b) else np.nan,
                             "diferencia_full_menos_TS": float(a[f"frac_{mass}"].mean() - b[f"frac_{mass}"].mean()) if len(a) and len(b) else np.nan})
    sens = pd.DataFrame(sens)

    clasif = preparar_dataset_masas_v2(base, nucleos, omp_diag, profile_col, o2col)
    clasif["z_m"] = -np.abs(pd.to_numeric(clasif["profundidad_m"], errors="coerce"))
    clasif = _add_o2_mgl_column(clasif, o2col)

    master.to_csv(out / "MASAS_AGUA_POR_PERFIL_V3_4.csv", index=False, encoding="utf-8-sig")
    hydro_summary.to_csv(out / "RESUMEN_HIDROGRAFICO_ROBUSTO_V3_4.csv", index=False, encoding="utf-8-sig")
    clasif.to_csv(out / "CLASIFICACION_HIDROGRAFICA_POR_BIN_V3_4.csv", index=False, encoding="utf-8-sig")
    if not omp_diag.empty: omp_diag.to_csv(out / "OMP_FULL_DIAGNOSTICO_V3_4.csv", index=False, encoding="utf-8-sig")
    if not omp_ts.empty: omp_ts.to_csv(out / "OMP_TS_SENSIBILIDAD_V3_4.csv", index=False, encoding="utf-8-sig")
    if not omp_stats.empty: omp_stats.to_csv(out / "OMP_RESIDUALES_RESUMEN_V3_4.csv", index=False, encoding="utf-8-sig")
    if not sens.empty: sens.to_csv(out / "OMP_SENSIBILIDAD_PROF_ESTANDAR_V3_4.csv", index=False, encoding="utf-8-sig")
    if not excl.empty: excl.to_csv(out / "OMP_PERFILES_EXCLUIDOS_V3_4.csv", index=False, encoding="utf-8-sig")

    with pd.ExcelWriter(out / "RESULTADOS_MASAS_AGUA_V3_4.xlsx", engine="openpyxl") as wr:
        master.to_excel(wr, sheet_name="Masas_por_perfil", index=False)
        hydro_summary.to_excel(wr, sheet_name="Resumen_robusto", index=False)
        clasif.to_excel(wr, sheet_name="Clasificacion_bins", index=False)
        if not omp_diag.empty: omp_diag.to_excel(wr, sheet_name="OMP_full", index=False)
        if not omp_ts.empty: omp_ts.to_excel(wr, sheet_name="OMP_TS_sens", index=False)
        if not omp_stats.empty: omp_stats.to_excel(wr, sheet_name="OMP_residuales", index=False)
        if not sens.empty: sens.to_excel(wr, sheet_name="OMP_prof_estandar", index=False)
        if not excl.empty: excl.to_excel(wr, sheet_name="OMP_excluidos", index=False)
        SWT.to_excel(wr, sheet_name="SWT_Torres2023")
        pd.DataFrame([CFG]).astype(str).to_excel(wr, sheet_name="Configuracion", index=False)

    fig_mapa_estaciones_vf(base, figout)
    fig_ts_tripanel_v31(base, nucleos, omp_diag, o2col, profile_col, figout)
    fig_perfiles_medianos_v31(base, nucleos, o2col, figout)
    fig_capas_superficiales_v3(nucleos, figout)
    fig_nucleos_v3(nucleos, figout)
    fig_secciones_hidrograficas_v31(base, o2col, figout)
    fig_omp_fracciones_secciones_v31(omp_diag, figout)
    fig_omp_residuales_v31(omp_diag, figout)
    fig_omp_sensibilidad_v3(omp_diag, omp_ts, figout)
    fig_omp_estandar_v3(omp_diag, omp_ts, figout)

    write_methodology_v31(out, path, o2col, omp_status)
    write_captions_v31(out)

    meta = {"version": VERSION, "input": str(path), "output": str(out), "source_reference": SOURCE_REFERENCE,
            "SWT": SWT.to_dict(), "OMP_weights": OMP_WEIGHTS, "OMP_status": omp_status,
            "sintesis_hidrografica": ["CSW", "SUW", "Aguas centrales", "AAIW", "NADW"],
            "endmembers_OMP": ["mSUW", "mAAIW", "mNADW"],
            "convencion_vertical": "z=-|profundidad| m; superficie=0; fondo=hacia valores más negativos",
            "figura_TS": "tripanel: profundidad / O2 mg L-1 / masa identificada por color",
            "nota": "Secciones hidrográficas y OMP continuas por estación. O2 mostrado en mg/L en todas las figuras.",
            "transectos_configurados": CFG["transectos"]}
    (out / "metadata_masas_agua_v3_4.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")

    print("\nRESUMEN HIDROGRÁFICO ROBUSTO")
    print(hydro_summary.to_string(index=False))
    print("\nDIAGNÓSTICO OMP:", omp_status)
    if not omp_stats.empty:
        print(omp_stats.to_string(index=False))
    print("\nResultados:", out)
    print("Figuras publicación:", figout)
    return master, omp_diag, omp_ts, clasif


if __name__ == "__main__":
    MASAS_PERFIL, OMP_FULL, OMP_TS, CLASIFICACION_BINS = main()
