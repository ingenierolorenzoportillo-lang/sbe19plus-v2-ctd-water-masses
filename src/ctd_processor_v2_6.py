# -*- coding: utf-8 -*-
"""CTD Crucero Auto para Spyder — versión 2.6 DISCOVERY-VF-ESTATICOS.

Basado en el asistente v1, pero el flujo operativo principal es automático:
infiere límites de lance, soak/atemperamiento, Loop Edit adaptativo y ajustes
dinámicos de oxígeno; solo envía perfiles ambiguos a una cola de revisión.

Entrada: Excel/CSV de estaciones + HEX originales + XMLCON de calibración.
Decodificador implementado: SBE 19plus / 19plus V2, modo perfil,
presión strain gauge, HEX raw (OutputFormat=0). Ver GUIA_CTD.md.
No requiere SBE Data Processing. No equivale a una certificación de datos.

Dependencias (instalar en el MISMO entorno de Spyder):
    python -m pip install numpy pandas scipy matplotlib openpyxl gsw seabirdscientific

Ejecutar el archivo completo con F5. El modo principal NO pregunta decisiones
perfil a perfil. Las figuras finales se guardan en PNG 300 dpi y los perfiles
ambiguos se marcan como REVISAR.
Los originales no se modifican. En modo incremental la sesión persiste y solo se reprocesan HEX nuevos o modificados.
"""
from __future__ import annotations

import sys
import re
import json
import math
import hashlib
import shutil
import sqlite3
import traceback
import textwrap
import unicodedata
import importlib.metadata
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone, timedelta
from pathlib import Path
import xml.etree.ElementTree as ET

try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, ScalarFormatter
    from scipy import signal
    import gsw
    import seabirdscientific.conversion as sbconv
    import seabirdscientific.cal_coefficients as sbcoef
    import seabirdscientific.processing as sbproc
except ImportError as exc:
    raise SystemExit(
        f"Falta una dependencia: {exc}\nIntérprete de Spyder: {sys.executable}\n"
        'En la consola IPython ejecuta:\n'
        '%pip install numpy pandas scipy matplotlib openpyxl gsw==3.6.20 seabirdscientific==2.7.10\n'
        'Después reinicia el kernel y ejecuta este archivo con F5.'
    ) from exc

VERSION = "2.6-DISCOVERY-VF-ESTATICOS"
BAD = -9.990e-29
FUENTES = {
    "manual_institucional": "INVEMAR, Programa GEO (noviembre de 2019). Manual de procesamiento del CTDO SBE 19 Plus. Documento aportado por el usuario.",
    "toolkit": "https://github.com/Sea-BirdScientific/seabirdscientific",
    "procesamiento": "https://github.com/Sea-BirdScientific/seabirdscientific/blob/main/documentation/processing.ipynb",
    "sbe19plus_v2": "https://www.seabird.com/products/sbe-19plus-v2-seacat-ctd",
    "formato_hex": "https://www.comm-tec.com/Docs/Manuali/SBE/19plusV2_011.pdf",
    "oxigeno": "https://www.seabird.com/hubfs/Resources/Application%20Notes/AN64-SBE43BackgroundInfoDeploymentMaintenance.pdf",
    "histeresis": "https://www.seabird.com/hubfs/Resources/Application%20Notes/AN64-3-SBEHysterisisCorrection.pdf",
    "teos10": "https://teos-10.github.io/GSW-Python/",
    "go_ship_manual": "https://www.pmel.noaa.gov/pubs/PDF/hood3609/hood3609.pdf",
    "qartod_ts": "https://cdn.ioos.noaa.gov/media/2020/03/QARTOD_TS_Manual_Update2_200324_final.pdf",
    "qartod_oxygen": "https://cdn.ioos.noaa.gov/media/2017/12/qartod_dissolved_oxygen_manual.pdf",
    "qartod_flags": "https://cdn.ioos.noaa.gov/media/2017/12/QARTOD-Data-Flags-Manual_Final_version1.1.pdf",
    "cf_conventions": "https://cfconventions.org/Data/cf-conventions/cf-conventions-1.13/cf-conventions.html",
}

# QARTOD/IOOS primary flags. These are kept separate from the internal bitmask QC.
QARTOD = {
    1: "pass",
    2: "not_evaluated",
    3: "suspect_or_high_interest",
    4: "fail",
    9: "missing",
}

# Flags internos de esta implementación, no una certificación QARTOD/WOCE.
QC = {1: "no finito", 2: "fuera del rango configurado", 4: "pico Hampel",
      8: "borde sin soporte de alineamiento/corrección", 16: "fuera del tramo elegido",
      32: "Loop Edit: velocidad o inversión", 64: "presión menor al mínimo elegido",
      128: "transitorio tras reiniciar corrección de masa térmica; revisar"}

@dataclass
class Variable:
    nombre: str
    unidad: str
    cnv: str
    origen: str = "derivada"
    graficar: bool = True

BASE_VARIABLES = {
    "presion_dbar": Variable("Presión", "dbar", "prdM", "sensor"),
    "temperatura_C": Variable("Temperatura in situ", "°C, ITS-90", "tv290C", "sensor"),
    "conductividad_S_m": Variable("Conductividad", "S m⁻¹", "c0S/m", "sensor"),
    "profundidad_m": Variable("Profundidad", "m", "depSM", graficar=False),
    "salinidad_practica": Variable("Salinidad práctica", "1, PSS-78", "sal00"),
    "salinidad_absoluta_g_kg": Variable("Salinidad absoluta", "g kg⁻¹", "saliA0"),
    "temperatura_conservativa_C": Variable("Temperatura conservativa", "°C", "t090C_CT"),
    "temperatura_potencial_C": Variable("Temperatura potencial, 0 dbar", "°C", "potemp090C"),
    "densidad_kg_m3": Variable("Densidad in situ", "kg m⁻³", "density00"),
    "sigma0_kg_m3": Variable("Anomalía de densidad potencial", "kg m⁻³", "sigma00"),
    "velocidad_sonido_m_s": Variable("Velocidad del sonido", "m s⁻¹", "soundVel"),
}

@dataclass
class Sensor:
    indice: int
    etiqueta: str
    sensor_id: str
    serie: str
    fecha_calibracion: str
    coef: dict
    tipo: str = "desconocido"
    senal: str = ""
    salida: str = ""
    motivo: str = ""

@dataclass
class Configuracion:
    ruta: Path
    nombre: str
    ajustes: dict
    sensores: list[Sensor]
    nvolt: int
    navg: int
    serial: str = ""
    avisos: list = field(default_factory=list)

@dataclass
class Lance:
    ruta: Path
    meta: dict
    config: Configuracion
    encabezado: str
    raw: pd.DataFrame
    datos: pd.DataFrame
    variables: dict
    dt: float
    hora_hex: str = ""
    fuente_hora: str = ""
    avisos: list = field(default_factory=list)
    plan: dict = field(default_factory=dict)

def norm(x):
    x = unicodedata.normalize("NFKD", str(x))
    return re.sub(r"[^a-z0-9]+", "_", "".join(c for c in x if not unicodedata.combining(c)).lower()).strip("_")

def slug(x):
    return norm(x)[:110] or "lance"

def si(pregunta, defecto=None):
    suf = " [y/n]" if defecto is None else (" [Y/n]" if defecto else " [y/N]")
    while True:
        a = input(pregunta + suf + ": ").strip().lower()
        if not a and defecto is not None:
            return defecto
        if a in ("y", "s", "si", "sí", "yes"):
            return True
        if a in ("n", "no"):
            return False
        print("Escribe y para sí o n para no.")

def numero(pregunta, defecto=None, minimo=None, maximo=None, entero=False):
    while True:
        a = input(pregunta + (f" [{defecto}]" if defecto is not None else "") + ": ").strip()
        try:
            v = float(a.replace(",", ".")) if a else float(defecto)
            if not np.isfinite(v) or (entero and not v.is_integer()):
                raise ValueError
            if minimo is not None and v < minimo or maximo is not None and v > maximo:
                raise ValueError
            return int(v) if entero else v
        except (ValueError, TypeError):
            print(f"Valor inválido. Intervalo permitido: {minimo} a {maximo}.")

def elegir(pregunta, opciones, defecto=1):
    print("\n" + pregunta)
    for i, o in enumerate(opciones, 1):
        print(f"  {i}. {o}")
    return numero("Opción", defecto, 1, len(opciones), True)

def ruta_input(texto, carpeta=False, extensiones=()):
    while True:
        s = input(texto + "\nPega la ruta, o pulsa Enter para abrir el selector: ").strip().strip('"').strip("'")
        if not s:
            try:
                import tkinter as tk
                from tkinter import filedialog
                root = tk.Tk(); root.withdraw(); root.attributes("-topmost", True)
                s = filedialog.askdirectory(title=texto) if carpeta else filedialog.askopenfilename(title=texto)
                root.destroy()
            except Exception as exc:
                print(f"No se pudo abrir el selector ({exc}). Pega la ruta en la consola.")
                continue
        if not s:
            if si("¿Cancelar la ejecución?", False):
                raise KeyboardInterrupt
            continue
        p = Path(s).expanduser()
        if carpeta and not p.is_dir() or not carpeta and not p.is_file():
            print("La ruta no existe o no es del tipo solicitado."); continue
        if extensiones and p.suffix.lower() not in extensiones:
            print("Extensiones admitidas: " + ", ".join(extensiones)); continue
        return p.resolve()

def tabla_visual(df, titulo, guardar=None, mostrar=True):
    print("\n" + titulo + "\n" + df.to_string(index=False))
    for ini in range(0, len(df), 18):
        sub = df.iloc[ini:ini+18].astype(str)
        vals = [[textwrap.fill(c, 43) for c in row] for row in sub.to_numpy()]
        lineas = [max([v.count("\n") + 1 for v in row] + [1]) for row in vals]
        fig, ax = plt.subplots(figsize=(min(19, max(10, len(df.columns)*2.8)), max(3,1.5 + sum(lineas)*.34)))
        ax.axis("off")
        tab = ax.table(cellText=vals, colLabels=list(df.columns), cellLoc="left", bbox=[0,0,1,.98])
        tab.auto_set_font_size(False); tab.set_fontsize(10)
        total = max(sum(lineas) + 2, 3)
        for (r, c), cell in tab.get_celld().items():
            cell.set_edgecolor("#D3DBE4")
            cell.set_height((1.7 if r == 0 else lineas[r-1]*1.25)/total)
            cell.set_facecolor("#183D59" if r == 0 else ("#EDF3F7" if r % 2 else "white"))
            if r == 0: cell.set_text_props(color="white", weight="bold")
        fig.suptitle(textwrap.fill(titulo, 100), x=.025,y=.985,ha="left",fontsize=13,weight="bold")
        fig.subplots_adjust(left=.025,right=.975,bottom=.025,top=.83)
        if guardar:
            fig.savefig(str(guardar) + f"_{ini//18+1}.png", dpi=180, bbox_inches="tight")
        if mostrar: mostrar_fig(fig)
        else: plt.close(fig)

def mostrar_fig(fig):
    plt.show(block=False)
    # draw() permite visualizar también cuando Spyder usa el panel Gráficos.
    fig.canvas.draw_idle()
    try:
        from IPython import get_ipython
        if get_ipython() is not None and "inline" in plt.get_backend().lower():
            from IPython.display import display
            display(fig)
    except ImportError:
        pass

def leer_tabla(p):
    if p.suffix.lower() in (".xlsx", ".xlsm"):
        df = pd.read_excel(p, sheet_name=0)
    else:
        for enc in ("utf-8-sig", "cp1252"):
            try:
                df = pd.read_csv(p, sep=None, engine="python", encoding=enc); break
            except UnicodeDecodeError:
                continue
    df = df.dropna(how="all").copy()
    nombres = [norm(c) for c in df.columns]
    if len(nombres) != len(set(nombres)):
        raise ValueError("Hay encabezados de columna duplicados en la tabla.")
    df.columns = nombres
    alias = {"nombre": "estacion", "lat": "latitud", "lon": "longitud", "long": "longitud",
             "fecha": "fecha_bitacora", "hora": "hora_inmersion", "archivo": "archivo_hex"}
    df.rename(columns={k:v for k,v in alias.items() if k in df and v not in df}, inplace=True)
    faltan = set(("estacion", "latitud", "longitud")) - set(df)
    if faltan: raise ValueError("Faltan columnas en la primera hoja: " + ", ".join(faltan))
    for c, lim in (("latitud", 90), ("longitud", 180)):
        df[c] = pd.to_numeric(df[c].astype(str).str.replace(",", ".", regex=False), errors="coerce")
        if df[c].isna().any() or (df[c].abs() > lim).any():
            raise ValueError(f"{c}: hay valores vacíos, no numéricos o fuera de ±{lim} grados.")
    if df.estacion.isna().any() or df.estacion.astype(str).str.strip().eq("").any():
        raise ValueError("Cada fila necesita un nombre de estación.")
    df["estacion"] = df.estacion.astype(str).str.strip().str.replace(r"\s+", " ", regex=True)
    if "fecha_bitacora" in df:
        def fecha(v):
            if pd.isna(v) or str(v).strip() == "": return ""
            if isinstance(v, (datetime, pd.Timestamp)): return v.strftime("%Y-%m-%d")
            s = str(v).strip()
            f = pd.to_datetime(s, dayfirst=not bool(re.match(r"^\d{4}-", s)), errors="raise")
            return f.strftime("%Y-%m-%d")
        df["fecha_bitacora"] = df.fecha_bitacora.map(fecha)
    if "hora_inmersion" in df:
        def hora(v):
            if pd.isna(v) or str(v).strip() == "": return ""
            if isinstance(v, (float, int)) and 0 <= v < 1:
                n = round(v*86400); return f"{n//3600:02d}:{n%3600//60:02d}:{n%60:02d}"
            if hasattr(v, "strftime"): return v.strftime("%H:%M:%S")
            s=str(v).strip()
            for formato in ("%H:%M:%S","%H:%M"):
                try:return datetime.strptime(s,formato).strftime("%H:%M:%S")
                except ValueError:pass
            raise ValueError(f"Hora de inmersión inválida: {s}. Usa HH:MM:SS o una celda de hora de Excel.")
        df["hora_inmersion"] = df.hora_inmersion.map(hora)
    return df.reset_index(drop=True).fillna("")

def vincular(df, carpeta):
    archivos = sorted(p for p in carpeta.rglob("*") if p.is_file() and p.suffix.lower() == ".hex")
    if not archivos: raise ValueError("No se encontraron archivos HEX en la carpeta y subcarpetas.")
    usados = set(); rutas = []
    for i, row in df.iterrows():
        nombre = str(row.get("archivo_hex", "")).strip().replace("\\", "/")
        candidatos = [p for p in archivos if p.relative_to(carpeta).as_posix().casefold() == nombre.casefold()]
        if not candidatos:
            candidatos = [p for p in archivos if p.name.casefold() == nombre.casefold()]
        if nombre and len(candidatos) != 1:
            raise ValueError(f"Fila {i+2}: '{nombre}' no existe o coincide con varios archivos. Usa su ruta relativa exacta.")
        if not nombre:
            disponibles = [p for p in archivos if p not in usados]
            print(f"\nAsociación manual para {row.estacion}. No se asociará por proximidad horaria.")
            for j, p in enumerate(disponibles, 1): print(f" {j}: {p.relative_to(carpeta)}")
            if not disponibles: raise ValueError("Hay más filas de estaciones que HEX disponibles.")
            p = disponibles[numero("Número del HEX correspondiente", minimo=1, maximo=len(disponibles), entero=True)-1]
        else: p = candidatos[0]
        if p in usados: raise ValueError(f"El archivo {p.name} está asociado a más de una fila.")
        usados.add(p); rutas.append(p)
    df["archivo_hex"] = [p.relative_to(carpeta).as_posix() for p in rutas]
    extra = [p.relative_to(carpeta).as_posix() for p in archivos if p not in usados]
    if extra:
        print("HEX sin fila de estación:\n" + "\n".join(extra))
        if not si("¿Continuar procesando solo los HEX asociados a la tabla?", False): raise KeyboardInterrupt
    tabla_visual(df[["archivo_hex", "estacion", "latitud", "longitud"]], "Relación de estaciones y archivos")
    if not si("¿Es correcta esta relación y las coordenadas están en grados decimales WGS84?", None):
        raise ValueError("Corrige la tabla de estaciones y vuelve a ejecutar.")
    plt.close("all")
    return df, rutas

def coeficientes(nodo):
    """Lee solo la rama de ecuación activa; no mezcla calibraciones antiguas y nuevas."""
    directos = {norm(c.tag): (c.text or "").strip() for c in nodo if len(c) == 0}
    ramas = [c for c in nodo if c.attrib.get("equation") is not None]
    if ramas:
        activo = directos.get("use2007equation", directos.get("useg_j"))
        if activo is None:
            if len(ramas) != 1: raise ValueError(f"Ecuación activa ambigua en {nodo.tag}.")
            rama = ramas[0]
        else:
            rr = [c for c in ramas if c.attrib["equation"] == str(int(float(activo)))]
            if len(rr) != 1: raise ValueError(f"No se encuentra la ecuación activa de {nodo.tag}.")
            rama = rr[0]
        directos.update({norm(c.tag):(c.text or "").strip() for c in rama.iter() if len(c)==0})
    else:
        for c in nodo.iter():
            if len(c)==0: directos.setdefault(norm(c.tag), (c.text or "").strip())
    return directos

def valor(s, *claves, defecto=None):
    for c in claves:
        v = s.coef.get(norm(c))
        if v not in (None, ""):
            n = float(str(v).replace("D", "E").replace("d", "e"))
            if not math.isfinite(n): raise ValueError(f"Coeficiente no finito: {s.etiqueta}.{c}")
            return n
    if defecto is not None: return defecto
    raise ValueError(f"Falta coeficiente {'/'.join(claves)} en {s.etiqueta} (serie {s.serie}).")

def clasificar(s):
    t = norm(s.etiqueta)
    if s.indice == 0 and t == "temperaturesensor": return "temperatura"
    if s.indice == 1 and t == "conductivitysensor": return "conductividad"
    if s.indice == 2 and t == "pressuresensor": return "presion"
    if t == "oxygensensor" and valor(s, "Use2007Equation", defecto=0) == 1: return "sbe43"
    if "sbe63" in t or ("oxygen" in t and "ta0" in s.coef and "a0" in s.coef): return "sbe63"
    if t in ("ph_sensor", "phsensor") and s.sensor_id == "43": return "ph18"
    if "eco" in t:  # conversión ECO solo si se reconocen coeficientes
        if "scalefactor" in s.coef and any(k in s.coef for k in ("vblank", "darkvoltage", "darkoutput")):
            return "eco"
    if t == "turbiditymeter" and "scalefactor" in s.coef and "darkvoltage" in s.coef: return "eco"
    if t == "notinuse": return "sin_uso"
    return "desconocido"

def leer_config(ruta):
    root = ET.parse(ruta).getroot()
    for e in root.iter(): e.tag = e.tag.split("}")[-1]
    ins = root.find(".//Instrument")
    if ins is None: raise ValueError("XMLCON sin nodo Instrument.")
    opts = {norm(c.tag):(c.text or "").strip() for c in ins if len(c)==0}
    nombre = opts.get("name", "desconocido")
    sensores = []
    for n in ins.findall("./SensorArray/Sensor"):
        if not len(n): raise ValueError("Sensor sin configuración en XMLCON.")
        e = n[0]; c = coeficientes(e)
        s = Sensor(int(n.attrib["index"]), e.tag, n.attrib.get("SensorID", ""),
                   c.get("serialnumber", ""), c.get("calibrationdate", ""), c)
        s.tipo = clasificar(s); sensores.append(s)
    sensores.sort(key=lambda x:x.indice)
    if len(set(s.indice for s in sensores)) != len(sensores): raise ValueError("Índices de sensor duplicados.")
    if not re.search(r"19\s*plus", nombre, re.I):
        invent = ", ".join(s.etiqueta for s in sensores)
        raise ValueError(f"Equipo identificado: {nombre}. Sensores: {invent}.\n"
                         "Este decodificador HEX admite SBE 19plus/19plus V2; otro equipo necesita su adaptador específico.")
    if int(opts.get("mode", "0")) != 0: raise ValueError("XMLCON en modo fondeo: este flujo es para perfiles.")
    if int(opts.get("pressuresensortype", "1")) != 1:
        raise ValueError("Presión Digiquartz u otro tipo detectada: requiere validar su formato con archivos reales; esta versión admite strain gauge.")
    nvolt = int(opts.get("externalvoltagechannels", "-1"))
    navg = int(opts.get("scanstoaverage", "0"))
    if not 0 <= nvolt <= 6 or navg < 1: raise ValueError("ExternalVoltageChannels o ScansToAverage inválidos/ausentes.")
    if [s.tipo for s in sensores[:3]] != ["temperatura", "conductividad", "presion"]:
        raise ValueError("La disposición T-C-P del XMLCON no corresponde al SBE 19plus.")
    if [s.indice for s in sensores] != list(range(len(sensores))):
        raise ValueError("SensorArray debe tener índices continuos desde 0; no se adivinarán canales.")
    for flag in ("surfaceparvoltageadded", "nmeadepthdataadded", "scantimeadded"):
        if int(opts.get(flag, "0")):
            raise ValueError(f"{flag} activo: este formato anexado aún no tiene adaptador validado.")
    for s in sensores[3:3+nvolt]: s.senal = f"volt_{s.indice-3}"
    seriales = sensores[3+nvolt:]
    serial = ""
    if seriales:
        if len(seriales)==1 and seriales[0].tipo == "sbe63":
            serial="sbe63"; seriales[0].senal="sbe63_fase"
        else:
            raise ValueError("Canales RS-232 no implementados: " + ", ".join(s.etiqueta for s in seriales))
    if len(sensores) < 3+nvolt: raise ValueError("Faltan sensores del SensorArray para los canales activos.")
    cfg = Configuracion(ruta, nombre, opts, sensores, nvolt, navg, serial)
    for s in sensores:
        if not s.fecha_calibracion: cfg.avisos.append(f"Sin fecha de calibración: {s.etiqueta}, índice {s.indice}.")
        if s.tipo=="desconocido": s.motivo="Conversión física no implementada; se conserva señal en voltios."
        # Fallar temprano ante calibraciones incompletas, sin reemplazarlas por cero.
        requeridos={"sbe43":("Soc","offset","A","B","C","E"),
                    "sbe63":("A0","A1","A2","B0","B1","C0","C1","C2","E","TA0","TA1","TA2","TA3"),
                    "ph18":("Slope","Offset"),"eco":("ScaleFactor",)}
        for k in requeridos.get(s.tipo,()): valor(s,k)
        if s.tipo=="sbe43" and valor(s,"Soc")<=0: raise ValueError("Soc de oxígeno debe ser positivo.")
        if s.tipo=="ph18" and valor(s,"Slope")==0: raise ValueError("Pendiente de pH igual a cero.")
    return cfg

def header_value(header, patrones):
    for pat in patrones:
        m = re.search(pat, header, re.I | re.M)
        if m: return m.group(1).strip()
    return ""

def fecha_header(header):
    # System Upload Time y reloj de consulta DS NO son el inicio de adquisición.
    eventos = re.findall(r"\bcast\s+\d+\s+(.+?)\s+samples\s+\d+\s+to\s+\d+", header, re.I)
    if len(eventos)>1: return "", "varios inicios de cast en el encabezado; requiere revisión"
    s = header_value(header, [r"^\s*[*#]+\s*start_time\s*=\s*([^\[\n]+)"])
    fuente = "start_time explícito"
    if not s and eventos: s=eventos[0]; fuente="inicio de cast, memoria del instrumento"
    if not s: return "", "no hay inicio inequívoco; hora de descarga ignorada"
    try: return pd.to_datetime(s).isoformat(), fuente
    except (ValueError, TypeError): return "", f"inicio sin interpretar: {s}"

def leer_hex(ruta, cfg):
    texto = ruta.read_text(encoding="cp1252")
    if "*END*" not in texto: raise ValueError("HEX sin delimitador *END*: formato no compatible.")
    header, body = texto.split("*END*",1)
    lineas = [(i+1,l.strip()) for i,l in enumerate(body.splitlines()) if l.strip()]
    if len(lineas)<10: raise ValueError("El archivo tiene menos de 10 scans.")
    # Validación completa: una fila inválida se informa, nunca se salta comprimiendo el tiempo.
    nmea_pos = int(cfg.ajustes.get("nmeapositiondataadded","0"))
    nmea_time = int(cfg.ajustes.get("nmeatimeadded","0"))
    ancho = 22 + cfg.nvolt*4 + (12 if cfg.serial=="sbe63" else 0) + 14*nmea_pos + 8*nmea_time
    for i,l in lineas:
        if len(l)!=ancho or not re.fullmatch(r"[0-9a-fA-F]+",l):
            raise ValueError(f"Fila de datos {i}: {len(l)} caracteres; se esperaban {ancho}.\n"
                             "HEX truncado, formato distinto o XMLCON incompatible. No se omitió la fila.")
    serie = header_value(header,[r"Temperature\s+SN\s*=\s*(\d+)", r"SERIAL\s+NO\.?\s*(\d+)"])
    if serie and serie.lstrip("0") != cfg.sensores[0].serie.lstrip("0"):
        raise ValueError(f"Serie T en HEX ({serie}) distinta a XMLCON ({cfg.sensores[0].serie}).")
    serie_c=header_value(header,[r"Conductivity\s+SN\s*=\s*(\d+)"])
    if serie_c and serie_c.lstrip("0")!=cfg.sensores[1].serie.lstrip("0"):
        raise ValueError("La serie de conductividad no coincide entre HEX y XMLCON.")
    modelo = header_value(header,[r"Sea-Bird\s+(SBE\S+)\s+Data", r"DeviceType=['\"]([^'\"]+)"])
    if modelo and "19plus" not in modelo.lower(): raise ValueError("Modelo en HEX incompatible: "+modelo)
    if re.search(r"mode\s*=\s*moored|<MooredMode>", header, re.I): raise ValueError("HEX de fondeo, no de perfil.")
    navg = header_value(header,[r"<ScansToAverage>\s*(\d+)",r"number of scans to average\s*=\s*(\d+)"])
    if navg and int(navg)!=cfg.navg: raise ValueError("ScansToAverage no coincide entre HEX y XMLCON.")
    activos=[]; vistos=0
    for i in range(6):
        a=header_value(header,[rf"<ExtVolt{i}>\s*(yes|no)",rf"Ext\s+Volt\s+{i}\s*=\s*(yes|no)"])
        if a: vistos+=1
        if a.lower()=="yes": activos.append(i)
    if vistos==6 and len(activos)!=cfg.nvolt: raise ValueError("Número de canales de voltaje distinto en HEX y XMLCON.")
    for etiqueta in ("SBE38","WETLABS","OPTODE","GTD","DualGTD"):
        v=header_value(header,[rf"<{etiqueta}>\s*(yes|no)",rf"\b{etiqueta}\s*=\s*(yes|no)"])
        if v.lower()=="yes": raise ValueError(f"Sensor serial {etiqueta} activo: formato no implementado.")
    v=header_value(header,[r"<SBE63>\s*(yes|no)",r"SBE63\s*=\s*(yes|no)"])
    if v and (v.lower()=="yes") != (cfg.serial=="sbe63"): raise ValueError("SBE63 no coincide entre HEX y XMLCON.")
    if re.search(r"pressure sensor\s*=\s*(quartz|digiquartz)",header,re.I): raise ValueError("HEX con presión Quartz.")
    filas=[s for _,s in lineas]
    def tomar(ini, n, factor=1, offset=0):
        return np.array([int(s[ini:ini+n],16) for s in filas],dtype=float)*factor+offset
    raw=pd.DataFrame({"scan":np.arange(1,len(filas)+1),"temp_counts":tomar(0,6),
        "cond_Hz":tomar(6,6,1/256),"pres_counts":tomar(12,6),"pres_temp_V":tomar(18,4,1/13107)})
    j=22
    for i in range(cfg.nvolt): raw[f"volt_{i}"]=tomar(j,4,1/13107); j+=4
    if cfg.serial=="sbe63":
        raw["sbe63_fase"]=tomar(j,6,1e-5,-10);raw["sbe63_temp_V"]=tomar(j+6,6,1e-6,-1);j+=12
    # Datos NMEA anexados son preservados en HEX; coordenadas de procesamiento = tabla confirmada.
    if nmea_pos: j+=14
    if nmea_time: j+=8
    if j!=ancho: raise AssertionError("Longitud de trama incoherente.")
    return header,raw

def convertir_base(raw,cfg):
    t,c,p=cfg.sensores[:3]
    tc=sbcoef.TemperatureCoefficients(*[valor(t,f"A{i}",f"TA{i}") for i in range(4)])
    T=sbconv.convert_temperature(raw.temp_counts.to_numpy(),tc,use_mv_r=True)
    T=T*valor(t,"Slope",defecto=1)+valor(t,"Offset","TOffset",defecto=0)
    pc=sbcoef.PressureCoefficients(**{k:valor(p,k) for k in sbcoef.PressureCoefficients.__annotations__})
    P=sbconv.convert_pressure(raw.pres_counts.to_numpy(),raw.pres_temp_V.to_numpy(),pc,units="dbar")
    P=P*valor(p,"Slope",defecto=1)+valor(p,"Offset","POffset",defecto=0)
    C=convertir_conductividad(raw,T,P,c)
    df=pd.DataFrame({"scan":raw.scan,"tiempo_s":(raw.scan-1)*cfg.navg/4,
                     "presion_dbar":P,"temperatura_C":T,"conductividad_S_m":C})
    vars={k:v for k,v in BASE_VARIABLES.items() if k in df}
    for s in cfg.sensores[3:]:
        if s.senal in raw:
            df[s.senal]=raw[s.senal]
            vars[s.senal]=Variable(f"{s.etiqueta}, señal", "µs" if s.tipo=="sbe63" else "V",f"raw{s.indice}","señal",s.tipo=="desconocido")
        if s.tipo=="sbe63":
            df["sbe63_temp_V"]=raw.sbe63_temp_V
            vars["sbe63_temp_V"]=Variable("SBE63 termistor, señal","V","sbe63TV","señal",False)
    return df,vars

def convertir_conductividad(raw,T,P,sensor):
    if valor(sensor,"ConductivityType",defecto=0)!=0:
        raise ValueError("Celda de conductividad de rango ampliado: falta adaptador validado.")
    cc=sbcoef.ConductivityCoefficients(**{k:valor(sensor,k,defecto=0 if k=="wbotc" else None)
                                       for k in sbcoef.ConductivityCoefficients.__annotations__})
    C=sbconv.convert_conductivity(raw.cond_Hz.to_numpy(),np.asarray(T),np.asarray(P),cc)
    return C*valor(sensor,"Slope","CSlope",defecto=1)+valor(sensor,"Offset",defecto=0)

def bloques(mask):
    e=np.diff(np.r_[False,np.asarray(mask,bool),False].astype(int))
    return list(zip(np.flatnonzero(e==1),np.flatnonzero(e==-1)))

def filtrar(x,tau,dt):
    x=np.asarray(x,float); out=x.copy()
    if tau<=0: return out
    for a,b in bloques(np.isfinite(x)):
        if b-a>=4: out[a:b]=sbproc.low_pass_filter(x[a:b],tau,dt)
        else: out[a:b]=np.nan
    return out

def adelantar(x,segundos,dt):
    x=np.asarray(x,float); out=np.full(len(x),np.nan); t=np.arange(len(x))*dt
    if segundos==0: return x.copy()
    for a,b in bloques(np.isfinite(x)):
        m=(t+segundos>=t[a])&(t+segundos<=t[b-1])
        out[m]=np.interp((t+segundos)[m],t[a:b],x[a:b])
    return out

def adelantar_flags(q,segundos,dt):
    """Propaga las banderas de los dos scans usados por la interpolación temporal."""
    q=np.asarray(q,np.uint16);out=np.zeros(len(q),np.uint16)
    u=np.arange(len(q))+segundos/dt
    a=np.floor(u).astype(int);b=np.ceil(u).astype(int)
    valido=(a>=0)&(b<len(q))
    out[~valido]=8;out[valido]=q[a[valido]]|q[b[valido]]
    return out

def histeresis(v,p,s,dt):
    """AN64-3, apéndice: distingue voltaje previo original y corregido."""
    h1,h2,h3=[valor(s,k) for k in ("H1","H2","H3")]
    if h2<=0 or h3<=0: raise ValueError("H2 y H3 deben ser positivos.")
    v=np.asarray(v,float); p=np.asarray(p,float)
    if not np.isfinite(v).all() or not np.isfinite(p).all():
        raise ValueError("Histéresis necesita serie continua; hay huecos en oxígeno/presión.")
    off=valor(s,"offset"); u=v+off; out=u.copy(); c=np.exp(-dt/h3)
    d=1+h1*np.expm1(p/h2)
    if np.any(d<=0): raise ValueError("Histéresis: coeficientes/presión producen un divisor no físico.")
    for i in range(1,len(v)): out[i]=(u[i]+out[i-1]*c*d[i]-u[i-1]*c)/d[i]
    return out-off

def hampel(x,n,k,piso):
    s=pd.Series(x); med=s.rolling(n,center=True,min_periods=n).median()
    mad=s.rolling(n,center=True,min_periods=n).apply(lambda a:np.median(np.abs(a-np.median(a))),raw=True)
    return ((s-med).abs()>np.maximum(k*1.4826*mad,piso)).to_numpy()

def derivar(df,lat,lon,vars):
    df=df.copy(); p=df.presion_dbar.to_numpy();t=df.temperatura_C.to_numpy();c=df.conductividad_S_m.to_numpy()
    with np.errstate(invalid="ignore",divide="ignore",over="ignore"):
        sp=gsw.SP_from_C(c*10,t,p)  # GSW requiere mS/cm; la conversión proporciona S/m.
        sa=gsw.SA_from_SP(sp,p,lon,lat); ct=gsw.CT_from_t(sa,t,p)
        df["profundidad_m"]=-gsw.z_from_p(p,lat)
        df["salinidad_practica"]=sp; df["salinidad_absoluta_g_kg"]=sa
        df["temperatura_conservativa_C"]=ct; df["temperatura_potencial_C"]=gsw.pt0_from_t(sa,t,p)
        df["densidad_kg_m3"]=gsw.rho(sa,ct,p); df["sigma0_kg_m3"]=gsw.sigma0(sa,ct)
        df["velocidad_sonido_m_s"]=gsw.sound_speed(sa,ct,p)
    vars.update(BASE_VARIABLES)
    return df

def oxigeno43(df,s,dt,tau,ventana):
    # La histéresis ya se aplica al voltaje original, antes del alineamiento.
    kw={"soc":valor(s,"Soc"),"v_offset":valor(s,"offset"),
        "tau_20":valor(s,"Tau20") if tau else 0,
        "a":valor(s,"A"),"b":valor(s,"B"),"c":valor(s,"C"),"e":valor(s,"E"),
        "d0":0,"d1":valor(s,"D1") if tau else 0,"d2":valor(s,"D2") if tau else 0,
        "h1":0,"h2":1,"h3":1}
    co=sbcoef.Oxygen43Coefficients(**kw)
    v=df[s.senal].to_numpy();t=df.temperatura_C.to_numpy();p=df.presion_dbar.to_numpy();sp=df.salinidad_practica.to_numpy()
    out=np.full(len(v),np.nan)
    lado=max(1,int(math.floor(ventana/2/dt))) if tau else 0
    for a,b in bloques(np.isfinite(v)&np.isfinite(t)&np.isfinite(p)&np.isfinite(sp)):
        if b-a<2*lado+1: continue
        out[a:b]=sbconv.convert_sbe43_oxygen(v[a:b],t[a:b],p[a:b],sp[a:b],co,
                     apply_tau_correction=tau,apply_hysteresis_correction=False,
                     window_size=ventana,sample_interval=dt)
        if lado: out[a:a+lado]=np.nan;out[b-lado:b]=np.nan
    return out

def auxiliares(df,cfg,vars,opts,dt):
    for s in cfg.sensores[3:]:
        tag=f"s{s.indice}"
        if s.tipo=="sbe43":
            ml=oxigeno43(df,s,dt,opts["tau_oxigeno"],opts["ventana_tau_s"])
        elif s.tipo=="sbe63":
            co=sbcoef.Oxygen63Coefficients(**{k:valor(s,k) for k in sbcoef.Oxygen63Coefficients.__annotations__})
            tc=sbcoef.Thermistor63Coefficients(**{k:valor(s,k) for k in sbcoef.Thermistor63Coefficients.__annotations__})
            ml=sbconv.convert_sbe63_oxygen(df[s.senal].to_numpy(),df.sbe63_temp_V.to_numpy(),
                    df.presion_dbar.to_numpy(),df.salinidad_practica.to_numpy(),co,tc)
        elif s.tipo=="ph18":
            col=f"ph_{tag}"; s.salida=col
            df[col]=sbconv.convert_sbe18_ph(df[s.senal].to_numpy(),df.temperatura_C.to_numpy(),
                          sbcoef.PH18Coefficients(valor(s,"Slope"),valor(s,"Offset")))
            vars[col]=Variable("pH (escala de calibración del sensor)","1",f"ph{tag}","sensor")
            continue
        elif s.tipo=="eco":
            col=f"eco_{tag}";s.salida=col
            df[col]=sbconv.convert_eco(df[s.senal].to_numpy(),sbcoef.ECOCoefficients(
                valor(s,"ScaleFactor"),valor(s,"Vblank","DarkVoltage","DarkOutput")))
            nombre,unidad=opts["unidades_eco"][str(s.indice)]
            vars[col]=Variable(nombre,unidad,f"eco{tag}","sensor")
            continue
        else: continue
        col=f"oxigeno_{tag}_ml_L"; s.salida=col
        df[col]=ml;vars[col]=Variable(f"Oxígeno disuelto ({tag})","mL L⁻¹",f"sbeox{tag}ML/L","sensor")
        mg=f"oxigeno_{tag}_mg_L";um=f"oxigeno_{tag}_umol_kg"
        df[mg]=sbconv.convert_oxygen_to_mg_per_l(ml)
        df[um]=sbconv.convert_oxygen_to_umol_per_kg(ml,df.sigma0_kg_m3.to_numpy())
        vars[mg]=Variable(f"Oxígeno disuelto ({tag})","mg L⁻¹",f"sbeox{tag}Mg/L",graficar=False)
        vars[um]=Variable(f"Oxígeno disuelto ({tag})","µmol kg⁻¹",f"sbeox{tag}Mm/Kg",graficar=False)
    return df

def mascara_loop(p,dt,lat,inicio,fondo,fin,vmin):
    """Velocidad antes/después y no recontar agua hasta superar extremo previo.
    Misma lógica FIXED del toolkit; límites del lance fijados por el usuario.
    """
    z=-gsw.z_from_p(np.asarray(p,float),lat); ok=np.zeros(len(z),bool)
    vantes=np.diff(z,prepend=np.nan)/dt;vdesp=np.diff(z,append=np.nan)/dt
    for a,b,signo in ((inicio,fondo,1),(fondo+1,fin,-1)):
        previo=-np.inf
        for i in range(a,b+1):
            if not np.isfinite(z[i]): continue
            avanza=signo*z[i]>previo
            vel=(signo*vantes[i]>=vmin) and (signo*vdesp[i]>=vmin)
            if avanza and vel:
                ok[i]=True;previo=signo*z[i]
    return ok

def procesar_lance(l,opts,etapas=False):
    df=l.datos.copy(); vars=l.variables.copy();plan=l.plan;dt=l.dt
    ini,fondo,fin=[plan[k]-1 for k in ("scan_inicio","scan_fondo","scan_fin")]
    # Series conservan N y orden; las exclusiones no alteran el paso de tiempo.
    registro={}; qc={}; df["presion_dbar"]-=plan["offset_presion_dbar"]
    if plan["offset_presion_dbar"]!=0:
        df["conductividad_S_m"]=convertir_conductividad(l.raw,df.temperatura_C,df.presion_dbar,l.config.sensores[1])
    for c in l.variables:
        x=df[c].to_numpy(dtype=float,copy=True);q=np.zeros(len(df),np.uint16);q[~np.isfinite(x)]|=1
        if c in opts["rangos"]:
            lo,hi=opts["rangos"][c]; fuera=(x<lo)|(x>hi);q[fuera]|=2;x[fuera]=np.nan
        qc[c]=q;df[c]=x
    registro["01_convertido"]=df.copy() if etapas else None
    usar_histeresis=bool(plan.get("histeresis",opts.get("histeresis",False)))
    if usar_histeresis:
        for s in l.config.sensores:
            if s.tipo=="sbe43": df[s.senal]=histeresis(df[s.senal],df.presion_dbar,s,dt)
    if opts["hampel"]:
        for c,piso in (("temperatura_C",.05),("conductividad_S_m",.02),("presion_dbar",.2)):
            flags=hampel(df[c],opts["hampel_n"],opts["hampel_k"],piso);qc[c][flags]|=4
            if opts["excluir_picos"]: df.loc[flags,c]=np.nan
    for c,tau in (("presion_dbar",opts["filtro_P_s"]),("temperatura_C",opts["filtro_TC_s"]),
                  ("conductividad_S_m",opts["filtro_TC_s"])):
        df[c]=filtrar(df[c],tau,dt)
    registro["02_filtrado"]=df.copy() if etapas else None
    for c,offset in (("temperatura_C",opts["adelanto_T_s"]),("conductividad_S_m",opts["adelanto_C_s"])):
        df[c]=adelantar(df[c],offset,dt)
        qc[c]=adelantar_flags(qc[c],offset,dt)
    for s in l.config.sensores[3:]:
        off=plan["adelantos_aux_s"].get(str(s.indice),0)
        df[s.senal]=adelantar(df[s.senal],off,dt)
        qc[s.senal]=adelantar_flags(qc[s.senal],off,dt)
        if s.tipo=="sbe63": df["sbe63_temp_V"]=adelantar(df.sbe63_temp_V,off,dt)
    registro["03_alineado"]=df.copy() if etapas else None
    if opts["masa_termica"]:
        tt=df.temperatura_C.to_numpy();cc=df.conductividad_S_m.to_numpy().copy()
        for a,b in bloques(np.isfinite(tt)&np.isfinite(cc)):
            cc[a:b]=sbproc.cell_thermal_mass(tt[a:b],cc[a:b],opts["alpha"],opts["beta_inv_s"],dt)
            qc["conductividad_S_m"][a:min(b,a+math.ceil(3*opts["beta_inv_s"]/dt))]|=128
        df["conductividad_S_m"]=cc
    registro["04_masa_termica"]=df.copy() if etapas else None
    for c in qc:
        qc[c][~np.isfinite(df[c]) & (qc[c]==0)]|=8
        df[f"qc_{c}"]=qc[c]
    # Derivar antes de quitar filas: tau/histéresis mantienen continuidad temporal.
    df=derivar(df,l.meta["latitud"],l.meta["longitud"],vars)
    opts_aux=dict(opts)
    opts_aux["tau_oxigeno"]=bool(plan.get("tau_oxigeno",opts.get("tau_oxigeno",False)))
    df=auxiliares(df,l.config,vars,opts_aux,dt)
    for c in vars:
        if c not in df: continue
        qcol=f"qc_{c}"
        if qcol not in df: df[qcol]=np.where(np.isfinite(df[c]),0,1).astype(np.uint16)
        if c=="salinidad_practica" or c.startswith("oxigeno_") and c.endswith("_ml_L"):
            lo,hi=opts["rango_SP"] if c=="salinidad_practica" else opts["rango_O2_ml_L"]
            mal=(df[c]<lo)|(df[c]>hi)
            df.loc[mal,qcol]=df.loc[mal,qcol].astype(int)|2
            df.loc[mal,c]=np.nan
            dependientes=[k for k in vars if (c=="salinidad_practica" and k in (
                "salinidad_absoluta_g_kg","temperatura_conservativa_C","temperatura_potencial_C",
                "densidad_kg_m3","sigma0_kg_m3","velocidad_sonido_m_s")) or
                (c.startswith("oxigeno_") and k.startswith(c[:-4]) and k!=c)]
            if c=="salinidad_practica": dependientes += [k for k in vars if k.startswith("oxigeno_")]
            for k in dependientes:
                df.loc[mal,k]=np.nan
                qk=f"qc_{k}"
                if qk not in df: df[qk]=0
                df.loc[mal,qk]=df.loc[mal,qk].astype(int)|2
    qs=np.zeros(len(df),np.uint16);n=np.arange(len(df))
    qs[(n<ini)|(n>fin)]|=16
    qs[(df.presion_dbar<opts["presion_min_dbar"]).to_numpy()]|=64
    qs[~np.isfinite(df.presion_dbar)]|=1
    if opts["loop"]:
        ok=mascara_loop(df.presion_dbar,dt,l.meta["latitud"],ini,fondo,fin,plan["velocidad_min_m_s"])
        qs[(n>=ini)&(n<=fin)&~ok]|=32
    df["qc_scan"]=qs
    df["direccion"]=np.where((n>=ini)&(n<=fondo),"bajada",np.where((n>fondo)&(n<=fin),"subida","fuera"))
    for c in list(vars):
        if c in df: df.loc[qs!=0,c]=np.nan
    # Derive final DESPUÉS de Loop Edit/máscara para que el producto publicado
    # respete la secuencia conceptual Sea-Bird. La derivación provisional previa
    # se mantiene solo para correcciones dinámicas continuas (p. ej. SBE43).
    df=derivar(df,l.meta["latitud"],l.meta["longitud"],vars)
    # Herencia de QC para variables que dependen de T/C/P; los flags sospechosos
    # se conservan aunque el usuario haya elegido no excluirlos.
    padres=(df.qc_temperatura_C.to_numpy(np.uint16)|df.qc_conductividad_S_m.to_numpy(np.uint16)|
            df.qc_presion_dbar.to_numpy(np.uint16))
    for c,v in vars.items():
        if c in df and v.origen=="derivada" and c!="profundidad_m":
            df[f"qc_{c}"]=df[f"qc_{c}"].to_numpy(np.uint16)|padres
    for s in l.config.sensores[3:]:
        if not s.salida or s.salida not in df:continue
        origen=df[f"qc_{s.senal}"].to_numpy(np.uint16)
        if s.tipo in ("sbe43","sbe63"):origen=origen|padres
        elif s.tipo=="ph18":origen=origen|df.qc_temperatura_C.to_numpy(np.uint16)
        if s.tipo in ("sbe43","sbe63"):
            destinos=[k for k in vars if k.startswith(f"oxigeno_s{s.indice}_")]
        else:destinos=[s.salida]
        for k in destinos:df[f"qc_{k}"]=df[f"qc_{k}"].to_numpy(np.uint16)|origen
    registro["05_qc_derivados"]=df.copy() if etapas else None
    return df,vars,registro

def promediar(df,vars,paso,nmin):
    salidas=[]
    for direccion in ("bajada","subida"):
        sub=df[(df.direccion==direccion)&(df.qc_scan==0)&df.presion_dbar.notna()].copy()
        if sub.empty: continue
        sub["bin_id"]=np.floor(sub.presion_dbar/paso+.5).astype(int)
        grupos=sub.groupby("bin_id");ids=np.arange(sub.bin_id.min(),sub.bin_id.max()+1)
        res=pd.DataFrame(index=ids);res.index.name="bin_id"
        res["bin_centro_dbar"]=ids*paso;res["direccion"]=direccion
        res["n_scans"]=grupos.size().reindex(ids,fill_value=0)
        res["scan_primero"]=grupos.scan.min();res["scan_ultimo"]=grupos.scan.max()
        res["tiempo_s"]=grupos.tiempo_s.mean()
        for col in vars:
            if col not in sub: continue
            n=grupos[col].count().reindex(ids,fill_value=0)
            media=grupos[col].mean().reindex(ids)
            res[col]=media.where(n>=nmin)
            res[f"n_{col}"]=n
            res[f"sd_{col}"]=grupos[col].std().reindex(ids).where(n>=2)
            if f"qc_{col}" in sub:
                # OR conserva los motivos presentes en las observaciones que forman el bin.
                res[f"qc_{col}"]=grupos[f"qc_{col}"].agg(lambda a:int(np.bitwise_or.reduce(a.to_numpy(np.uint16)))).reindex(ids,fill_value=1)
        res["flag"]=np.where(res[["presion_dbar","temperatura_C","conductividad_S_m"]].notna().all(axis=1),0,BAD)
        res.loc[res.flag!=0,[c for c in vars if c in res]]=np.nan
        salidas.append(res.reset_index())
    return pd.concat(salidas,ignore_index=True) if salidas else pd.DataFrame()



def _qartod_init(values):
    """Flag primario QARTOD: 1=pass y 9=missing como estado inicial."""
    x = np.asarray(values, float)
    return np.where(np.isfinite(x), 1, 9).astype(np.int8)


def _qartod_set(flags, mask, code):
    """Promueve flags sin degradar FAIL ni reemplazar MISSING."""
    f = np.asarray(flags, dtype=np.int8).copy()
    m = np.asarray(mask, bool) & (f != 9)
    if int(code) == 4:
        f[m] = 4
    elif int(code) == 3:
        f[m & (f != 4)] = 3
    elif int(code) == 2:
        f[m & (f == 1)] = 2
    return f


def _qartod_aggregate(*arrays):
    """Agrega flags por fila: FAIL > SUSPECT > PASS; 9 solo si todo falta."""
    aa = [np.asarray(a, dtype=np.int8) for a in arrays if a is not None]
    if not aa:
        return np.array([], dtype=np.int8)
    A = np.vstack(aa)
    out = np.full(A.shape[1], 2, dtype=np.int8)
    all_missing = np.all(A == 9, axis=0)
    any_fail = np.any(A == 4, axis=0)
    any_sus = np.any(A == 3, axis=0)
    any_pass = np.any(A == 1, axis=0)
    out[any_pass] = 1
    out[any_sus] = 3
    out[any_fail] = 4
    out[all_missing] = 9
    return out


def _spike_suspect_mask(x, k=8.0, floor=0.0, window=7):
    """Diagnóstico robusto de spikes; solo se usa para SUSPECT, nunca borra datos."""
    x = np.asarray(x, float)
    if np.isfinite(x).sum() < max(7, window):
        return np.zeros(len(x), bool)
    n = int(window)
    n += (n % 2 == 0)
    med = pd.Series(x).rolling(n, center=True, min_periods=max(3, n//2)).median().to_numpy()
    resid = x - med
    rr = resid[np.isfinite(resid)]
    if len(rr) < 5:
        return np.zeros(len(x), bool)
    mad = float(np.nanmedian(np.abs(rr - np.nanmedian(rr))))
    scale = max(1.4826 * mad, float(floor) / max(float(k), 1.0), 1e-12)
    threshold = max(float(k) * scale, float(floor))
    return np.isfinite(resid) & (np.abs(resid) > threshold)


def aplicar_qc_publicacion_bins(bins, l, opts):
    """QC post-proceso inspirado en QARTOD/IOOS, aplicado a bins de presión.

    Mantiene intacto el bitmask QC interno. Añade flags primarios interoperables:
    1 pass, 2 not evaluated, 3 suspect/high interest, 4 fail, 9 missing.
    Los umbrales de spike/N2 son diagnósticos configurables de campaña, no
    límites universales ni una certificación GO-SHIP.
    """
    if bins.empty or not bool(AUTO_CFG.get("qc_publicacion", True)):
        return bins, pd.DataFrame()
    b = bins.copy()
    n = len(b)

    p = pd.to_numeric(b.get("presion_dbar", np.nan), errors="coerce").to_numpy(float)
    t = pd.to_numeric(b.get("temperatura_C", np.nan), errors="coerce").to_numpy(float)
    sp = pd.to_numeric(b.get("salinidad_practica", np.nan), errors="coerce").to_numpy(float)

    qp = _qartod_init(p)
    qt = _qartod_init(t)
    qs = _qartod_init(sp)

    # Gross-range FAIL: rango físico/instrumental amplio. SBE19plus V2 T=-5..35 C;
    # P se conserva amplio porque el full-scale real depende del sensor instalado.
    qp = _qartod_set(qp, np.isfinite(p) & ((p < 0.0) | (p > 11000.0)), 4)
    qt = _qartod_set(qt, np.isfinite(t) & ((t < -5.0) | (t > 35.0)), 4)
    qs = _qartod_set(qs, np.isfinite(sp) & ((sp < 0.0) | (sp > 50.0)), 4)

    # Spikes en producto ya binneado: SUSPECT solamente para no confundir
    # estructuras oceanográficas reales (frentes/termoclinas) con errores.
    qt = _qartod_set(qt, _spike_suspect_mask(
        t, AUTO_CFG["qc_spike_k_mad"], AUTO_CFG["qc_spike_floor_T_C"]), 3)
    qs = _qartod_set(qs, _spike_suspect_mask(
        sp, AUTO_CFG["qc_spike_k_mad"], AUTO_CFG["qc_spike_floor_SP"]), 3)

    # Temperatura por debajo del punto de congelación TEOS-10.
    if all(c in b for c in ("salinidad_absoluta_g_kg", "temperatura_C", "presion_dbar")):
        try:
            SA = b["salinidad_absoluta_g_kg"].to_numpy(float)
            tf = np.asarray(gsw.t_freezing(SA, p, 0), float)
            delta = tf - t
            qt = _qartod_set(qt, np.isfinite(delta) & (delta > float(AUTO_CFG["qc_freezing_suspect_C"])), 3)
            qt = _qartod_set(qt, np.isfinite(delta) & (delta > float(AUTO_CFG["qc_freezing_fail_C"])), 4)
            b["temperatura_congelacion_TEOS10_C"] = tf
        except Exception:
            pass

    # Inversiones de densidad: diagnóstico de alta sensibilidad; se marca
    # SUSPECT y se preserva el dato para revisión científica.
    b["qc_densidad_inversion"] = np.full(n, 2, dtype=np.int8)
    if all(c in b for c in ("salinidad_absoluta_g_kg", "temperatura_conservativa_C", "presion_dbar")):
        try:
            SA = b["salinidad_absoluta_g_kg"].to_numpy(float)
            CT = b["temperatura_conservativa_C"].to_numpy(float)
            good = np.isfinite(SA) & np.isfinite(CT) & np.isfinite(p)
            idx = np.flatnonzero(good)
            if len(idx) >= 4:
                order = idx[np.argsort(p[idx])]
                N2, _ = gsw.Nsquared(SA[order], CT[order], p[order], float(l.meta["latitud"]))
                N2 = np.asarray(N2, float)
                inv = np.flatnonzero(np.isfinite(N2) & (N2 < float(AUTO_CFG["qc_N2_suspect_s2"])))
                qinv = np.full(n, 1, dtype=np.int8)
                for j in inv:
                    if j < len(order)-1:
                        qinv[order[j]] = 3
                        qinv[order[j+1]] = 3
                        qt[order[j]] = 3 if qt[order[j]] == 1 else qt[order[j]]
                        qt[order[j+1]] = 3 if qt[order[j+1]] == 1 else qt[order[j+1]]
                        qs[order[j]] = 3 if qs[order[j]] == 1 else qs[order[j]]
                        qs[order[j+1]] = 3 if qs[order[j+1]] == 1 else qs[order[j+1]]
                b["qc_densidad_inversion"] = qinv
        except Exception:
            pass

    # Dominio de validez de la aproximación de 75 términos TEOS-10.
    qteos = np.full(n, 2, dtype=np.int8)
    if all(c in b for c in ("salinidad_absoluta_g_kg", "temperatura_conservativa_C", "presion_dbar")):
        try:
            SA = b["salinidad_absoluta_g_kg"].to_numpy(float)
            CT = b["temperatura_conservativa_C"].to_numpy(float)
            good = np.isfinite(SA) & np.isfinite(CT) & np.isfinite(p)
            qteos[good] = 1
            inside = np.asarray(gsw.infunnel(SA[good], CT[good], p[good]), bool)
            gi = np.flatnonzero(good)
            qteos[gi[~inside]] = 3
        except Exception:
            pass

    b["qartod_presion"] = qp
    b["qartod_temperatura"] = qt
    b["qartod_salinidad"] = qs
    b["qartod_teos10_funnel"] = qteos

    # Oxígeno: gross range + spike + estado dinámico del lance.
    o2cols = [c for c in b.columns if c.startswith("oxigeno_") and c.endswith("_ml_L")]
    qo = None
    if o2cols:
        o2 = b[o2cols[0]].to_numpy(float)
        qo = _qartod_init(o2)
        lo, hi = opts.get("rango_O2_ml_L", [0.0, 20.0])
        qo = _qartod_set(qo, np.isfinite(o2) & ((o2 < lo) | (o2 > hi)), 4)
        qo = _qartod_set(qo, _spike_suspect_mask(
            o2, AUTO_CFG["qc_spike_k_mad"], AUTO_CFG["qc_spike_floor_O2_ml_L"]), 3)
        if l.plan.get("estado_O2") == "REVISAR":
            qo = _qartod_set(qo, np.isfinite(o2), 3)
        b["qartod_oxigeno"] = qo
    else:
        b["qartod_oxigeno"] = np.full(n, 2, dtype=np.int8)

    b["qartod_resumen_ctd"] = _qartod_aggregate(qp, qt, qs, qteos)
    b["qartod_resumen"] = _qartod_aggregate(qp, qt, qs, qteos, qo)

    rows = []
    for col in ["qartod_presion", "qartod_temperatura", "qartod_salinidad",
                "qartod_oxigeno", "qartod_teos10_funnel", "qc_densidad_inversion",
                "qartod_resumen_ctd", "qartod_resumen"]:
        if col not in b:
            continue
        vals = b[col].to_numpy()
        for flag, meaning in QARTOD.items():
            count = int(np.sum(vals == flag))
            if count:
                rows.append({"variable_qc": col, "flag": flag, "meaning": meaning, "n_bins": count})
    return b, pd.DataFrame(rows)


def exportar_netcdf_cf_perfil(bins, l, ident, out, version):
    """Exporta un perfil NetCDF CF-1.13 si xarray está disponible."""
    if not bool(AUTO_CFG.get("exportar_netcdf_cf", True)) or bins.empty:
        return None, "desactivado"
    try:
        import xarray as xr
    except Exception as exc:
        return None, f"xarray no disponible: {exc}"

    sub = bins[bins.direccion == "bajada"].copy() if "direccion" in bins else bins.copy()
    sub = sub[np.isfinite(pd.to_numeric(sub["presion_dbar"], errors="coerce"))].copy()
    sub = sub.sort_values("presion_dbar")
    if sub.empty:
        return None, "sin datos de bajada"

    p = sub["presion_dbar"].to_numpy(float)
    ds = xr.Dataset(coords={"pressure": ("pressure", p)})
    ds["pressure"].attrs.update({"standard_name": "sea_water_pressure", "long_name": "Sea water pressure", "units": "dbar", "positive": "down"})
    ds = ds.assign_coords(latitude=float(l.meta["latitud"]), longitude=float(l.meta["longitud"]))
    ds["latitude"].attrs.update({"standard_name": "latitude", "units": "degrees_north"})
    ds["longitude"].attrs.update({"standard_name": "longitude", "units": "degrees_east"})

    mapping = {
        "temperatura_C": ("sea_water_temperature", "degree_Celsius"),
        "conductividad_S_m": ("sea_water_electrical_conductivity", "S m-1"),
        "salinidad_practica": ("sea_water_practical_salinity", "1"),
        "salinidad_absoluta_g_kg": ("sea_water_absolute_salinity", "g kg-1"),
        "temperatura_conservativa_C": ("sea_water_conservative_temperature", "degree_Celsius"),
        "densidad_kg_m3": ("sea_water_density", "kg m-3"),
    }
    for col, (std, unit) in mapping.items():
        if col in sub:
            ds[col] = ("pressure", sub[col].to_numpy(float))
            ds[col].attrs.update({"standard_name": std, "units": unit})

    if "profundidad_m" in sub:
        ds["depth"] = ("pressure", sub["profundidad_m"].to_numpy(float))
        ds["depth"].attrs.update({"standard_name": "depth", "units": "m", "positive": "down"})
        ds["z"] = ("pressure", -np.abs(sub["profundidad_m"].to_numpy(float)))
        ds["z"].attrs.update({"long_name": "Vertical coordinate relative to sea surface", "units": "m", "positive": "up"})

    o2u = [c for c in sub if c.startswith("oxigeno_") and c.endswith("_umol_kg")]
    if o2u:
        ds[o2u[0]] = ("pressure", sub[o2u[0]].to_numpy(float))
        ds[o2u[0]].attrs.update({"standard_name": "moles_of_oxygen_per_unit_mass_in_sea_water", "units": "umol kg-1"})

    qcols = [c for c in sub if c.startswith("qartod_")]
    qmap = {
        "temperatura_C": "qartod_temperatura",
        "salinidad_practica": "qartod_salinidad",
    }
    if o2u:
        qmap[o2u[0]] = "qartod_oxigeno"
    for q in qcols:
        ds[q] = ("pressure", sub[q].to_numpy(np.int8))
        ds[q].attrs.update({
            "long_name": q.replace("_", " "),
            "standard_name": "status_flag",
            "flag_values": np.array([1,2,3,4,9], dtype=np.int8),
            "flag_meanings": "pass not_evaluated suspect_or_high_interest fail missing",
        })
    for var, q in qmap.items():
        if var in ds and q in ds:
            ds[var].attrs["ancillary_variables"] = q

    ds.attrs.update({
        "Conventions": "CF-1.13",
        "title": f"CTD profile {l.meta['estacion']}",
        "source": "SBE 19plus/19plus V2 HEX processed with XMLCON calibration",
        "processing_software": f"CTD Crucero Auto {version}",
        "history": datetime.now(timezone.utc).isoformat() + " automated CTD processing and QC",
        "source_hex": l.ruta.name,
        "xmlcon_sha256": sha256(l.config.ruta),
        "qartod_flag_scheme": "1 pass; 2 not evaluated; 3 suspect/high interest; 4 fail; 9 missing",
        "publication_note": "Automated QC does not replace independent bottle/Winkler validation or final scientific review.",
    })
    ndir = out / "netcdf_cf"
    ndir.mkdir(exist_ok=True)
    path = ndir / f"{ident}.nc"
    ds.to_netcdf(path, engine="scipy")
    return path, "OK"


def diagnostico(l,guardar=None,plan=None):
    fig,axs=plt.subplots(1,3,figsize=(14,4.8),layout="constrained")
    for ax,c,label in zip(axs,("presion_dbar","temperatura_C","conductividad_S_m"),
                          ("Presión [dbar]","Temperatura [°C]","Conductividad [S/m]")):
        ax.plot(l.datos.scan,l.datos[c],color="#146B8B",lw=.9)
        ax.set(xlabel="Scan original (primer scan = 1)",ylabel=label)
        ax.grid(alpha=.22);ax.xaxis.set_major_locator(MaxNLocator(5));ax.tick_params(labelsize=10)
        if plan:
            for key,color in (("scan_inicio","#2C8B43"),("scan_fondo","#C85225"),("scan_fin","#763DA7")):
                ax.axvline(plan[key],color=color,lw=1.2,ls="--",label=key.replace("scan_",""))
    axs[0].invert_yaxis();axs[0].set_title('"V de la Verdad"')
    if plan: axs[0].legend(fontsize=9)
    fig.suptitle(textwrap.fill(f"{l.meta['estacion']} — {l.ruta.name}\nDiagnóstico de adquisición y selección del lance",115),fontsize=13)
    if guardar: fig.savefig(guardar,dpi=180,bbox_inches="tight")
    mostrar_fig(fig)
    return fig

def sugerir_limites(df,dt=.25):
    p=filtrar(df.presion_dbar.to_numpy(),1,dt)
    if not np.isfinite(p).any(): raise ValueError("No hay presión finita para reconocer el lance.")
    b=int(np.nanargmax(p));maxp=p[b]
    umbral=max(1.0,float(np.nanpercentile(p,2))+1.0)
    prev=np.flatnonzero(np.isfinite(p[:b])&(p[:b]<=umbral))
    a=int(prev[-1]) if len(prev) else 0
    post=np.flatnonzero(np.isfinite(p[b+1:])&(p[b+1:]<=umbral))
    fin=b+1+int(post[0]) if len(post) else len(p)-1
    if maxp-umbral<2: a=0;fin=len(p)-1
    return a+1,b+1,fin+1

def indicio_bomba(l):
    um=header_value(l.encabezado,[r"<MinimumCondFreq>\s*([0-9.]+)",r"minimum cond freq\s*=\s*([0-9.]+)"])
    de=header_value(l.encabezado,[r"<PumpDelay>\s*([0-9.]+)",r"pump delay\s*=\s*([0-9.]+)"])
    if not um or not de:return None
    mascara=l.raw.cond_Hz.to_numpy()>=float(um)
    duracion=max(1,math.ceil(float(de)/l.dt))
    for a,b in bloques(mascara):
        if b-a>duracion:return a+duracion+1
    return None

def previsualizar_oxigeno(l,opts,carpeta):
    sensores=[s for s in l.config.sensores if s.tipo in ("sbe43","sbe63")]
    for s in sensores:
        print(f"\n{l.ruta.name}: {s.etiqueta}, serie {s.serie}.")
        print("Un adelanto positivo evalúa la señal posterior en el tiempo de presión actual.")
        print("La coincidencia bajada/subida también depende del movimiento y variabilidad del agua.")
        candidatos=[0,3,4,5,6,7] if s.tipo=="sbe43" and not opts["tau_oxigeno"] else [0,.5,1,2,3,4]
        if si("¿Mostrar comparación de adelantos antes de elegir?",True):
            fig,axs=plt.subplots(2,3,figsize=(13,10),layout="constrained")
            for ax,off in zip(axs.flat,candidatos):
                l.plan["adelantos_aux_s"][str(s.indice)]=off
                test,_,_=procesar_lance(l,opts)
                z_total=-np.abs(test.profundidad_m.to_numpy(dtype=float))
                for d,color in (("bajada","#155F84"),("subida","#D17B2A")):
                    sub=test[test.direccion==d]
                    z_plot=-np.abs(sub.profundidad_m.to_numpy(dtype=float))
                    ax.plot(sub[s.salida],z_plot,label=d,color=color,lw=1)
                ax.set(xlabel="Oxígeno [mL/L]",ylabel="Coordenada vertical, z [m]",title=f"Adelanto {off:g} s")
                if np.isfinite(z_total).any():
                    ax.set_ylim(float(np.nanmin(z_total)),0.0)
                ax.grid(alpha=.2);ax.xaxis.set_major_locator(MaxNLocator(4))
            axs[0,0].legend();fig.suptitle(textwrap.fill(l.meta["estacion"]+" — selección de adelanto de oxígeno",90),fontsize=14)
            fig.savefig(carpeta/f"{slug(l.ruta.stem)}_ajuste_oxigeno_s{s.indice}.png",dpi=180,bbox_inches="tight")
            mostrar_fig(fig)
        pred=.5 if opts["tau_oxigeno"] else (4.0 if s.tipo=="sbe43" else 0.0)
        l.plan["adelantos_aux_s"][str(s.indice)]=numero("Adelanto elegido [s]; es una decisión revisable",pred,-30,60)
        l.plan.setdefault("justificacion_oxigeno",{})[str(s.indice)]=input("Motivo (ej.: comparación visual, valor validado en campaña): ").strip() or "No documentado"
        plt.close("all")

def cuestionario(lances):
    print("\nAJUSTES GENERALES. Enter acepta el valor entre corchetes; y/n confirma decisiones.")
    o={}
    o["maniobra"]=("manual","winche","mixta/desconocida")[elegir("¿Cómo se hizo la maniobra?",["A mano","Con winche","Mixta o desconocida"])-1]
    o["bomba_estandar"]=si("¿Se confirma bomba de 2000 rpm y conducto T-C estándar? Si no lo sabes, responde n",None)
    o["estabilizacion"]=input("Describe la estabilización (ej.: 3 min a 5 m; o desconocida): ").strip() or "desconocida"
    o["incidencias"]=input("Incidencias generales (pausas, fallas de bomba, cambios de montaje), o Enter: ").strip()
    tz=elegir("¿En qué hora estaba configurado el reloj del CTD?",["Desconocida: conservar fecha literal, sin llamarla UTC","UTC","Colombia (UTC−5)","Otro desfase UTC"])
    o["ctd_utc_offset_h"]=None if tz==1 else (0 if tz==2 else (-5 if tz==3 else numero("Desfase UTC en horas",0,-12,14)))
    o["zona_bitacora"]="America/Bogota (UTC-5)"
    print("Se mantendrán por separado la hora del CTD y la inmersión según bitácora.")
    print("El programa aplica coeficientes del XMLCON; no puede obtener una nueva calibración de un perfil.")
    o["calibracion_confirmada"]=si("¿El XMLCON corresponde al equipo, sensores y calibración vigentes en esta campaña?",None)
    if not o["calibracion_confirmada"]: raise ValueError("Se requiere el XMLCON correspondiente antes de convertir las señales.")
    print("Filtros Sea-Bird/INVEMAR: constantes temporales, no tamaños de ventana en scans.")
    o["filtro_P_s"]=numero("Constante de filtro de presión [s]; 0 desactiva",1,0,30)
    o["filtro_TC_s"]=numero("Constante de filtro de temperatura y conductividad [s]",.5,0,30)
    o["adelanto_T_s"]=numero("Adelanto de temperatura respecto a presión [s]",.5,-30,30)
    o["adelanto_C_s"]=numero("Adelanto de conductividad [s]; estándar SBE19plus = 0",0,-30,30)
    print("Masa térmica: α=0.04, τ=8 s para bomba 2000 rpm y conducto T-C; otros montajes necesitan otros parámetros.")
    o["masa_termica"]=si("¿Aplicar corrección de masa térmica de conductividad?",o["bomba_estandar"])
    o["alpha"]=numero("Amplitud α",.04,.00001,1) if o["masa_termica"] else 0
    o["beta_inv_s"]=numero("Constante 1/β [s]",8,.01,120) if o["masa_termica"] else 0
    print("Loop Edit marca frenadas/inversiones. 0.5 m/s (manual) o 1 m/s (winche) del manual GEO pueden excluir gran parte de un descenso lento.")
    o["loop"]=si("¿Aplicar Loop Edit?",True)
    pred=.5 if o["maniobra"]=="manual" else (1 if o["maniobra"]=="winche" else .1)
    o["velocidad_min_m_s"]=numero("Velocidad mínima inicial [m/s]; se revisará su efecto por lance",pred,0,3) if o["loop"] else 0
    o["presion_min_dbar"]=numero("Presión mínima utilizable [dbar], para excluir exposición superficial",.5,0,100)
    o["bin_dbar"]=numero("Tamaño del promedio por presión [dbar], NO metros",1,.05,100)
    o["nmin_bin"]=numero("Mínimo de observaciones válidas por variable y bin; 1 conserva bins de una observación",1,1,1000,True)
    o["guardar_subida"]=si("¿Conservar también la subida en las bases finales?",True)
    print("Se conservan todas las observaciones válidas del bin; no se trunca a 100 ni se rellenan bins vacíos.")
    o["rangos"]={"temperatura_C":[-5,45],"conductividad_S_m":[0,10],"presion_dbar":[-5,11000]}
    o["rango_SP"]=[0,50];o["rango_O2_ml_L"]=[0,20]
    print("Rangos amplios de control, editables por campaña:",o["rangos"],"SP: 0–50; O₂: 0–20 mL/L.")
    if si("¿Editar esos rangos?",False):
        for nombre,rng in list(o["rangos"].items())+[("SP",o["rango_SP"]),("O2_mL_L",o["rango_O2_ml_L"])]:
            rng[0]=numero(nombre+": mínimo",rng[0]);rng[1]=numero(nombre+": máximo",rng[1],rng[0]+1e-12)
    o["hampel"]=si("¿Detectar picos con Hampel? Puede marcar gradientes reales; requiere revisión",False)
    o["hampel_n"]=11;o["hampel_k"]=8;o["excluir_picos"]=False
    if o["hampel"]:
        o["hampel_n"]=numero("Ventana impar [scans]",11,3,501,True)
        if o["hampel_n"]%2==0:o["hampel_n"]+=1
        o["hampel_k"]=numero("Umbral de desviaciones robustas",8,3,30)
        o["excluir_picos"]=si("¿Excluir esos picos? Con n se guardan y quedan marcados",False)
    todos=[s for l in lances for s in l.config.sensores]
    hay43=any(s.tipo=="sbe43" for s in todos)
    o["tau_oxigeno"]=False;o["ventana_tau_s"]=2.;o["histeresis"]=False
    if hay43:
        print("SBE43: Tau corrige la respuesta dinámica usando dV/dt; requiere Tau20, D1, D2.")
        print("Con Tau activo, vuelve a evaluar el adelanto: 3–7 s y Tau completos pueden sobrecorregir el retardo.")
        completos=all(all(k in s.coef for k in ("tau20","d1","d2")) for s in todos if s.tipo=="sbe43")
        o["tau_oxigeno"]=si("¿Aplicar Tau al SBE43?",completos) if completos else False
        if not completos: print("Tau no disponible: faltan coeficientes en XMLCON; quedará registrado.")
        if o["tau_oxigeno"]:
            dtmax=max(l.dt for l in lances)
            o["ventana_tau_s"]=numero("Ventana para estimar dV/dt [s]",max(2,2*dtmax),2*dtmax,max(30,10*dtmax))
        hcomp=all(all(k in s.coef for k in ("h1","h2","h3")) for s in todos if s.tipo=="sbe43")
        print("Histéresis SBE43: depende de la historia de presión; recomendada por Sea-Bird para perfiles >1000 m.")
        profundo=any(l.datos.presion_dbar.max()>1000 for l in lances)
        o["histeresis"]=si("¿Aplicar histéresis con H1, H2 y H3 del XMLCON?",profundo) if hcomp else False
        if not hcomp: print("Histéresis no disponible: faltan coeficientes; no se inventarán.")
    o["unidades_eco"]={}
    for s in todos:
        if s.tipo=="eco" and str(s.indice) not in o["unidades_eco"]:
            nom="Turbidez" if "turb" in s.etiqueta.lower() else "Fluorescencia equivalente"
            unit="NTU" if "turb" in s.etiqueta.lower() else "µg L⁻¹ (equivalente de calibración)"
            print(f"\nSensor {s.etiqueta}, serie {s.serie}. Ecuación: ScaleFactor × (V − voltaje oscuro).")
            nom=input(f"Nombre de la variable [{nom}]: ").strip() or nom
            unit=input(f"Unidad según su certificado [{unit}]: ").strip() or unit
            o["unidades_eco"][str(s.indice)]=[nom,unit]
    o["etapas"]=si("¿Guardar también los datos intermedios de cada etapa?",True)
    o["mostrar_finales"]=si("¿Mostrar las figuras finales en Spyder además de guardarlas?",True)
    o["dpi"]=numero("Resolución PNG [dpi]; también se guardará PDF vectorial",400,150,600,True)
    return o

def preparar_lances(lances,opts,carpeta):
    for l in lances:
        print(f"\nREVISIÓN DE {l.meta['estacion']} / {l.ruta.name}: {len(l.datos)} scans, Δt={l.dt:g} s")
        print("La propuesta identifica el descenso principal; revisa el atemperamiento, la bomba y posibles lances repetidos.")
        a,b,c=sugerir_limites(l.datos,l.dt)
        l.plan={"scan_inicio":a,"scan_fondo":b,"scan_fin":c,"offset_presion_dbar":0.,
                "velocidad_min_m_s":opts["velocidad_min_m_s"],"adelantos_aux_s":{}}
        posible_bomba=indicio_bomba(l)
        if posible_bomba is not None:
            print(f"Indicio de bomba: umbral de conductividad + PumpDelay se cumplen hacia el scan {posible_bomba}.")
            print("Es una estimación para registros que incluyen el arranque, no un estado medido de la bomba.")
            l.plan["indicio_bomba_scan"]=posible_bomba
        while True:
            fig=diagnostico(l,plan=l.plan)
            l.plan["scan_inicio"]=numero("Primer scan del descenso ya estabilizado",l.plan["scan_inicio"],1,len(l.datos)-2,True)
            l.plan["scan_fondo"]=numero("Último scan de bajada / fondo",l.plan["scan_fondo"],l.plan["scan_inicio"]+1,len(l.datos),True)
            l.plan["scan_fin"]=numero("Último scan utilizable de subida (igual al fondo si solo hay bajada)",
                                     max(l.plan["scan_fin"],l.plan["scan_fondo"]),l.plan["scan_fondo"],len(l.datos),True)
            plt.close(fig);fig=diagnostico(l,plan=l.plan)
            if si("¿Confirmas estos límites en el gráfico?",None):
                fig.savefig(carpeta/f"{slug(l.ruta.stem)}_seleccion.png",dpi=180,bbox_inches="tight")
                plt.close(fig);break
            plt.close(fig)
        l.plan["incidencias"]=input("Observaciones de este lance, o Enter: ").strip()
        l.plan["estabilizacion_confirmada"]=si("¿El inicio elegido es posterior a la estabilización y al funcionamiento de la bomba?",None)
        if not l.plan["estabilizacion_confirmada"]:
            l.avisos.append("Estabilización/bomba no confirmadas: perfil exploratorio, requiere revisar el inicio.")
            if not si("¿Continuar como resultado exploratorio dejando esta limitación documentada?",False):
                raise ValueError("Revisa la maniobra y el scan inicial antes de continuar.")
        pmodo=elegir("Corrección adicional de cero de presión (la conversión ya resta 14.7 psi)",
                     ["Sin corrección adicional","Restar un offset conocido en dbar","Estimar mediana de un tramo confirmado EN AIRE"])
        if pmodo==2: l.plan["offset_presion_dbar"]=numero("Offset a RESTAR [dbar]; negativo si procede",0,-20,20)
        if pmodo==3:
            print("No elijas un tramo sumergido: se eliminaría presión hidrostática real.")
            a=numero("Primer scan en aire",minimo=1,maximo=len(l.datos)-1,entero=True)
            b=numero("Último scan en aire",minimo=a+1,maximo=len(l.datos),entero=True)
            z=l.datos.presion_dbar.iloc[a-1:b]
            print(f"Mediana = {z.median():.4f} dbar; desviación = {z.std():.4f} dbar.")
            if si("¿Confirmas que este tramo estuvo en aire y se debe restar su mediana?",False):
                l.plan["offset_presion_dbar"]=float(z.median());l.plan["scans_cero_presion"]=[a,b]
        for s in l.config.sensores[3:]:
            if s.tipo not in ("sbe43","sbe63","sin_uso"):
                print(f"{s.etiqueta}: sin adelanto universal para este sensor.")
                l.plan["adelantos_aux_s"][str(s.indice)]=numero("Adelanto de esta señal [s]",0,-30,60)
        # Estimar pérdidas con la presión filtrada antes de fijar el umbral.
        if opts["loop"]:
            p=filtrar(l.datos.presion_dbar-l.plan["offset_presion_dbar"],opts["filtro_P_s"],l.dt)
            a,b,c=[l.plan[k]-1 for k in ("scan_inicio","scan_fondo","scan_fin")]
            while True:
                v=l.plan["velocidad_min_m_s"]
                ok=mascara_loop(p,l.dt,l.meta["latitud"],a,b,c,v)
                prop=100*ok[a:b+1].mean()
                print(f"Con mínimo {v:g} m/s quedaría {prop:.1f}% de los scans de BAJADA, antes de otros controles.")
                if si("¿Mantener esta velocidad mínima?",None): break
                l.plan["velocidad_min_m_s"]=numero("Nuevo mínimo [m/s]",.1,0,3)
        previsualizar_oxigeno(l,opts,carpeta)

def resumen(lances,opts,carpeta):
    general=pd.DataFrame([
        ("Equipo",lances[0].config.nombre),("Filtros P / T-C [s]",f"{opts['filtro_P_s']} / {opts['filtro_TC_s']}"),
        ("Adelanto T / C [s]",f"{opts['adelanto_T_s']} / {opts['adelanto_C_s']}"),
        ("Masa térmica",f"{opts['masa_termica']}; α={opts['alpha']}, τ={opts['beta_inv_s']} s"),
        ("Oxígeno SBE43",f"Tau={opts['tau_oxigeno']}, ventana={opts['ventana_tau_s']} s; histéresis={opts['histeresis']}"),
        ("Bin",f"{opts['bin_dbar']} dbar; mínimo {opts['nmin_bin']} observaciones/variable"),
        ("Edición de picos",f"Hampel={opts['hampel']}; excluir={opts['excluir_picos']}"),
        ("TEOS-10","SP → SA → CT; σ₀ referido a 0 dbar"),
        ("Productos","CNV/CSV de trazabilidad; figuras finales PNG 300 dpi; Excel consolidado"),
        ("Zona horaria CTD",str(opts["ctd_utc_offset_h"])+" horas UTC; None = desconocida")],columns=["Ajuste","Valor"])
    tabla_visual(general,"Procesamiento propuesto — revisa antes de iniciar",carpeta/"resumen_general")
    filas=[]
    for l in lances:
        p=l.plan
        filas.append([l.ruta.name,f"{p['scan_inicio']}–{p['scan_fondo']}–{p['scan_fin']}",
                      l.dt,p["offset_presion_dbar"],p["velocidad_min_m_s"],str(p["adelantos_aux_s"])])
    detalle=pd.DataFrame(filas,columns=["HEX","Inicio / fondo / fin","Δt [s]","Offset P [dbar]","Mín. [m/s]","Adelantos auxiliares [s]"])
    tabla_visual(detalle,"Ajustes por lance",carpeta/"resumen_lances")
    if not si("¿APROBAR este plan y comenzar el procesamiento?",None): return False
    plt.close("all");return True

def guardar_cnv(path,df,vars,l,info,binario=False):
    cols=[c for c in vars if c in df]
    nombres=[(vars[c].cnv,f"{vars[c].nombre} [{vars[c].unidad}]") for c in cols]
    if "bin_centro_dbar" in df:
        cols += ["bin_centro_dbar","n_scans"];nombres += [("pbin","Centro de bin [dbar]"),("nbin","Scans por bin [1]")]
    elif "scan" in df:
        cols=["scan","tiempo_s"]+cols;nombres=[("scan","Scan original [1]"),("timeS","Tiempo transcurrido [s]")]+nombres
    vals=df[cols].to_numpy(float)
    if "flag" in df: flag=df.flag.to_numpy()
    elif "qc_scan" in df: flag=np.where(df.qc_scan==0,0,BAD)
    else: flag=np.zeros(len(df))
    vals=np.column_stack([vals,flag]);nombres.append(("flag","0 = válido; bad_flag = bin/scan inválido"))
    vals[~np.isfinite(vals)]=BAD
    lines=[f"* CTD Python {VERSION}; no generado por SBE Data Processing",
           f"* fuente_hex = {l.ruta.name}",f"* estacion = {l.meta['estacion']}",
           f"* latitude = {l.meta['latitud']}",f"* longitude = {l.meta['longitud']}",
           f"* configuracion = {l.config.ruta.name}",f"* proceso = {info}",
           f"* hora_inicial_hex = {l.hora_hex or 'desconocida'}; {l.fuente_hora}",
           "* Para parámetros, calibración, unidades y QC completos, consultar manifiesto.json y CSV.",
           "# file_type = ascii",f"# nquan = {len(nombres)}",f"# nvalues = {len(df)}",f"# bad_flag = {BAD:.8e}"]
    for i,(nombre,desc) in enumerate(nombres): lines.append(f"# name {i} = {nombre}: {desc}")
    if "bin_centro_dbar" in df: lines.append(f"# interval = decibars: {l.plan.get('bin_dbar',1)}")
    else: lines.append(f"# interval = seconds: {l.dt}")
    lines.append("*END*")
    with path.open("w",encoding="utf-8",newline="\n") as f:
        f.write("\n".join(lines)+"\n");np.savetxt(f,vals,fmt="%.10e")

def estilo_ejes(ax):
    ax.grid(alpha=.18,linewidth=.6)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    f=ScalarFormatter(useOffset=False);f.set_powerlimits((-3,4));ax.xaxis.set_major_formatter(f)
    ax.tick_params(labelsize=10,direction="out")
    ax.spines[["top","right"]].set_visible(False)

def guardar_fig(fig,base,dpi=300,mostrar=False):
    # Producto gráfico final de crucero: PNG de 300 dpi.
    fig.savefig(base.with_suffix(".png"),dpi=300,bbox_inches="tight",facecolor="white")
    if mostrar:
        mostrar_fig(fig)
        if "inline" in plt.get_backend().lower():plt.close(fig)
    else:plt.close(fig)

def figura_perfiles(df,vars,titulo,base,opts):
    sub=df[df.direccion=="bajada"]
    cols=[c for c,v in vars.items() if v.graficar and c in sub and c!="profundidad_m"]
    if not cols:return

    # SOLO PARA VISUALIZACION:
    # la base conserva profundidad_m positiva hacia abajo; en las figuras se
    # representa la coordenada vertical oceanografica z: 0 en superficie y
    # valores negativos hacia el fondo.
    z_plot=-np.abs(sub.profundidad_m.to_numpy(dtype=float))

    ncol=min(3,len(cols));nrow=math.ceil(len(cols)/ncol)
    fig,axs=plt.subplots(nrow,ncol,figsize=(4.2*ncol,4.4*nrow),squeeze=False,layout="constrained")
    colors=["#125D81","#A84525","#3A7551","#795AA6"]
    for i,(ax,col) in enumerate(zip(axs.flat,cols)):
        v=vars[col];ax.plot(sub[col],z_plot,color=colors[i%4],lw=1.35,marker=".",ms=2)
        ax.set_xlabel(textwrap.fill(v.nombre,32)+f"\n[{v.unidad}]",fontsize=11)
        ax.set_ylabel("Coordenada vertical, z [m]",fontsize=11)
        ax.set_title(f"({chr(97+i)})",loc="left",fontsize=12,weight="bold")
        if np.isfinite(z_plot).any():
            ax.set_ylim(float(np.nanmin(z_plot)),0.0)
        estilo_ejes(ax)
        if not sub[col].notna().any():ax.text(.5,.5,"Sin valores válidos",ha="center",transform=ax.transAxes)
    for ax in list(axs.flat)[len(cols):]:ax.set_visible(False)
    fig.suptitle(textwrap.fill(titulo+" — perfil de bajada",100),fontsize=14,weight="bold")
    guardar_fig(fig,base,opts["dpi"],opts["mostrar_finales"])

def figura_ts(df,color,titulo,base,opts,limites=None):
    sa="salinidad_absoluta_g_kg";ct="temperatura_conservativa_C"
    if color not in df:return False
    sub=df[(df.direccion=="bajada")&df[[sa,ct,color]].notna().all(axis=1)]
    if sub.empty:return False
    x=sub[sa].to_numpy();y=sub[ct].to_numpy();v=sub[color].to_numpy()

    # SOLO PARA VISUALIZACION:
    # profundidad_m permanece positiva en las bases. Cuando se usa como color
    # en el diagrama T-S se representa como z: 0 en superficie y negativo abajo.
    if color=="profundidad_m":
        v=-np.abs(v.astype(float))
        if limites is not None:
            zz=np.abs(np.asarray(limites,dtype=float))
            limites=(-float(np.nanmax(zz)),-float(np.nanmin(zz)))

    dx=max(np.ptp(x)*.08,.015);dy=max(np.ptp(y)*.08,.03)
    xx=np.linspace(max(0,x.min()-dx),x.max()+dx,130);yy=np.linspace(y.min()-dy,y.max()+dy,130)
    X,Y=np.meshgrid(xx,yy);sig=gsw.sigma0(X,Y)
    fig,ax=plt.subplots(figsize=(8.2,6.5),layout="constrained")
    lev=MaxNLocator(nbins=7).tick_values(float(np.nanmin(sig)),float(np.nanmax(sig)))
    lev=lev[(lev>sig.min())&(lev<sig.max())]
    if len(lev):
        cs=ax.contour(X,Y,sig,levels=lev,colors="#757D85",linewidths=.65,alpha=.75)
        ax.clabel(cs,inline=True,fmt="%.2f",fontsize=8)
    kw={"vmin":limites[0],"vmax":limites[1]} if limites else {}
    sc=ax.scatter(x,y,c=v,cmap="viridis" if color=="profundidad_m" else "cividis",s=18,linewidths=0,rasterized=True,**kw)
    cb=fig.colorbar(sc,ax=ax,pad=.025,shrink=.93)
    cb.set_label("Coordenada vertical, z [m]" if color=="profundidad_m" else "Oxígeno disuelto [mL L⁻¹]",fontsize=11)
    ax.set(xlabel="Salinidad absoluta, $S_A$ [g kg$^{-1}$]",ylabel=r"Temperatura conservativa, $\Theta$ [°C]",
           xlim=(xx[0],xx[-1]),ylim=(yy[0],yy[-1]))
    ax.set_title(textwrap.fill(titulo,72)+"\n"+r"TEOS-10; isolíneas de $\sigma_0$ [kg m$^{-3}$]",fontsize=12,pad=12)
    estilo_ejes(ax);ax.xaxis.set_major_locator(MaxNLocator(5));ax.yaxis.set_major_locator(MaxNLocator(6))
    guardar_fig(fig,base,opts["dpi"],opts["mostrar_finales"])
    return True

def sha256(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()

def json_default(x):
    if isinstance(x,(np.integer,np.floating,np.bool_)):return x.item()
    if isinstance(x,(Path,datetime)):return str(x)
    raise TypeError(type(x).__name__)

def json_dump(p,obj):
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=json_default),encoding="utf-8")

def ejecutar(lances,opts,out,tabla_path,metadatos):
    for sub in ("00_entradas","datos","cnv","figuras","etapas","control_calidad","control_calidad_publicacion"):
        (out/sub).mkdir(exist_ok=True)
    shutil.copy2(tabla_path,out/"00_entradas"/tabla_path.name)
    for l in lances:
        destino=out/"00_entradas"/(sha256(l.config.ruta)[:10]+"_"+l.config.ruta.name)
        if not destino.exists():shutil.copy2(l.config.ruta,destino)
    metadatos.to_csv(out/"00_entradas"/"estaciones_asociadas.csv",index=False,encoding="utf-8-sig")
    manifiesto={"version":VERSION,"inicio_utc":datetime.now(timezone.utc).isoformat(),"ajustes":opts,
        "fuentes":FUENTES,"flags_internos":QC,"python":sys.version,
        "paquetes":{p:importlib.metadata.version(p) for p in ("numpy","pandas","scipy","matplotlib","gsw","seabirdscientific")},
        "entrada_tabla_sha256":sha256(tabla_path),"lances":[],
        "bin_definicion":"Centro k*Δp; intervalos [centro-Δp/2, centro+Δp/2). Medias por variable, sin interpolación ni truncamiento.",
        "estado":"EN_PROCESO"}
    json_dump(out/"manifiesto.json",manifiesto)
    bases=[];resultados=[];vars_total={};errores=[]
    for i,l in enumerate(lances,1):
        ident=f"{i:03d}_{slug(l.meta['estacion'])}_{slug(l.ruta.stem)}"
        print(f"\nProcesando {i}/{len(lances)}: {l.ruta.name}")
        item={"id":ident,"hex":str(l.ruta),"sha256_hex":sha256(l.ruta),"sha256_xmlcon":sha256(l.config.ruta),
              "metadatos":l.meta,"hora_hex":l.hora_hex,"fuente_hora":l.fuente_hora,"dt_s":l.dt,
              "plan":l.plan,"avisos":l.avisos,"sensores":[asdict(s) for s in l.config.sensores]}
        try:
            l.plan["bin_dbar"]=opts["bin_dbar"]
            df,vars,etapas=procesar_lance(l,opts,opts["etapas"])
            bins=promediar(df,vars,opts["bin_dbar"],opts["nmin_bin"])
            if bins.empty:raise ValueError("No quedan bins. Revisa límites, filtros, calibraciones y velocidad mínima.")
            if not opts["guardar_subida"]:bins=bins[bins.direccion=="bajada"].copy()
            bins, qc_pub = aplicar_qc_publicacion_bins(bins, l, opts)
            if not qc_pub.empty:
                qc_pub.to_csv(out/"control_calidad_publicacion"/f"{ident}_qartod_conteos.csv",index=False,encoding="utf-8-sig")
            buenos=bins[bins.direccion=="bajada"][["salinidad_absoluta_g_kg","temperatura_conservativa_C"]].dropna()
            if buenos.empty:raise ValueError("No hay datos T-S válidos de bajada; este lance no se presentará como procesado.")
            for c,v in l.meta.items():
                if c not in bins and c not in ("archivo_xmlcon",):bins[c]=v
            bins["id_perfil"]=ident;bins["inicio_hex_literal"]=l.hora_hex
            if l.hora_hex and opts["ctd_utc_offset_h"] is not None:
                start=pd.Timestamp(l.hora_hex)
                if start.tzinfo is None:start=start.tz_localize(timezone(timedelta(hours=opts["ctd_utc_offset_h"])))
                bins["fecha_hora_utc_media_bin"]=(start.tz_convert("UTC")+pd.to_timedelta(bins.tiempo_s,unit="s")).astype(str)
            l.raw.to_csv(out/"datos"/f"{ident}_senales.csv",index=False,encoding="utf-8-sig")
            (out/"datos"/f"{ident}_encabezado_hex.txt").write_text(l.encabezado,encoding="utf-8")
            df.to_csv(out/"datos"/f"{ident}_scans_procesados_qc.csv",index=False,encoding="utf-8-sig")
            bins.to_csv(out/"datos"/f"{ident}_perfiles.csv",index=False,encoding="utf-8-sig")
            nc_path, nc_status = exportar_netcdf_cf_perfil(bins, l, ident, out, VERSION)
            for direc,sub in bins.groupby("direccion"):
                guardar_cnv(out/"cnv"/f"{ident}_{direc}.cnv",sub,vars,l,"promedio por presión; "+direc)
            if opts["etapas"]:
                for etapa,datos in etapas.items():
                    if datos is not None:datos.to_csv(out/"etapas"/f"{ident}_{etapa}.csv",index=False,encoding="utf-8-sig")
            det=[]
            for c in ("qc_scan",)+tuple(k for k in df if k.startswith("qc_") and k!="qc_scan"):
                for bit,motivo in QC.items():
                    n=int(((df[c].to_numpy(dtype=np.uint16)&bit)!=0).sum())
                    if n:det.append([c,bit,motivo,n])
            pd.DataFrame(det,columns=["variable_qc","bit","motivo","n_scans"]).to_csv(
                out/"control_calidad"/f"{ident}_conteos_qc.csv",index=False,encoding="utf-8-sig")
            item.update(
                estado=l.plan.get("estado_operativo", "PROCESADO_QC_AUTOMATICO"),
                estado_CTD=l.plan.get("estado_CTD", "NO_EVALUADO"),
                estado_O2=l.plan.get("estado_O2", "NO_EVALUADO"),
                estado_publicacion="PROVISIONAL_HASTA_VALIDACION_EXTERNA_Y_REVISION_CIENTIFICA",
                netcdf_cf=str(nc_path) if nc_path is not None else "",
                netcdf_cf_estado=nc_status,
                n_scans=len(df), n_bins=len(bins), n_bins_ts_bajada=len(buenos),
                fraccion_scans_aceptados=float((df.qc_scan==0).mean())
            )
            bases.append(bins);resultados.append((l,bins,vars,ident));vars_total.update(vars)
        except Exception as exc:
            item.update(estado="ERROR",error=str(exc));errores.append([l.ruta.name,str(exc)])
            (out/"control_calidad"/f"{ident}_error.txt").write_text(traceback.format_exc(),encoding="utf-8")
            print("ERROR en este lance:",exc)
        manifiesto["lances"].append(item);json_dump(out/"manifiesto.json",manifiesto)
    if not bases:
        manifiesto["estado"]="SIN_RESULTADOS";json_dump(out/"manifiesto.json",manifiesto)
        raise ValueError("Ningún lance produjo datos válidos. Consulta control_calidad y manifiesto.json.")
    base=pd.concat(bases,ignore_index=True,sort=False)
    base.to_csv(out/"base_ctd_procesada.csv",index=False,encoding="utf-8-sig")
    dic=pd.DataFrame([{"columna":k,**asdict(v)} for k,v in vars_total.items()])
    dic.to_csv(out/"diccionario_variables.csv",index=False,encoding="utf-8-sig")
    qc_all=pd.DataFrame([{k:v for k,v in it.items() if k in ("id","hex","estado","estado_CTD","estado_O2","estado_publicacion","netcdf_cf_estado","n_scans","n_bins","n_bins_ts_bajada","error")} for it in manifiesto["lances"]])
    qc_all.to_csv(out/"resumen_procesamiento.csv",index=False,encoding="utf-8-sig")
    with sqlite3.connect(out/"base_ctd_procesada.sqlite") as con:
        base.to_sql("perfiles",con,index=False,if_exists="replace")
        dic.to_sql("variables",con,index=False,if_exists="replace")
        qc_all.to_sql("control_perfiles",con,index=False,if_exists="replace")
        pd.DataFrame([{"json":json.dumps(manifiesto,ensure_ascii=False,default=json_default)}]).to_sql("procedencia",con,index=False,if_exists="replace")
    z=base.loc[base.direccion=="bajada","profundidad_m"].dropna()
    zlim=(float(z.min()),float(z.max())) if len(z) else None
    for l,bins,vars,ident in resultados:
        titulo=f"{l.meta['estacion']} | {l.meta['latitud']:.5f}°, {l.meta['longitud']:.5f}°"
        try:
            figura_perfiles(bins,vars,titulo,out/"figuras"/f"{ident}_panel_variables",opts)
            figura_ts(bins,"profundidad_m",l.meta["estacion"],out/"figuras"/f"{ident}_TS_profundidad",opts,zlim)
            ocols=[k for k in vars if k.startswith("oxigeno_") and k.endswith("_ml_L")]
            for col in ocols:
                valores=base.loc[base.direccion=="bajada",col].dropna()
                clim=(float(valores.min()),float(valores.max())) if len(valores) else None
                figura_ts(bins,col,l.meta["estacion"],out/"figuras"/f"{ident}_TS_{col}",opts,clim)
            if not ocols:print(f"{l.ruta.name}: no hay sensor de oxígeno convertible; se omite T–S con oxígeno.")
        except Exception as exc:
            errores.append([l.ruta.name,"Figuras: "+str(exc)])
            print("No se completaron todas las figuras de",l.ruta.name,exc)
    try:
        figura_ts(base,"profundidad_m","Campaña — perfiles de bajada",out/"figuras"/"campana_TS_profundidad",opts,zlim)
        for col in [c for c in vars_total if c.startswith("oxigeno_") and c.endswith("_ml_L")]:
            figura_ts(base,col,"Campaña — perfiles de bajada",out/"figuras"/f"campana_TS_{col}",opts)
    except Exception as exc:errores.append(["campaña","Figuras: "+str(exc)])
    manifiesto["estado"]="PARCIAL_CON_ERRORES" if errores else "COMPLETO_QC_AUTOMATICO_PROVISIONAL"
    manifiesto["fin_utc"]=datetime.now(timezone.utc).isoformat();manifiesto["errores"]=errores
    json_dump(out/"manifiesto.json",manifiesto)
    with sqlite3.connect(out/"base_ctd_procesada.sqlite") as con:
        pd.DataFrame([{"json":json.dumps(manifiesto,ensure_ascii=False,default=json_default)}]).to_sql("procedencia",con,index=False,if_exists="replace")
    tabla_visual(qc_all,"Resultado de la ejecución",mostrar=True)
    print(f"\nEstado: {manifiesto['estado']}\nResultados: {out}\nBase CSV y SQLite generadas.")
    if errores:print("Hay errores documentados; no todos los lances/productos se completaron.")
    return base,qc_all,manifiesto

def main_manual_legacy():
    """Modo manual heredado de v1. No se usa en el flujo automático de crucero."""
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.labelsize":11,
                         "axes.titlesize":12,"pdf.fonttype":42,"ps.fonttype":42,"savefig.facecolor":"white"})
    print(textwrap.dedent("""
    ASISTENTE DE PROCESAMIENTO CTD PARA SPYDER
    Necesitas:
      1. Excel (.xlsx, primera hoja) o CSV de estaciones. Columnas:
         archivo_hex | estacion | latitud | longitud
         Opcionales: fecha_bitacora | hora_inmersion
         Ejemplo: CTD_001.hex | Estación 01/Nocturno | 12.002967 | -73.662794
         Fecha día/mes/año; coordenadas decimales WGS84; oeste negativo.
         Si falta archivo_hex, se pedirá asociar cada fila manualmente.
      2. Carpeta con archivos originales .hex (incluye subcarpetas).
      3. Configuración .xmlcon correspondiente al equipo y campaña.
      4. Carpeta donde guardar resultados; conocimiento de la maniobra y
         disponibilidad para revisar gráficas y confirmar parámetros.

    El XMLCON identifica sensores: no es necesario conocerlos de antemano.
    Soporte HEX de esta versión: SBE19plus/V2, modo perfil, presión strain gauge.
    Sensores no convertibles se informan; nunca se adivinan sus ecuaciones.
    Una tabla de estaciones CSV es diferente de un archivo oceanográfico CNV.
    El procesamiento usa calibraciones existentes; no sustituye calibración
    de laboratorio, validación con muestras ni revisión científica.
    El modo automático guarda Excel consolidado, PNG 300 dpi, CNV/CSV de trazabilidad y control de calidad.
    Ctrl+C cancela sin modificar originales. Responde y/n; Enter acepta valores.
    """))
    if not si("¿Leíste los requisitos y quieres continuar?",None):return None
    tabla_path=ruta_input("Selecciona el Excel/CSV de estaciones",extensiones=(".xlsx",".xlsm",".csv"))
    df=leer_tabla(tabla_path)
    hex_dir=ruta_input("Selecciona la carpeta con los HEX",carpeta=True)
    xml_path=ruta_input("Selecciona el XMLCON de esta campaña",extensiones=(".xmlcon",))
    df,rutas=vincular(df,hex_dir);cfg=leer_config(xml_path)
    if importlib.metadata.version("seabirdscientific")!="2.7.10":
        print("La versión del toolkit difiere de 2.7.10; sus funciones podrían haber cambiado.")
        if not si("¿Continuar y registrar la versión instalada?",False):return None
    invent=pd.DataFrame([[s.indice,s.etiqueta,s.serie,s.fecha_calibracion,s.tipo,
                         s.motivo or "Ecuación identificada"] for s in cfg.sensores],
                       columns=["Canal XML","Sensor","Serie","Calibración","Tipo","Conversión"])
    tabla_visual(invent,"Sensores identificados en el XMLCON")
    desconocidos=[s for s in cfg.sensores if s.tipo=="desconocido"]
    if desconocidos and not si("Hay sensores sin conversión: ¿continuar conservando sus voltajes como señales sin convertir?",False):return None
    print("\nVerificando todos los HEX; ninguna fila corrupta será descartada silenciosamente.")
    lances=[]
    for i,p in enumerate(rutas):
        header,raw=leer_hex(p,cfg);datos,vars=convertir_base(raw,cfg)
        hora,fuente=fecha_header(header)
        avisos=cfg.avisos.copy()
        if not hora:avisos.append(fuente)
        l=Lance(p,df.iloc[i].to_dict(),cfg,header,raw,datos,vars,cfg.navg/4,hora,fuente,avisos)
        lances.append(l)
        print(f"  {p.name}: {len(raw)} scans, Δt={l.dt:g} s; inicio: {hora or 'no identificado'} ({fuente}).")
    opts=cuestionario(lances)
    padre=ruta_input("Selecciona la carpeta de salida",carpeta=True)
    out=padre/("CTD_procesado_"+datetime.now().strftime("%Y%m%d_%H%M%S_%f"));out.mkdir()
    diag=out/"diagnosticos";diag.mkdir()
    invent.to_csv(diag/"sensores_detectados.csv",index=False,encoding="utf-8-sig")
    while True:
        preparar_lances(lances,opts,diag)
        if resumen(lances,opts,diag):break
        if not si("¿Volver al cuestionario para modificar el plan?",True):return None
        opts=cuestionario(lances)
    json_dump(out/"plan_confirmado.json",{"ajustes":opts,"lances":{l.ruta.name:l.plan for l in lances}})
    return ejecutar(lances,opts,out,tabla_path,df)


# =============================================================================
# MODO AUTOMÁTICO DE CRUCERO - v2
# =============================================================================
# Filosofía: inferir -> cuantificar -> procesar -> revisar SOLO excepciones.
# Los coeficientes de calibración absoluta SIEMPRE provienen del XMLCON/certificado.
# Los datos del propio perfil se usan para decisiones dinámicas (cast, loops, lag, etc.).

AUTO_CFG = {
    # Rutas: None abre selector. En crucero conviene escribirlas una sola vez aquí.
    "tabla_estaciones": None,
    "carpeta_hex": None,
    "xmlcon": None,
    "carpeta_salida": None,

    # Productos
    "bin_dbar": 1.0,
    "presion_min_dbar": 0.5,
    "nmin_bin": 1,
    "guardar_subida_final": False,
    "guardar_intermedios": True,
    "mostrar_solo_revision": True,
    "dpi": 300,

    # Metadatos/configuración que NO deben inferirse de los perfiles.
    # None conserva la hora literal del CTD sin llamarla UTC.
    "ctd_utc_offset_h": None,
    # True para bomba/ducto T-C estándar de SBE19plus/V2 en agua salada.
    "montaje_TC_estandar": True,

    # Operación en crucero. True = procesa solo HEX nuevos/cambiados y consolida.
    "modo_incremental": False, # Por defecto: nueva carpeta por corrida; nunca mezcla campañas
    "nombre_sesion": "CTD_CRUCERO_SESION_V23",
    "modo_entrada": None, # None: preguntar; 1 completa, 2 faltan HEX, 3 sin tabla
    "tolerancia_min": 15, # SOLO nivel de revisión, nunca excluye candidatos
    "hora_ctd_offset_h": None, # Hora de reloj del CTD respecto a UTC, p.ej. -5
    "hora_bitacora_offset_h": -5, # Colombia; CAMBIAR para otros países


    # Umbrales de revisión. No cambian datos: solo deciden si un lance requiere ojo humano.
    "confianza_min_cast": 0.72,
    "retencion_min_loop": 0.80,
    "confianza_min_o2": 0.60,

    # Optimizaciones
    "optimizar_C_campana": True,
    "optimizar_O2": True,
    "auto_cero_presion": True,
    "auto_histeresis_profunda": True,
    "umbral_histeresis_dbar": 1000.0,

    # Optimización T-C. Si la función objetivo es plana o el mínimo cae en el
    # borde superior, se conserva el preset y se documenta la baja identificabilidad.
    "C_lag_min_s": 0.0,
    "C_lag_max_s": 0.15,
    "C_lag_step_s": 0.025,
    "C_mejora_min_rel": 0.01,
    "C_min_lances_campana": 3,

    # SBE43. Se explora un intervalo suficientemente amplio y se penalizan
    # soluciones en el borde, funciones objetivo planas y desacuerdo con la campaña.
    "O2_lag_step_s": 0.25,
    "O2_lag_max_s_sin_tau": 10.0,
    "O2_lag_max_s_con_tau": 4.0,
    "O2_mejora_tau_min_rel": 0.05,
    "O2_mejora_min_rel": 0.01,
    "O2_tolerancia_parsimonia_rel": 0.01,
    "O2_min_lances_campana": 3,
    # QC operacional del ajuste dinámico O2. No son umbrales universales de QARTOD:
    # controlan únicamente cuándo un lance se manda a revisión científica.
    "O2_penalizacion_optimo_revisar_rel": 0.15,
    "O2_empeora_vs_cero_revisar_rel": 0.03,

    # QC post-proceso sobre bins de 1 dbar. Los umbrales son conservadores y
    # editables por campaña; los tests diagnósticos marcan SUSPECT, no eliminan datos.
    "qc_publicacion": True,
    "qc_spike_k_mad": 8.0,
    "qc_spike_floor_T_C": 0.20,
    "qc_spike_floor_SP": 0.10,
    "qc_spike_floor_O2_ml_L": 0.15,
    "qc_N2_suspect_s2": -1e-5,
    "qc_freezing_suspect_C": 0.05,
    "qc_freezing_fail_C": 0.20,

    # Validación externa opcional. Para afirmaciones de exactitud absoluta de alta
    # calidad se requieren referencias independientes (salinómetro/Winkler).
    "tabla_referencias": None,
    "referencia_tolerancia_dbar": 2.0,

    # Salidas interoperables. NetCDF CF es opcional: si xarray no está instalado,
    # el procesamiento continúa y la auditoría lo documenta como no generado.
    "exportar_netcdf_cf": True,

    # Identificador para no reutilizar silenciosamente parámetros de campaña
    # estimados por una versión anterior del algoritmo.
    "algoritmo_dinamico_version": "2.2",
}


def _rmed(x, n):
    x = np.asarray(x, float)
    n = max(3, int(n))
    n += (n % 2 == 0)
    return pd.Series(x).rolling(n, center=True, min_periods=max(2, n // 3)).median().to_numpy()


def _rscale(x, floor=1e-12):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 4:
        return 1.0
    q = np.nanpercentile(x, [25, 75])
    s = (q[1] - q[0]) / 1.349
    if not np.isfinite(s) or s < floor:
        s = np.nanstd(x)
    return float(max(s if np.isfinite(s) else 1.0, floor))


def _persistent_first(mask, nmin, desde=0, hasta=None):
    mask = np.asarray(mask, bool).copy()
    if hasta is None:
        hasta = len(mask)
    mask[:max(0, desde)] = False
    mask[max(0, hasta):] = False
    for a, b in bloques(mask):
        if b - a >= nmin:
            return a
    return None


def opciones_auto(lances):
    """Ajustes de campaña que no requieren preguntar lance por lance."""
    marino = []
    for l in lances:
        c = l.datos.conductividad_S_m.to_numpy(float)
        marino.append(np.nanmedian(c) > 0.5)
    masa = bool(np.mean(marino) >= 0.5 and AUTO_CFG["montaje_TC_estandar"])

    o = {
        "maniobra": "inferida_desde_presion",
        "bomba_estandar": bool(AUTO_CFG["montaje_TC_estandar"]),
        "estabilizacion": "detección automática por presión + régimen hidráulico",
        "incidencias": "",
        "ctd_utc_offset_h": AUTO_CFG["ctd_utc_offset_h"],
        "zona_bitacora": "America/Bogota (UTC-5)",
        "calibracion_confirmada": True,
        "filtro_P_s": 1.0,
        "filtro_TC_s": 0.5,
        "adelanto_T_s": 0.5,
        "adelanto_C_s": 0.0,
        "masa_termica": masa,
        "alpha": 0.04 if masa else 0.0,
        "beta_inv_s": 8.0 if masa else 0.0,
        "loop": True,
        "velocidad_min_m_s": 0.0,   # se infiere por lance
        "presion_min_dbar": float(AUTO_CFG["presion_min_dbar"]),
        "bin_dbar": float(AUTO_CFG["bin_dbar"]),
        "nmin_bin": int(AUTO_CFG["nmin_bin"]),
        "guardar_subida": bool(AUTO_CFG["guardar_subida_final"]),
        "rangos": {"temperatura_C": [-5, 45], "conductividad_S_m": [0, 10], "presion_dbar": [-5, 11000]},
        "rango_SP": [0, 50],
        "rango_O2_ml_L": [0, 20],
        "hampel": False,
        "hampel_n": 11,
        "hampel_k": 8,
        "excluir_picos": False,
        # Se pueden sobreescribir por lance en plan.
        "tau_oxigeno": False,
        "ventana_tau_s": 2.0,
        "histeresis": False,
        "unidades_eco": {},
        "etapas": bool(AUTO_CFG["guardar_intermedios"]),
        "mostrar_finales": False,
        "dpi": 300,
    }

    todos = [s for l in lances for s in l.config.sensores]
    for s in todos:
        if s.tipo == "eco" and str(s.indice) not in o["unidades_eco"]:
            nom = "Turbidez" if "turb" in s.etiqueta.lower() else "Fluorescencia equivalente"
            unit = "NTU" if "turb" in s.etiqueta.lower() else "unidad de calibración"
            o["unidades_eco"][str(s.indice)] = [nom, unit]
    return o


def estimar_cero_presion(l):
    """Estima offset SOLO con evidencia fuerte de tramo en aire al inicio."""
    if not AUTO_CFG["auto_cero_presion"]:
        return 0.0, 0.0, "desactivado"

    um = header_value(
        l.encabezado,
        [r"<MinimumCondFreq>\s*([0-9.]+)", r"minimum cond freq\s*=\s*([0-9.]+)"]
    )
    if not um or "cond_Hz" not in l.raw:
        return 0.0, 0.25, "sin evidencia independiente de aire"

    cond = l.raw.cond_Hz.to_numpy(float)
    p = l.datos.presion_dbar.to_numpy(float)
    seco = np.isfinite(cond) & (cond < float(um))
    nmin = max(8, int(round(5 / l.dt)))
    bloques_secos = [(a, b) for a, b in bloques(seco) if a == 0 and b - a >= nmin]

    if not bloques_secos:
        return 0.0, 0.30, "sin tramo seco inicial persistente"

    a, b = bloques_secos[0]
    b = min(b, a + max(nmin, int(round(60 / l.dt))))
    z = p[a:b]
    z = z[np.isfinite(z)]
    if len(z) < nmin:
        return 0.0, 0.30, "tramo de aire insuficiente"

    med = float(np.nanmedian(z))
    sd = float(np.nanstd(z))
    span = float(np.nanmax(z) - np.nanmin(z))

    if abs(med) <= 3 and sd <= 0.15 and span <= 0.6:
        conf = float(np.clip(0.75 + 0.2 * (1 - min(sd / 0.15, 1)), 0, 0.98))
        return med, conf, f"aire confirmado por conductividad; sd={sd:.3f} dbar"

    return 0.0, 0.35, f"tramo seco no estable: med={med:.3f}, sd={sd:.3f}"


def inferir_cast(l, opts):
    """Detecta soak/downcast/fondo/upcast con presión, persistencia y bomba inferida."""
    offset, conf0, mot0 = estimar_cero_presion(l)
    p0 = l.datos.presion_dbar.to_numpy(float) - offset
    p = filtrar(p0, opts["filtro_P_s"], l.dt)

    if np.isfinite(p).sum() < 20:
        raise ValueError("Presión insuficiente para inferencia automática.")

    bottom = int(np.nanargmax(p))
    pmax = float(p[bottom])
    if pmax < max(0.75, float(opts["presion_min_dbar"]) + 0.35):
        raise ValueError(f"Excursión vertical insuficiente para perfil: {pmax:.2f} dbar")

    nvel = max(3, int(round(2 / l.dt)))
    ps = _rmed(p, nvel)
    dp = np.gradient(ps, l.dt)
    dp = _rmed(dp, nvel)

    pos = dp[:bottom][np.isfinite(dp[:bottom]) & (dp[:bottom] > 0)]
    med_dp = float(np.nanmedian(pos)) if len(pos) else 0.05

    w = max(6, int(round(8 / l.dt)))
    trend = np.full(len(p), np.nan)
    trend[w:] = ps[w:] - ps[:-w]
    frac_pos = pd.Series(dp > 0).rolling(w, min_periods=max(2, w // 2)).mean().to_numpy()
    delta_req = max(0.25, min(0.7, 0.18 * pmax), 0.30 * max(med_dp, 0.02) * w * l.dt)
    desc = (trend > delta_req) & (frac_pos >= 0.70) & (p > opts["presion_min_dbar"])

    pump = indicio_bomba(l)
    desde = max(0, pump - 1 if pump is not None else 0)
    cand = _persistent_first(desc, max(3, int(round(4 / l.dt))), desde=desde, hasta=bottom)

    if cand is None:
        a, _, _ = sugerir_limites(l.datos, l.dt)
        start = a - 1
        fallback = True
    else:
        start = max(desde, cand - w)
        fallback = False
        while start < bottom and (not np.isfinite(p[start]) or p[start] < opts["presion_min_dbar"]):
            start += 1

    surf = max(opts["presion_min_dbar"], min(2.0, 0.02 * pmax))
    low = np.isfinite(p) & (p <= surf)
    post = _persistent_first(low, max(3, int(round(5 / l.dt))), desde=bottom + 1)
    end = int(post) if post is not None else len(p) - 1
    if end <= bottom:
        end = len(p) - 1

    sec = dp[start:bottom + 1]
    monot = float(np.mean(sec[np.isfinite(sec)] > 0)) if np.isfinite(sec).any() else 0
    pre = dp[max(0, start - w):start]
    quiet = float(np.nanmedian(np.abs(pre))) if np.isfinite(pre).any() else np.nan
    quiet_score = 1.0 if np.isfinite(quiet) and quiet < max(0.08, 0.25 * med_dp) else 0.45
    pump_score = 0.75 if pump is None else (1.0 if pump - 1 <= start else 0.20)
    return_score = 1.0 if post is not None else 0.55
    depth_score = min(1.0, pmax / 50.0)

    conf = (
        0.10 * conf0 +
        0.35 * monot +
        0.15 * quiet_score +
        0.15 * pump_score +
        0.10 * return_score +
        0.15 * depth_score
    )
    if fallback:
        conf *= 0.75
    conf = float(np.clip(conf, 0, 0.99))

    return {
        "scan_inicio": int(start + 1),
        "scan_fondo": int(bottom + 1),
        "scan_fin": int(end + 1),
        "offset_presion_dbar": float(offset),
        "velocidad_min_m_s": 0.0,
        "adelantos_aux_s": {},
        "indicio_bomba_scan": int(pump) if pump is not None else None,
        "pmax_dbar": pmax,
        "confianza_cast": conf,
        "monotonia_down": monot,
        "cero_presion_confianza": conf0,
        "cero_presion_metodo": mot0,
        "deteccion_fallback": fallback,
    }


def inferir_loop(l, opts, plan):
    p = filtrar(
        l.datos.presion_dbar.to_numpy(float) - plan["offset_presion_dbar"],
        opts["filtro_P_s"], l.dt
    )
    a, b, c = [plan[k] - 1 for k in ("scan_inicio", "scan_fondo", "scan_fin")]
    z = -gsw.z_from_p(p, float(l.meta["latitud"]))
    v = np.gradient(_rmed(z, max(3, int(round(1 / l.dt)))), l.dt)

    sec = v[a:b + 1]
    pos = sec[np.isfinite(sec) & (sec > 0)]
    med = float(np.nanmedian(pos)) if len(pos) else 0.1
    q10 = float(np.nanpercentile(pos, 10)) if len(pos) >= 10 else med * 0.1

    candidatos = []
    for frac in (0.25, 0.20, 0.15, 0.10, 0.05):
        candidatos.append(max(0.02, min(frac * med, max(q10, 0.02))))
    candidatos = sorted(set(round(x, 4) for x in candidatos), reverse=True)

    elegido = candidatos[-1] if candidatos else 0.02
    ret = 0.0
    for vv in candidatos:
        ok = mascara_loop(p, l.dt, l.meta["latitud"], a, b, c, vv)
        r = float(ok[a:b + 1].mean()) if b >= a else 0
        if r >= 0.80:
            elegido, ret = vv, r
            break
        if r > ret:
            elegido, ret = vv, r

    plan["velocidad_min_m_s"] = float(elegido)
    plan["velocidad_mediana_down_m_s"] = med
    plan["retencion_loop"] = ret
    return plan


def _savgol_safe(x, window, poly=2):
    """Suavizado Savitzky-Golay robusto frente a NaN; conserva NaN originales."""
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    if ok.sum() < 7:
        return x.copy()
    idx = np.arange(len(x))
    xf = x.copy()
    xf[~ok] = np.interp(idx[~ok], idx[ok], x[ok])
    w = max(5, int(window))
    if w % 2 == 0:
        w += 1
    if w >= len(xf):
        w = len(xf) - 1 if len(xf) % 2 == 0 else len(xf)
    if w < 5:
        return x.copy()
    poly = min(int(poly), w - 2)
    y = signal.savgol_filter(xf, window_length=w, polyorder=poly, mode="interp")
    y[~ok] = np.nan
    return y


def _score_C_lag(l, opts, plan, lag):
    """
    Función objetivo para el pequeño alineamiento temporal de conductividad.

    Un T-C mal alineado genera anomalías de salinidad de alta frecuencia en
    regiones de fuerte gradiente térmico. La versión 2.0 usaba una mediana
    móvil centrada que podía devolver residuo exactamente cero en perfiles
    monótonos. Aquí se cuantifican: (1) residuo respecto a Savitzky-Golay,
    (2) ruido de la derivada de SP y (3) correlación del residuo con dT/dt.
    """
    p = filtrar(
        l.datos.presion_dbar.to_numpy(float) - plan["offset_presion_dbar"],
        opts["filtro_P_s"], l.dt
    )
    t = adelantar(
        filtrar(l.datos.temperatura_C.to_numpy(float), opts["filtro_TC_s"], l.dt),
        opts["adelanto_T_s"], l.dt
    )
    c = adelantar(
        filtrar(l.datos.conductividad_S_m.to_numpy(float), opts["filtro_TC_s"], l.dt),
        lag, l.dt
    )

    if opts["masa_termica"]:
        cc = c.copy()
        for aa, bb in bloques(np.isfinite(t) & np.isfinite(c)):
            if bb - aa >= 4:
                cc[aa:bb] = sbproc.cell_thermal_mass(
                    t[aa:bb], c[aa:bb], opts["alpha"], opts["beta_inv_s"], l.dt
                )
        c = cc

    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        sp = gsw.SP_from_C(c * 10, t, p)

    a, b = plan["scan_inicio"] - 1, plan["scan_fondo"] - 1
    ss = np.asarray(sp[a:b + 1], float)
    tt = np.asarray(t[a:b + 1], float)
    pp = np.asarray(p[a:b + 1], float)

    good = np.isfinite(ss) & np.isfinite(tt) & np.isfinite(pp)
    good &= pp > max(2.0, float(opts["presion_min_dbar"]))
    if good.sum() < 120:
        return np.nan

    # Mantiene el orden temporal; interpola únicamente para derivadas/suavizado.
    idx = np.arange(len(ss))
    def fill(x):
        x = np.asarray(x, float)
        ok = np.isfinite(x) & good
        if ok.sum() < 5:
            return np.full_like(x, np.nan)
        y = x.copy()
        y[~ok] = np.interp(idx[~ok], idx[ok], x[ok])
        return y

    sf = fill(ss); tf = fill(tt); pf = fill(pp)
    if not np.isfinite(sf).any():
        return np.nan

    # Tendencia de 5 s y derivada suavizada de 8 s.
    w5 = max(7, int(round(5.0 / l.dt)))
    w8 = max(9, int(round(8.0 / l.dt)))
    strend = _savgol_safe(sf, w5, 2)
    resid = sf - strend

    dsdt = np.gradient(sf, l.dt)
    dsdt_s = _savgol_safe(dsdt, w8, 2)
    dtdt = np.gradient(tf, l.dt)
    dpdt = np.gradient(pf, l.dt)

    # Se concentra donde el T-C alignment es observable: fuerte gradiente de T
    # y movimiento descendente real. Excluye bordes del filtro/alineamiento.
    abs_dt = np.abs(dtdt)
    lim_t = np.nanpercentile(abs_dt[good], 60)
    med_dp = np.nanmedian(dpdt[good & (dpdt > 0)])
    min_dp = max(0.01, 0.05 * med_dp) if np.isfinite(med_dp) else 0.01
    edge = max(w8, int(round(3.0 / l.dt)))
    interior = np.zeros(len(ss), bool)
    if len(ss) > 2 * edge:
        interior[edge:-edge] = True
    else:
        interior[:] = True
    focus = good & interior & (abs_dt >= lim_t) & (dpdt > min_dp)
    if focus.sum() < 30:
        focus = good & interior & (dpdt > min_dp)
    if focus.sum() < 20:
        return np.nan

    sscale = _rscale(sf[good], 1e-4)
    spike = np.nanmedian(np.abs(resid[focus])) / sscale
    deriv_noise = np.nanmedian(np.abs(dsdt[focus] - dsdt_s[focus])) * l.dt / sscale
    second = np.gradient(dsdt, l.dt)
    rough = np.nanmedian(np.abs(second[focus])) * (l.dt ** 2) / sscale

    rr = resid[focus]
    gg = dtdt[focus]
    corr = 0.0
    finite_corr = np.isfinite(rr) & np.isfinite(gg)
    if finite_corr.sum() >= 20 and np.nanstd(rr[finite_corr]) > 0 and np.nanstd(gg[finite_corr]) > 0:
        corr = abs(float(np.corrcoef(rr[finite_corr], gg[finite_corr])[0, 1]))

    return float(spike + 0.45 * deriv_noise + 0.10 * rough + 0.05 * corr)


def optimizar_C_campana(lances, opts, diag_prev=None):
    """
    Estima un único pequeño adelanto de C para la campaña. Si la señal objetiva
    no identifica un mínimo defendible, conserva el preset (0 s) y lo documenta.
    """
    rows = []
    if opts.get("adelanto_C_fijado", False):
        opts["C_metodo"] = "parametro_campana_reutilizado"
        opts["C_confianza"] = float(opts.get("C_confianza", 1.0))
        return float(opts.get("adelanto_C_s", 0.0)), pd.DataFrame()
    if not AUTO_CFG["optimizar_C_campana"]:
        opts["C_metodo"] = "optimizacion_desactivada"
        opts["C_confianza"] = 1.0
        return float(opts.get("adelanto_C_s", 0.0)), pd.DataFrame()

    cand = np.arange(
        float(AUTO_CFG["C_lag_min_s"]),
        float(AUTO_CFG["C_lag_max_s"]) + float(AUTO_CFG["C_lag_step_s"]) / 2,
        float(AUTO_CFG["C_lag_step_s"])
    )

    cast_info = []
    for l in lances:
        if l.plan.get("confianza_cast", 0) < 0.65 or l.plan.get("pmax_dbar", 0) < 30:
            continue
        vals = []
        for lag in cand:
            score = _score_C_lag(l, opts, l.plan, float(lag))
            vals.append(score)
            rows.append({"hex": l.ruta.name, "lag_C_s": float(lag), "score": score,
                         "algorithm_version": AUTO_CFG["algoritmo_dinamico_version"]})

        arr = np.asarray(vals, float)
        valid = np.isfinite(arr)
        if valid.sum() < max(4, len(cand) // 2):
            continue
        i = int(np.nanargmin(arr))
        i0 = int(np.nanargmin(np.abs(cand - 0.0)))
        base = arr[i0] if np.isfinite(arr[i0]) else np.nanmedian(arr[valid])
        best = arr[i]
        rel_range = (np.nanmax(arr[valid]) - best) / max(abs(best), 1e-12)
        imp = (base - best) / max(abs(base), 1e-12) if np.isfinite(base) else 0.0
        boundary = i == len(cand) - 1
        flat = rel_range < 0.005
        cast_info.append({
            "hex": l.ruta.name, "best_lag": float(cand[i]), "best_score": float(best),
            "improvement": float(imp), "boundary": bool(boundary), "flat": bool(flat)
        })

    diag_new = pd.DataFrame(rows)
    prev = pd.DataFrame() if diag_prev is None else diag_prev.copy()
    if not prev.empty:
        if "algorithm_version" in prev:
            prev = prev[prev.algorithm_version.astype(str) == str(AUTO_CFG["algoritmo_dinamico_version"])]
        else:
            prev = pd.DataFrame()
    diag = pd.concat([prev, diag_new], ignore_index=True) if not prev.empty else diag_new
    if not diag.empty:
        diag = diag.sort_values(["hex", "lag_C_s"]).drop_duplicates(["hex", "lag_C_s"], keep="last")
    if diag.empty:
        opts.update(C_metodo="sin_datos_suficientes_preset", C_confianza=0.20, C_n_lances=0)
        return float(opts.get("adelanto_C_s", 0.0)), diag

    # Normaliza cada curva por score a lag=0 para que todos los lances pesen
    # de forma comparable independientemente de su rango hidrográfico.
    curves = []
    for h, g in diag.groupby("hex"):
        g = g.copy().sort_values("lag_C_s")
        base_rows = g[np.isclose(g.lag_C_s, 0.0)]
        if base_rows.empty or not np.isfinite(base_rows.score.iloc[0]) or base_rows.score.iloc[0] <= 0:
            continue
        base = float(base_rows.score.iloc[0])
        g["score_rel"] = g.score / base
        curves.append(g)
    if not curves:
        opts.update(C_metodo="objetivo_no_identificable_preset", C_confianza=0.20, C_n_lances=0)
        return 0.0, diag

    dn = pd.concat(curves, ignore_index=True)
    camp = dn.groupby("lag_C_s").agg(
        score_rel_med=("score_rel", "median"),
        score_rel_iqr=("score_rel", lambda x: float(np.nanpercentile(x, 75) - np.nanpercentile(x, 25))),
        n_lances=("hex", "nunique")
    ).reset_index()
    for _, r in camp.iterrows():
        diag.loc[np.isclose(diag.lag_C_s, r.lag_C_s), "score_rel_camp_med"] = r.score_rel_med

    ncasts = int(dn.hex.nunique())
    bestrow = camp.sort_values(["score_rel_med", "lag_C_s"]).iloc[0]
    raw_best_lag = float(bestrow.lag_C_s)
    best_rel = float(bestrow.score_rel_med)
    improvement = max(0.0, 1.0 - best_rel)
    boundary = np.isclose(raw_best_lag, cand[-1])
    near = camp[camp.score_rel_med <= best_rel * 1.01 + 1e-12]
    chosen = float(near.sort_values("lag_C_s").iloc[0].lag_C_s)
    near_frac = len(near) / max(len(camp), 1)
    median_rel = float(np.nanmedian(camp.score_rel_med))
    contrast = max(0.0, (median_rel - best_rel) / max(abs(median_rel), 1e-12))
    evidence = max(improvement, contrast)

    identifiable = (
        ncasts >= int(AUTO_CFG["C_min_lances_campana"]) and
        evidence >= float(AUTO_CFG["C_mejora_min_rel"]) and
        not boundary and near_frac <= 0.60
    )
    if not identifiable:
        lag = 0.0
        metodo = "objetivo_plano_o_borde_conserva_preset"
    else:
        lag = chosen
        metodo = "minimo_robusto_campana"

    conf = 0.20
    conf += 0.30 * min(ncasts / max(int(AUTO_CFG["C_min_lances_campana"]), 1), 1.0)
    conf += 0.30 * min(evidence / 0.05, 1.0)
    conf += 0.20 * (1.0 - min(near_frac / 0.60, 1.0))
    if boundary:
        conf *= 0.45
    if not identifiable:
        conf = min(conf, 0.49)

    opts["C_metodo"] = metodo
    opts["C_confianza"] = float(np.clip(conf, 0, 0.98))
    opts["C_n_lances"] = ncasts
    opts["C_mejora_rel"] = improvement
    opts["C_contraste_rel"] = contrast
    opts["C_minimo_borde"] = bool(boundary)
    opts["C_lag_objetivo_s"] = raw_best_lag

    lag = float(np.clip(round(lag / float(AUTO_CFG["C_lag_step_s"])) * float(AUTO_CFG["C_lag_step_s"]),
                        float(AUTO_CFG["C_lag_min_s"]), float(AUTO_CFG["C_lag_max_s"])))
    return lag, diag


def _core_o2(l, opts, plan, s, hyst):
    p = l.datos.presion_dbar.to_numpy(float) - plan["offset_presion_dbar"]
    t = l.datos.temperatura_C.to_numpy(float)
    c = l.datos.conductividad_S_m.to_numpy(float)
    v = l.datos[s.senal].to_numpy(float)

    if hyst:
        v = histeresis(v, p, s, l.dt)

    p = filtrar(p, opts["filtro_P_s"], l.dt)
    t = filtrar(t, opts["filtro_TC_s"], l.dt)
    c = filtrar(c, opts["filtro_TC_s"], l.dt)
    t = adelantar(t, opts["adelanto_T_s"], l.dt)
    c = adelantar(c, opts["adelanto_C_s"], l.dt)

    if opts["masa_termica"]:
        cc = c.copy()
        for a, b in bloques(np.isfinite(t) & np.isfinite(c)):
            if b - a >= 4:
                cc[a:b] = sbproc.cell_thermal_mass(
                    t[a:b], c[a:b], opts["alpha"], opts["beta_inv_s"], l.dt
                )
        c = cc

    sp = gsw.SP_from_C(c * 10, t, p)
    return p, t, sp, v


def _curve_bin(p, y, a, b, paso=2.0):
    pp = np.asarray(p[a:b], float)
    yy = np.asarray(y[a:b], float)
    ok = np.isfinite(pp) & np.isfinite(yy) & (pp >= 2)
    if ok.sum() < 10:
        return pd.Series(dtype=float)

    bid = np.floor(pp[ok] / paso + 0.5).astype(int)
    s = pd.DataFrame({"b": bid, "y": yy[ok]}).groupby("b").y.median()
    s.index = s.index.astype(float) * paso
    return s


def _score_o2(p, o2, plan):
    a, b, c = plan["scan_inicio"] - 1, plan["scan_fondo"] - 1, plan["scan_fin"] - 1
    down = _curve_bin(p, o2, a, b + 1, 2.0)
    up = _curve_bin(p, o2, b, c + 1, 2.0)
    common = down.index.intersection(up.index)

    if len(common) < 20:
        return np.nan, len(common)

    yd = down.loc[common].to_numpy(float)
    yu = up.loc[common].to_numpy(float)
    scale = _rscale(np.r_[yd, yu], 0.02)
    amp = np.nanmedian(np.abs(yd - yu)) / scale

    x = common.to_numpy(float)
    gd = np.gradient(yd, x)
    gu = np.gradient(yu, x)
    shape = np.nanmedian(np.abs(gd - gu)) / _rscale(np.r_[gd, gu], 1e-4)

    sec = np.diff(yd, n=2)
    rough = np.nanmedian(np.abs(sec)) / _rscale(yd, 0.02) if len(sec) else 0
    invalid = np.mean((np.r_[yd, yu] < 0) | (np.r_[yd, yu] > 20))

    return float(amp + 0.35 * shape + 0.05 * rough + 3 * invalid), len(common)


def _resumen_objetivo_o2(subset, lag_max):
    """Resume identificabilidad de una curva score(lag) para un estado Tau dado."""
    subset = subset[np.isfinite(subset.score)].sort_values("lag_s").copy()
    if subset.empty:
        return None
    raw = subset.sort_values(["score", "lag_s"]).iloc[0]
    best_score = float(raw.score)
    tol = float(AUTO_CFG["O2_tolerancia_parsimonia_rel"])
    near = subset[subset.score <= best_score * (1.0 + tol) + 1e-12]
    chosen = near.sort_values("lag_s").iloc[0]
    zero = subset[np.isclose(subset.lag_s, 0.0)]
    score0 = float(zero.score.iloc[0]) if not zero.empty else float(np.nanmedian(subset.score))
    improvement = (score0 - best_score) / max(abs(score0), 1e-12)
    boundary = float(raw.lag_s) >= float(lag_max) - 0.5 * float(AUTO_CFG["O2_lag_step_s"])
    near_frac = len(near) / max(len(subset), 1)
    q50 = float(np.nanmedian(subset.score))
    contrast = (q50 - best_score) / max(abs(q50), 1e-12)
    return {
        "row": chosen,
        "raw_best_lag": float(raw.lag_s),
        "best_score": best_score,
        "score0": score0,
        "improvement": float(improvement),
        "boundary": bool(boundary),
        "near_frac": float(near_frac),
        "contrast": float(contrast),
    }


def optimizar_o2_lance(l, opts, plan):
    """
    Evalúa la función objetivo del SBE43 en un rango amplio. La selección final
    puede ser reemplazada después por un parámetro robusto de campaña.
    """
    sensores = [s for s in l.config.sensores if s.tipo == "sbe43"]
    if not sensores:
        return pd.DataFrame()

    s = sensores[0]
    hkeys = all(k in s.coef for k in ("h1", "h2", "h3"))
    hyst = bool(
        AUTO_CFG["auto_histeresis_profunda"] and hkeys and
        plan["pmax_dbar"] >= AUTO_CFG["umbral_histeresis_dbar"]
    )
    taukeys = all(k in s.coef for k in ("tau20", "d1", "d2"))

    rows = []
    step = max(float(AUTO_CFG["O2_lag_step_s"]), float(l.dt))
    estados_tau = [False, True] if taukeys else [False]
    for tau in estados_tau:
        try:
            p, t, sp, v = _core_o2(l, opts, plan, s, hyst)
        except Exception as exc:
            plan.setdefault("avisos_auto", []).append("O2 core: " + str(exc))
            continue

        hi = float(AUTO_CFG["O2_lag_max_s_con_tau"] if tau else AUTO_CFG["O2_lag_max_s_sin_tau"])
        for lag in np.arange(0, hi + step / 2, step):
            va = adelantar(v, float(lag), l.dt)
            tmp = pd.DataFrame({
                s.senal: va,
                "temperatura_C": t,
                "presion_dbar": p,
                "salinidad_practica": sp
            })
            try:
                ml = oxigeno43(tmp, s, l.dt, tau, 2.0)
                score, n = _score_o2(p, ml, plan)
            except Exception:
                score, n = np.nan, 0

            rows.append({
                "hex": l.ruta.name,
                "sensor": s.serie,
                "tau": bool(tau),
                "histeresis": hyst,
                "lag_s": float(lag),
                "score": score,
                "common_bins": int(n),
                "algorithm_version": AUTO_CFG["algoritmo_dinamico_version"]
            })

    diag = pd.DataFrame(rows)
    valid = diag[np.isfinite(diag.score)].copy() if not diag.empty else pd.DataFrame()

    plan["histeresis"] = hyst
    if valid.empty:
        plan["tau_oxigeno"] = False
        plan["adelantos_aux_s"][str(s.indice)] = 0.0
        plan["o2_score"] = np.nan
        plan["o2_confianza"] = 0.10
        plan["o2_metodo"] = "sin_objetivo_valido"
        plan["o2_minimo_borde"] = True
        return diag

    r_no = _resumen_objetivo_o2(valid[~valid.tau], float(AUTO_CFG["O2_lag_max_s_sin_tau"])) if (~valid.tau).any() else None
    r_tau = _resumen_objetivo_o2(valid[valid.tau], float(AUTO_CFG["O2_lag_max_s_con_tau"])) if valid.tau.any() else None

    elegido = r_no if r_no is not None else r_tau
    tau = False if r_no is not None else True
    if r_tau is not None and (r_no is None or r_tau["best_score"] < (1.0 - float(AUTO_CFG["O2_mejora_tau_min_rel"])) * r_no["best_score"]):
        elegido = r_tau
        tau = True

    row = elegido["row"]
    lag = float(row.lag_s)
    # La evidencia puede ser una mejora respecto a lag=0 o, si el mínimo real
    # está en 0 s, un contraste claro frente al resto de candidatos.
    evidencia = max(float(elegido["improvement"]), float(elegido["contrast"]))
    if evidencia < float(AUTO_CFG["O2_mejora_min_rel"]):
        lag = 0.0
        metodo = "objetivo_poco_identificable_conserva_0s"
    elif np.isclose(float(elegido["raw_best_lag"]), 0.0):
        lag = 0.0
        metodo = "minimo_objetivo_individual_0s"
    else:
        metodo = "minimo_objetivo_individual"

    bins = float(row.common_bins)
    conf = 0.20
    conf += 0.25 * min(bins / 100.0, 1.0)
    conf += 0.25 * min(max(evidencia, 0.0) / 0.10, 1.0)
    conf += 0.20 * min(max(elegido["contrast"], 0.0) / 0.10, 1.0)
    conf += 0.10 * (1.0 - min(elegido["near_frac"] / 0.50, 1.0))
    if elegido["boundary"]:
        conf *= 0.45
    conf = float(np.clip(conf, 0, 0.98))

    plan["tau_oxigeno"] = bool(tau)
    plan["adelantos_aux_s"][str(s.indice)] = lag
    plan["o2_score"] = float(row.score)
    plan["o2_confianza"] = conf
    plan["o2_metodo"] = metodo
    plan["o2_minimo_borde"] = bool(elegido["boundary"])
    plan["o2_lag_individual_s"] = float(row.lag_s)
    plan["o2_lag_raw_min_s"] = float(elegido["raw_best_lag"])
    plan["o2_mejora_rel"] = float(elegido["improvement"])
    plan["o2_evidencia_rel"] = float(evidencia)
    plan["o2_near_frac"] = float(elegido["near_frac"])

    if len(sensores) > 1:
        plan.setdefault("avisos_auto", []).append("Más de un SBE43: revisar decisión Tau individual.")
    return diag


def optimizar_o2_campana(lances, opts, carpeta, diag_prev=None):
    """
    Combina las curvas objetivo de lances de alta calidad y, cuando hay soporte
    suficiente, usa un único Tau/lag por sensor para evitar sobreajuste por perfil.
    """
    diags = []
    for l in lances:
        d = optimizar_o2_lance(l, opts, l.plan) if AUTO_CFG["optimizar_O2"] else pd.DataFrame()
        if not d.empty:
            d.to_csv(carpeta / f"{slug(l.ruta.stem)}_objetivo_O2.csv", index=False, encoding="utf-8-sig")
            diags.append(d)

    if not diags:
        return pd.DataFrame(), pd.DataFrame()

    current_diag = pd.concat(diags, ignore_index=True)
    prev = pd.DataFrame() if diag_prev is None else diag_prev.copy()
    if not prev.empty:
        if "algorithm_version" in prev:
            prev = prev[prev.algorithm_version.astype(str) == str(AUTO_CFG["algoritmo_dinamico_version"])]
        else:
            prev = pd.DataFrame()
    all_diag = pd.concat([prev, current_diag], ignore_index=True) if not prev.empty else current_diag
    all_diag = all_diag.sort_values(["hex", "sensor", "tau", "lag_s"]).drop_duplicates(
        ["hex", "sensor", "tau", "lag_s"], keep="last"
    )
    camp_rows = []

    for sensor, ds in all_diag.groupby("sensor"):
        # Referencia común por lance = score sin Tau y lag=0. Así Tau y no-Tau
        # quedan en la misma escala relativa.
        norm_parts = []
        for h, g in ds.groupby("hex"):
            base = g[(~g.tau) & np.isclose(g.lag_s, 0.0) & np.isfinite(g.score)]
            if base.empty or float(base.score.iloc[0]) <= 0:
                continue
            gg = g[np.isfinite(g.score)].copy()
            gg["score_rel"] = gg.score / float(base.score.iloc[0])
            norm_parts.append(gg)
        if not norm_parts:
            continue
        nd = pd.concat(norm_parts, ignore_index=True)
        curve = nd.groupby(["tau", "lag_s"]).agg(
            score_rel_med=("score_rel", "median"),
            score_rel_iqr=("score_rel", lambda x: float(np.nanpercentile(x,75)-np.nanpercentile(x,25))),
            n_lances=("hex", "nunique")
        ).reset_index()
        curve.insert(0, "sensor", sensor)
        camp_rows.append(curve)

        ncasts = int(nd.hex.nunique())
        if ncasts < int(AUTO_CFG["O2_min_lances_campana"]) and not opts.get("O2_campana_fijado", False):
            continue

        summaries = {}
        for tau in sorted(curve.tau.unique()):
            sub = curve[curve.tau == tau].sort_values("lag_s")
            raw = sub.sort_values(["score_rel_med", "lag_s"]).iloc[0]
            best = float(raw.score_rel_med)
            near = sub[sub.score_rel_med <= best * (1.0 + float(AUTO_CFG["O2_tolerancia_parsimonia_rel"])) + 1e-12]
            chosen = near.sort_values("lag_s").iloc[0]
            lag_max = float(AUTO_CFG["O2_lag_max_s_con_tau"] if tau else AUTO_CFG["O2_lag_max_s_sin_tau"])
            summaries[bool(tau)] = {
                "lag": float(chosen.lag_s), "raw_lag": float(raw.lag_s), "score": best,
                "boundary": float(raw.lag_s) >= lag_max - 0.5 * float(AUTO_CFG["O2_lag_step_s"]),
                "near_frac": len(near) / max(len(sub),1)
            }

        # Si ya existe un parámetro estable aprendido con esta versión, úsalo
        # y deja que las nuevas curvas sirvan para detectar lances anómalos.
        if opts.get("O2_campana_fijado", False) and str(opts.get("O2_sensor", sensor)) == str(sensor):
            pick_tau = bool(opts.get("O2_tau_campana", False))
            fixed_lag = float(opts.get("O2_lag_campana_s", 0.0))
            subfix = curve[curve.tau == pick_tau]
            if not subfix.empty:
                ix = (subfix.lag_s - fixed_lag).abs().idxmin()
                rfix = subfix.loc[ix]
                chosen = {"lag": float(rfix.lag_s), "raw_lag": float(rfix.lag_s),
                          "score": float(rfix.score_rel_med), "boundary": False, "near_frac": 0.0}
            else:
                chosen = None
        else:
            pick_tau = False
            chosen = summaries.get(False)
            if True in summaries and (chosen is None or summaries[True]["score"] < (1.0 - float(AUTO_CFG["O2_mejora_tau_min_rel"])) * chosen["score"]):
                pick_tau = True
                chosen = summaries[True]
        if chosen is None:
            continue

        improvement = max(0.0, 1.0 - float(chosen["score"]))
        selected_curve = curve[curve.tau == pick_tau]
        medrel = float(np.nanmedian(selected_curve.score_rel_med)) if not selected_curve.empty else float(chosen["score"])
        contrast = max(0.0, (medrel - float(chosen["score"])) / max(abs(medrel), 1e-12))
        evidencia = max(improvement, contrast)
        # Si la campaña no identifica un mínimo defendible, conserva 0 s.
        if (not opts.get("O2_campana_fijado", False) and
                evidencia < float(AUTO_CFG["O2_mejora_min_rel"])):
            chosen = dict(chosen)
            chosen["lag"] = 0.0

        if opts.get("O2_campana_fijado", False):
            conf = float(opts.get("O2_confianza_campana", 1.0))
        else:
            conf = 0.20
            conf += 0.30 * min(ncasts / max(int(AUTO_CFG["O2_min_lances_campana"]),1), 1.0)
            conf += 0.30 * min(evidencia / 0.10, 1.0)
            conf += 0.20 * (1.0 - min(chosen["near_frac"] / 0.50, 1.0))
            if chosen["boundary"]:
                conf *= 0.45
            if evidencia < float(AUTO_CFG["O2_mejora_min_rel"]):
                conf = min(conf, 0.44)
        conf = float(np.clip(conf, 0, 0.98))

        # Aplica el consenso solo a lances que usan este sensor. La histéresis
        # sigue siendo individual porque depende de la profundidad alcanzada.
        # La compatibilidad del lance se evalúa DENTRO del Tau seleccionado:
        #   G = (J0-Jcamp)/J0 : mejora respecto a no adelantar.
        #   P = Jcamp/Jmin-1 : penalización respecto al óptimo individual.
        for l in lances:
            sensors_l = [s for s in l.config.sensores if s.tipo == "sbe43" and s.serie == str(sensor)]
            if not sensors_l:
                continue
            s = sensors_l[0]
            d_tau = ds[(ds.hex == l.ruta.name) & (ds.tau == pick_tau) & np.isfinite(ds.score)].copy()
            dl = d_tau[np.isclose(d_tau.lag_s, chosen["lag"])]
            d0 = d_tau[np.isclose(d_tau.lag_s, 0.0)]
            best_indiv_tau = float(d_tau.score.min()) if not d_tau.empty else np.nan
            score0 = float(d0.score.iloc[0]) if not d0.empty else np.nan
            common_score = float(dl.score.iloc[0]) if not dl.empty else np.nan
            penalty = (common_score / best_indiv_tau - 1.0
                       if np.isfinite(common_score) and np.isfinite(best_indiv_tau) and best_indiv_tau > 0 else np.nan)
            gain0 = ((score0 - common_score) / score0
                     if np.isfinite(common_score) and np.isfinite(score0) and score0 > 0 else np.nan)
            ratio = common_score / best_indiv_tau if np.isfinite(penalty) else np.nan
            compat = 1.0
            if np.isfinite(penalty):
                compat = float(np.clip(1.0 - max(0.0, penalty) / 0.30, 0.0, 1.0))
            if np.isfinite(gain0) and gain0 < 0:
                compat *= float(np.clip(1.0 + gain0 / 0.10, 0.0, 1.0))
            l.plan["tau_oxigeno"] = bool(pick_tau)
            l.plan["adelantos_aux_s"][str(s.indice)] = float(chosen["lag"])
            l.plan["o2_lag_campana_s"] = float(chosen["lag"])
            l.plan["o2_score_campana_en_lance"] = common_score
            l.plan["o2_score_cero_lance"] = score0
            l.plan["o2_score_optimo_lance_tau"] = best_indiv_tau
            l.plan["o2_ratio_campana_vs_optimo_lance"] = ratio
            l.plan["o2_penalizacion_optimo_rel"] = penalty
            l.plan["o2_mejora_vs_cero_rel"] = gain0
            l.plan["o2_minimo_borde_campana"] = bool(chosen["boundary"])
            l.plan["o2_metodo"] = "consenso_robusto_campana"
            l.plan["o2_confianza"] = float(min(conf, 0.50 + 0.48 * compat))
            if np.isfinite(common_score):
                l.plan["o2_score"] = common_score

        # Expone el parámetro de campaña para persistencia incremental.
        opts["O2_sensor"] = str(sensor)
        opts["O2_tau_campana"] = bool(pick_tau)
        opts["O2_lag_campana_s"] = float(chosen["lag"])
        opts["O2_confianza_campana"] = float(conf)
        opts["O2_evidencia_campana_rel"] = float(evidencia)
        opts["O2_campana_estable"] = bool(
            opts.get("O2_campana_fijado", False) or
            (ncasts >= int(AUTO_CFG["O2_min_lances_campana"]) and
             conf >= 0.60 and not chosen["boundary"])
        )

    camp = pd.concat(camp_rows, ignore_index=True) if camp_rows else pd.DataFrame()
    if not camp.empty:
        camp.to_csv(carpeta / "diagnostico_O2_campana.csv", index=False, encoding="utf-8-sig")
    return all_diag, camp


def guardar_diag_auto(l, carpeta):
    p = l.plan
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.8), layout="constrained")
    for ax, c, label in zip(
        axs,
        ("presion_dbar", "temperatura_C", "conductividad_S_m"),
        ("Presión [dbar]", "Temperatura [°C]", "Conductividad [S/m]")
    ):
        ax.plot(l.datos.scan, l.datos[c], lw=.9)
        for key, lab in (
            ("scan_inicio", "inicio auto"),
            ("scan_fondo", "fondo auto"),
            ("scan_fin", "fin auto")
        ):
            ax.axvline(p[key], ls="--", lw=1.1, label=lab)
        ax.set(xlabel="Scan", ylabel=label)
        ax.grid(alpha=.2)

    axs[0].invert_yaxis()
    axs[0].legend(fontsize=8)
    fig.suptitle(
        f"{l.meta['estacion']} — inferencia automática | confianza del lance={p['confianza_cast']:.2f}"
    )
    fig.savefig(
        carpeta / f"{slug(l.ruta.stem)}_diagnostico_auto.png",
        dpi=180, bbox_inches="tight"
    )
    if AUTO_CFG["mostrar_solo_revision"] and p.get("estado_auto") == "REVISAR":
        mostrar_fig(fig)
    else:
        plt.close(fig)


def preparar_auto(lances, opts, carpeta):
    carpeta.mkdir(parents=True, exist_ok=True)

    # 1. Geometría del lance + Loop Edit adaptativo.
    for l in lances:
        l.plan = inferir_cast(l, opts)
        inferir_loop(l, opts, l.plan)

    # 2. Ajuste C común de campaña. Acumula evidencia entre corridas incrementales.
    c_acum_path = carpeta / "diagnostico_adelanto_C_acumulado.csv"
    c_prev = pd.read_csv(c_acum_path) if c_acum_path.exists() else pd.DataFrame()
    lagC, diagC = optimizar_C_campana(lances, opts, c_prev)
    opts["adelanto_C_s"] = lagC
    if not diagC.empty:
        diagC.to_csv(c_acum_path, index=False, encoding="utf-8-sig")
        diagC.to_csv(carpeta / "diagnostico_adelanto_C_campana.csv", index=False, encoding="utf-8-sig")

    # 3. O2: curvas por lance + consenso robusto de campaña. También acumula
    # evidencia entre corridas para que el modo tiempo real aprenda después de
    # varios lances aunque se procese uno por vez.
    o2_acum_path = carpeta / "diagnostico_O2_acumulado.csv"
    o2_prev = pd.read_csv(o2_acum_path) if o2_acum_path.exists() else pd.DataFrame()
    diagO2, diagO2camp = optimizar_o2_campana(lances, opts, carpeta, o2_prev)
    if not diagO2.empty:
        diagO2.to_csv(o2_acum_path, index=False, encoding="utf-8-sig")

    filas = []
    for l in lances:
        p = l.plan

        # Estado CTD: geometría, Loop Edit y fallback. Separado de O2 para que
        # una discrepancia biogeoquímica no convierta automáticamente todo el CTD en dudoso.
        razones_ctd = []
        if p["confianza_cast"] < AUTO_CFG["confianza_min_cast"]:
            razones_ctd.append("cast ambiguo")
        if p["retencion_loop"] < AUTO_CFG["retencion_min_loop"]:
            razones_ctd.append("Loop Edit retiene < umbral operativo")
        if p.get("deteccion_fallback"):
            razones_ctd.append("inicio por fallback")
        estado_ctd = "REVISAR" if razones_ctd else "AUTO_OK"

        # Estado O2: solo aplica si existe SBE43. Los umbrales son de cola de
        # revisión, no límites universales ni criterios de certificación.
        sbe43 = [s for s in l.config.sensores if s.tipo == "sbe43"]
        razones_o2 = []
        if not sbe43:
            estado_o2 = "NO_APLICA"
            borde_o2 = False
        else:
            borde_o2 = bool(p.get("o2_minimo_borde_campana", p.get("o2_minimo_borde", False)))
            if p.get("o2_confianza", 0.0) < AUTO_CFG["confianza_min_o2"]:
                razones_o2.append("O2 dinámico poco definido")
            if borde_o2:
                razones_o2.append("mínimo O2 en borde del rango")
            penalty = p.get("o2_penalizacion_optimo_rel", np.nan)
            if np.isfinite(penalty) and penalty > float(AUTO_CFG["O2_penalizacion_optimo_revisar_rel"]):
                razones_o2.append("O2 consenso penaliza >15% frente al óptimo del lance")
            gain0 = p.get("o2_mejora_vs_cero_rel", np.nan)
            if np.isfinite(gain0) and gain0 < -float(AUTO_CFG["O2_empeora_vs_cero_revisar_rel"]):
                razones_o2.append("O2 consenso empeora >3% frente a lag=0")
            estado_o2 = "REVISAR" if razones_o2 else "AUTO_OK"

        estado_operativo = "REVISAR" if (estado_ctd == "REVISAR" or estado_o2 == "REVISAR") else "AUTO_OK"
        p["estado_CTD"] = estado_ctd
        p["estado_O2"] = estado_o2
        p["estado_operativo"] = estado_operativo
        p["estado_auto"] = estado_operativo  # compatibilidad con figuras/lógica previa
        p["motivo_revision_CTD"] = "; ".join(razones_ctd)
        p["motivo_revision_O2"] = "; ".join(razones_o2)
        p["motivo_revision"] = "; ".join(razones_ctd + razones_o2)
        guardar_diag_auto(l, carpeta)

        filas.append({
            "estacion": l.meta["estacion"],
            "hex": l.ruta.name,
            "estado_CTD": estado_ctd,
            "estado_O2": estado_o2,
            "estado_operativo": estado_operativo,
            "estado": estado_operativo,
            "confianza_cast": p["confianza_cast"],
            "scan_inicio": p["scan_inicio"],
            "scan_fondo": p["scan_fondo"],
            "scan_fin": p["scan_fin"],
            "pmax_dbar": p["pmax_dbar"],
            "offset_P_dbar": p["offset_presion_dbar"],
            "vmed_down_m_s": p["velocidad_mediana_down_m_s"],
            "vmin_loop_m_s": p["velocidad_min_m_s"],
            "retencion_loop": p["retencion_loop"],
            "align_C_s": opts["adelanto_C_s"],
            "C_metodo": opts.get("C_metodo", ""),
            "C_confianza": opts.get("C_confianza", np.nan),
            "O2_tau": p.get("tau_oxigeno", False),
            "O2_histeresis": p.get("histeresis", False),
            "O2_advance_effective_s": next(iter(p["adelantos_aux_s"].values()), np.nan),
            "O2_lag_s": next(iter(p["adelantos_aux_s"].values()), np.nan),
            "O2_lag_individual_s": p.get("o2_lag_individual_s", np.nan),
            "O2_lag_campana_s": p.get("o2_lag_campana_s", np.nan),
            "O2_score": p.get("o2_score", np.nan),
            "O2_score_lag0": p.get("o2_score_cero_lance", np.nan),
            "O2_score_optimo_tau": p.get("o2_score_optimo_lance_tau", np.nan),
            "O2_mejora_vs_cero_rel": p.get("o2_mejora_vs_cero_rel", np.nan),
            "O2_penalizacion_optimo_rel": p.get("o2_penalizacion_optimo_rel", np.nan),
            "O2_confianza": p.get("o2_confianza", np.nan),
            "O2_minimo_borde": borde_o2,
            "O2_ratio_campana_vs_optimo": p.get("o2_ratio_campana_vs_optimo_lance", np.nan),
            "motivo_revision_CTD": p["motivo_revision_CTD"],
            "motivo_revision_O2": p["motivo_revision_O2"],
            "motivo_revision": p["motivo_revision"]
        })

    return pd.DataFrame(filas), diagC, diagO2, diagO2camp




def _leer_tabla_generica(path):
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xlsm", ".xls"):
        return pd.read_excel(p, sheet_name=0)
    return pd.read_csv(p, sep=None, engine="python")


def _alias_referencias(df):
    x = df.copy()
    x.columns = [norm(c) for c in x.columns]
    aliases = {
        "station": "estacion", "est": "estacion",
        "pressure": "presion_dbar", "pressure_dbar": "presion_dbar", "presion": "presion_dbar",
        "salinity": "salinidad_referencia", "practical_salinity": "salinidad_referencia",
        "sp_ref": "salinidad_referencia", "salinidad_ref": "salinidad_referencia",
        "oxygen_umol_kg": "oxigeno_referencia_umol_kg", "o2_umol_kg": "oxigeno_referencia_umol_kg",
        "oxygen_ml_l": "oxigeno_referencia_ml_l", "o2_ml_l": "oxigeno_referencia_ml_l",
    }
    x.rename(columns={k:v for k,v in aliases.items() if k in x and v not in x}, inplace=True)
    return x


def validar_referencias_externas(base, path, out):
    """Compara CTD con referencias discretas sin modificar coeficientes."""
    if not path:
        return pd.DataFrame(), pd.DataFrame(), "NO_CONFIGURADA"
    pth = Path(path).expanduser()
    if not pth.exists():
        return pd.DataFrame(), pd.DataFrame(), "RUTA_NO_EXISTE"
    try:
        ref = _alias_referencias(_leer_tabla_generica(pth))
    except Exception as exc:
        return pd.DataFrame(), pd.DataFrame(), f"ERROR_LECTURA: {exc}"
    req = {"estacion", "presion_dbar"}
    if not req.issubset(ref.columns):
        return pd.DataFrame(), pd.DataFrame(), "FALTAN_COLUMNAS_ESTACION_PRESION"
    if not any(c in ref for c in ("salinidad_referencia", "oxigeno_referencia_umol_kg", "oxigeno_referencia_ml_l")):
        return pd.DataFrame(), pd.DataFrame(), "SIN_VARIABLES_REFERENCIA"

    tol = float(AUTO_CFG.get("referencia_tolerancia_dbar", 2.0))
    b = base.copy()
    if "direccion" in b:
        b = b[b.direccion.astype(str) == "bajada"].copy()
    b["__est"] = b["estacion"].map(norm) if "estacion" in b else ""
    rows = []
    o2u_cols = [c for c in b if c.startswith("oxigeno_") and c.endswith("_umol_kg")]
    o2m_cols = [c for c in b if c.startswith("oxigeno_") and c.endswith("_ml_L")]

    for _, r in ref.iterrows():
        est = norm(r.get("estacion", ""))
        pr = pd.to_numeric(pd.Series([r.get("presion_dbar")]), errors="coerce").iloc[0]
        if not est or not np.isfinite(pr):
            continue
        cand = b[(b.__est == est) & np.isfinite(pd.to_numeric(b.presion_dbar, errors="coerce"))].copy()
        if cand.empty:
            rows.append({"estacion": r.get("estacion", ""), "presion_ref_dbar": pr, "match": False, "motivo": "estación no encontrada"})
            continue
        dist = np.abs(cand.presion_dbar.to_numpy(float) - float(pr))
        j = int(np.argmin(dist))
        rr = cand.iloc[j]
        ok = float(dist[j]) <= tol
        outrow = {
            "estacion": r.get("estacion", ""), "presion_ref_dbar": float(pr),
            "presion_ctd_dbar": float(rr.presion_dbar), "delta_presion_dbar": float(dist[j]),
            "match": bool(ok), "motivo": "" if ok else "fuera de tolerancia"
        }
        if ok and "salinidad_referencia" in ref and pd.notna(r.get("salinidad_referencia")):
            spref = float(r["salinidad_referencia"])
            spctd = float(rr.get("salinidad_practica", np.nan))
            outrow.update({"SP_ref": spref, "SP_CTD": spctd, "residuo_SP": spctd-spref})
            try:
                lat = float(rr.get("latitud", np.nan)); lon = float(rr.get("longitud", np.nan))
                saref = float(gsw.SA_from_SP(spref, float(pr), lon, lat))
                sactd = float(rr.get("salinidad_absoluta_g_kg", np.nan))
                outrow.update({"SA_ref_g_kg": saref, "SA_CTD_g_kg": sactd, "residuo_SA_g_kg": sactd-saref})
            except Exception:
                pass
        if ok and "oxigeno_referencia_umol_kg" in ref and pd.notna(r.get("oxigeno_referencia_umol_kg")) and o2u_cols:
            vref = float(r["oxigeno_referencia_umol_kg"]); vctd = float(rr.get(o2u_cols[0], np.nan))
            outrow.update({"O2_ref_umol_kg": vref, "O2_CTD_umol_kg": vctd,
                           "residuo_O2_umol_kg": vctd-vref,
                           "residuo_O2_pct": 100*(vctd-vref)/vref if vref != 0 else np.nan})
        if ok and "oxigeno_referencia_ml_l" in ref and pd.notna(r.get("oxigeno_referencia_ml_l")) and o2m_cols:
            vref = float(r["oxigeno_referencia_ml_l"]); vctd = float(rr.get(o2m_cols[0], np.nan))
            outrow.update({"O2_ref_ml_L": vref, "O2_CTD_ml_L": vctd,
                           "residuo_O2_ml_L": vctd-vref,
                           "residuo_O2_pct": 100*(vctd-vref)/vref if vref != 0 else np.nan})
        rows.append(outrow)

    detalle = pd.DataFrame(rows)
    summaries = []
    for c in ["residuo_SP", "residuo_SA_g_kg", "residuo_O2_umol_kg", "residuo_O2_ml_L", "residuo_O2_pct"]:
        if c not in detalle:
            continue
        mask = detalle["match"].astype(bool) if "match" in detalle else np.zeros(len(detalle), bool)
        x = pd.to_numeric(detalle.loc[mask, c], errors="coerce").dropna().to_numpy(float)
        if len(x):
            summaries.append({
                "variable_residuo": c, "n": len(x), "bias": float(np.mean(x)),
                "rmse": float(np.sqrt(np.mean(x*x))), "mad": float(np.median(np.abs(x-np.median(x)))),
                "max_abs": float(np.max(np.abs(x)))
            })
    resumen = pd.DataFrame(summaries)
    detalle.to_csv(out/"validacion_referencias.csv", index=False, encoding="utf-8-sig")
    resumen.to_csv(out/"resumen_validacion_referencias.csv", index=False, encoding="utf-8-sig")
    return detalle, resumen, "OK"


def _fecha_campana(base, decisiones):
    fechas = []
    if "fecha_bitacora" in base:
        # leer_tabla normaliza esta columna a YYYY-MM-DD; interpretar ISO primero
        # evita la ambigüedad dayfirst (p. ej., 2026-09-01).
        f = pd.to_datetime(base["fecha_bitacora"], format="%Y-%m-%d", errors="coerce")
        if f.notna().sum() == 0:
            f = pd.to_datetime(base["fecha_bitacora"], errors="coerce", dayfirst=True)
        fechas.extend(list(f.dropna()))
    if not fechas and decisiones is not None and "hex" in decisiones:
        for name in decisiones.hex.astype(str):
            m = re.search(r"(\d{2})[-_](\d{2})[-_](\d{4})", name)
            if m:
                try:
                    fechas.append(pd.Timestamp(year=int(m.group(3)), month=int(m.group(2)), day=int(m.group(1))))
                except Exception:
                    pass
    return max(fechas) if fechas else pd.NaT


def generar_auditoria_estandares(base, resumen, decisiones, cfg, validacion, resumen_validacion, ref_status, out):
    """Genera auditoría explícita; nunca declara certificación GO-SHIP automáticamente."""
    rows = []
    def add(marco, criterio, estado, evidencia, accion=""):
        rows.append({"marco": marco, "criterio": criterio, "estado": estado, "evidencia": evidencia, "accion": accion})

    add("Sea-Bird", "Integridad HEX/XMLCON y seriales", "CUMPLE",
        "El decodificador valida ancho de trama, HEX completo, configuración, seriales T/C y ScansToAverage antes de procesar.")
    add("Sea-Bird", "Secuencia de procesamiento", "CUMPLE_CON_NOTA",
        "Conversión → histéresis SBE43 (si aplica) → filtros → alineamiento → Cell Thermal Mass → Loop Edit/máscara → Derive final → Bin Average. Se usa una derivación provisional continua solo para correcciones dinámicas de O2.")
    add("Producto", "Downcast y binning", "CUMPLE",
        f"Producto principal de bajada; bins centrados cada {AUTO_CFG['bin_dbar']} dbar; n mínimo por bin={AUTO_CFG['nmin_bin']}.")
    add("TEOS-10", "SP, SA, CT, sigma0, densidad y profundidad", "CUMPLE",
        "Derivados con GSW-Python; se añade control gsw.infunnel para la aproximación de 75 términos.")
    add("QARTOD/IOOS", "Flags primarios interoperables", "CUMPLE",
        "Se generan flags 1/2/3/4/9 separados del bitmask interno; gross range, spikes diagnósticos, freezing, densidad/inversiones y estado O2.")
    add("QARTOD/IOOS", "Climatología regional", "NO_EVALUADO",
        "No se configuró una climatología externa. El código no inventa límites climatológicos.",
        "Opcional: incorporar WOA/atlas regional o climatología propia documentada.")

    nc_ok = False
    if resumen is not None and "netcdf_cf_estado" in resumen:
        vals = resumen.netcdf_cf_estado.astype(str)
        nc_ok = len(vals) > 0 and bool((vals == "OK").all())
    add("CF-1.13", "NetCDF y metadatos/flags", "GENERADO" if nc_ok else "PARCIAL",
        "NetCDF CF con standard_name, unidades, ancillary_variables y flag_values/flag_meanings." if nc_ok else
        "Exportación CF no se completó en todos los perfiles; CSV/Excel siguen disponibles.",
        "Antes de depósito público, ejecutar cf-checker/compliance-checker." if nc_ok else "Instalar xarray y revisar log de exportación.")

    has_sal = (validacion is not None and not validacion.empty and
               any(c in validacion for c in ("residuo_SP", "residuo_SA_g_kg")))
    has_o2 = (validacion is not None and not validacion.empty and
              any(c in validacion for c in ("residuo_O2_umol_kg", "residuo_O2_ml_L")))
    add("GO-SHIP/alta exactitud", "Validación salinidad con muestras discretas", "EVALUADO" if has_sal else "PENDIENTE",
        "Existe comparación CTD-salinómetro." if has_sal else "No se suministró salinidad discreta independiente.",
        "Usar salinómetro trazable/IAPSO si se requiere demostrar exactitud absoluta." if not has_sal else "")
    sbe43_present = any(ss.tipo == "sbe43" for ss in cfg.sensores)
    add("GO-SHIP/alta exactitud", "Validación O2 con Winkler", "EVALUADO" if has_o2 else ("PENDIENTE" if sbe43_present else "NO_APLICA"),
        "Existe comparación SBE43-Winkler." if has_o2 else ("SBE43 presente sin referencia Winkler cargada." if sbe43_present else "Sin SBE43."),
        "Tomar/analizar muestras Winkler para calibración/verificación absoluta del O2." if sbe43_present and not has_o2 else "")

    fecha = _fecha_campana(base, decisiones)
    cal_evidence = []
    stale = False
    if pd.notna(fecha):
        for ss in cfg.sensores:
            if not ss.fecha_calibracion:
                continue
            fc = pd.to_datetime(ss.fecha_calibracion, errors="coerce", dayfirst=True)
            if pd.isna(fc):
                continue
            months = (fecha.year-fc.year)*12 + (fecha.month-fc.month) + (fecha.day-fc.day)/30.44
            cal_evidence.append(f"{ss.etiqueta} SN {ss.serie}: {months:.1f} meses")
            stale = stale or months > 18
    add("Fabricante", "Trazabilidad de calibración", "ADVERTENCIA" if stale else ("DOCUMENTADA" if cal_evidence else "NO_EVALUADO"),
        "; ".join(cal_evidence) if cal_evidence else "Fecha de campaña/calibración insuficiente para calcular antigüedad.",
        "Sea-Bird recomienda típicamente recalibración anual; con intervalos largos, reforzar validación con muestras y/o calibración post-crucero." if stale else "")

    add("GO-SHIP", "Calidad climática/reference quality", "NO_DEMOSTRADO",
        "El procesamiento sigue buenas prácticas, pero software+XMLCON no demuestran por sí solos 0.002 °C, 0.002 g/kg, 3 dbar y 1% O2. El SBE19plus V2 especifica T inicial ±0.005 °C y el SBE43 ±2% de saturación.",
        "No etiquetar el conjunto como GO-SHIP/climate-quality sin calibración trazable y validación discreta suficiente.")

    op_bad = int((decisiones.estado_operativo == "REVISAR").sum()) if decisiones is not None and "estado_operativo" in decisiones else np.nan
    add("Publicación", "QC reproducible y trazabilidad",
        "PROVISIONAL" if (not has_sal or (sbe43_present and not has_o2)) else "APTO_PARA_REVISION_CIENTIFICA_FINAL",
        f"Algoritmo {VERSION}; {op_bad} perfiles en cola operativa de revisión. Referencias externas: {ref_status}.",
        "Resolver perfiles REVISAR y documentar validaciones externas antes de afirmar exactitud absoluta.")

    audit = pd.DataFrame(rows)
    audit.to_csv(out/"AUDITORIA_ESTANDARES_CTD.csv", index=False, encoding="utf-8-sig")
    return audit


def generar_metodos_publicacion(base, decisiones, cfg, audit, resumen_validacion, out):
    sensors = []
    for ss in cfg.sensores:
        if ss.tipo != "sin_uso":
            sensors.append(f"- {ss.etiqueta}: serie {ss.serie or 'NA'}, calibración {ss.fecha_calibracion or 'NA'}, tipo {ss.tipo}")
    c_lag = np.nan
    o2_lag = np.nan
    o2_tau = None
    if decisiones is not None and not decisiones.empty:
        if "align_C_s" in decisiones:
            c_lag = float(pd.to_numeric(decisiones.align_C_s, errors="coerce").median())
        if "O2_advance_effective_s" in decisiones:
            o2_lag = float(pd.to_numeric(decisiones.O2_advance_effective_s, errors="coerce").median())
        if "O2_tau" in decisiones and not decisiones.O2_tau.mode().empty:
            o2_tau = bool(decisiones.O2_tau.mode().iloc[0])
    valtxt = "No se cargaron referencias discretas independientes."
    if resumen_validacion is not None and not resumen_validacion.empty:
        valtxt = resumen_validacion.to_string(index=False)
    lines = [
        "# Métodos CTD — texto de trazabilidad para publicación", "",
        f"Generado automáticamente por **{VERSION}**. Este documento describe el procesamiento; no constituye certificación GO-SHIP.", "",
        "## Instrumentación", "", *sensors, "",
        "## Procesamiento", "",
        f"Los archivos HEX originales se decodificaron usando el XMLCON correspondiente. Se verificaron integridad de trama, configuración y seriales. La secuencia aplicada fue: conversión a unidades físicas; corrección de histéresis SBE43 cuando correspondió; filtrado; alineamiento temporal; corrección de masa térmica de la celda de conductividad; detección automática de atemperamiento, descenso, fondo y subida; Loop Edit adaptativo; cálculo final de variables derivadas; selección del downcast y promedio en bins de {AUTO_CFG['bin_dbar']} dbar.", "",
        "La Salinidad Práctica (SP) se calculó en PSS-78 y las propiedades termodinámicas se calcularon con TEOS-10/GSW, incluyendo Salinidad Absoluta (SA), Temperatura Conservativa (CT), sigma0, densidad y profundidad. La profundidad positiva se conserva internamente; las figuras representan z negativo bajo la superficie.", "",
        f"El alineamiento de conductividad se estimó a escala de campaña minimizando artefactos de salinidad; valor representativo: **{c_lag:.3f} s**. Para el SBE43, el adelanto efectivo se estimó por una función objetivo de coherencia descenso/subida y se evaluó Tau de forma explícita; valor representativo: **{o2_lag:.2f} s**, Tau={o2_tau}. La histéresis se aplica a perfiles profundos cuando los coeficientes H1/H2/H3 están disponibles. Los coeficientes absolutos de calibración no se ajustan a partir del propio CTD.", "",
        "## Control de calidad", "",
        "Se conserva un bitmask interno de procedencia y se añaden flags primarios QARTOD/IOOS (1=pass, 2=not evaluated, 3=suspect/high interest, 4=fail, 9=missing) para presión, temperatura, salinidad y oxígeno. Se ejecutan controles de rango, spikes diagnósticos, temperatura de congelación TEOS-10, inversiones de densidad/N² y dominio gsw.infunnel. Los tests diagnósticos no eliminan automáticamente estructuras oceanográficas reales: los valores sospechosos se conservan y se marcan para revisión.", "",
        "## Validación externa", "", valtxt, "",
        "## Limitaciones", "",
        "La exactitud absoluta no puede inferirse exclusivamente de la suavidad del perfil o del ajuste descenso/subida. Para demostrar calidad de referencia comparable a GO-SHIP se requieren calibraciones trazables y referencias independientes, en particular salinidad de botella/salinómetro y O2 Winkler. El producto debe someterse a revisión científica final antes de publicación o depósito.", ""
    ]
    path = out/"METODOS_CTD_PUBLICACION.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def exportar_excel_crucero(base, resumen, manifiesto, decisiones, invent, out, audit=None, validacion=None, resumen_validacion=None):
    xlsx = out / "CTD_BASE_FINAL.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
        base.to_excel(writer, sheet_name="Base_CTD", index=False)
        resumen.to_excel(writer, sheet_name="Resumen_QC", index=False)
        decisiones.to_excel(writer, sheet_name="Decisiones_auto", index=False)
        invent.to_excel(writer, sheet_name="Sensores_XMLCON", index=False)
        if audit is not None and not audit.empty:
            audit.to_excel(writer, sheet_name="Auditoria_estandares", index=False)
        if validacion is not None and not validacion.empty:
            validacion.to_excel(writer, sheet_name="Validacion_referencias", index=False)
        if resumen_validacion is not None and not resumen_validacion.empty:
            resumen_validacion.to_excel(writer, sheet_name="Resumen_validacion", index=False)

        proc = []
        for k, v in manifiesto.items():
            if k == "lances":
                continue
            if isinstance(v, (dict, list)):
                v = json.dumps(v, ensure_ascii=False, default=json_default)
            proc.append({"clave": k, "valor": str(v)})
        pd.DataFrame(proc).to_excel(writer, sheet_name="Procedencia", index=False)

        wb = writer.book
        for ws in wb.worksheets:
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
            for col in ws.columns:
                letra = col[0].column_letter
                mx = max(
                    (len(str(c.value)) if c.value is not None else 0 for c in col[:200]),
                    default=8
                )
                ws.column_dimensions[letra].width = min(45, max(10, mx + 2))
    return xlsx


def _ruta_auto(cfgkey, texto, carpeta=False, extensiones=()):
    val = AUTO_CFG.get(cfgkey)
    if val:
        p = Path(val).expanduser().resolve()
        if (carpeta and p.is_dir()) or ((not carpeta) and p.is_file()):
            return p
    return ruta_input(texto, carpeta=carpeta, extensiones=extensiones)


def _cargar_json_seguro(path, defecto):
    try:
        if Path(path).exists():
            return json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        pass
    return defecto


def _merge_por_clave(old, new, clave):
    if old is None or len(old) == 0:
        return new.copy()
    if new is None or len(new) == 0:
        return old.copy()
    old = old.copy()
    new = new.copy()
    if clave in old.columns and clave in new.columns:
        keys = set(new[clave].astype(str))
        old = old[~old[clave].astype(str).isin(keys)]
    return pd.concat([old, new], ignore_index=True, sort=False)


def _preparar_incremental(out, rutas, df):
    """Devuelve índices nuevos/cambiados; puede forzar reproceso tras fijar parámetros de campaña."""
    state_path = out / "estado_procesamiento.json"
    state = _cargar_json_seguro(state_path, {"archivos": {}, "parametros_campana": {}})
    indices = []
    hashes = {}
    force_all = bool(state.get("forzar_reproceso_completo_proxima_ejecucion", False))
    version_changed = bool(state.get("archivos")) and state.get("software_version") != VERSION
    if version_changed:
        force_all = True
        print(f"\nReproceso completo automático: cambió la versión de procesamiento a {VERSION}.")
    for i, ruta in enumerate(rutas):
        rel = str(df.iloc[i].get("archivo_hex", ruta.name))
        h = sha256(ruta)
        hashes[rel] = h
        if force_all or state.get("archivos", {}).get(rel) != h:
            indices.append(i)
    if force_all and not version_changed:
        print("\nReproceso completo automático: los parámetros dinámicos de campaña se estabilizaron en la corrida anterior.")
    return state, state_path, indices, hashes


def _archivar_manifiesto_previo(out):
    man = out / "manifiesto.json"
    if not man.exists():
        return
    hist = out / "historial_manifiestos"
    hist.mkdir(exist_ok=True)
    dst = hist / ("manifiesto_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f") + ".json")
    shutil.copy2(man, dst)


def _cargar_base_previa(out):
    bp = out / "base_ctd_procesada.csv"
    rp = out / "resumen_procesamiento.csv"
    dp = out / "DECISIONES_AUTOMATICAS.csv"
    old_base = pd.read_csv(bp) if bp.exists() else pd.DataFrame()
    old_res = pd.read_csv(rp) if rp.exists() else pd.DataFrame()
    old_dec = pd.read_csv(dp) if dp.exists() else pd.DataFrame()
    return old_base, old_res, old_dec

# -----------------------------------------------------------------------------
# v2.5: reconocimiento de someros, conciliación y primera pasada provisional
# -----------------------------------------------------------------------------

def _v23_datetime(v):
    if v is None or str(v).strip() in ("", "NaT", "nan", "None"):
        return pd.NaT
    try:
        ts = pd.Timestamp(v)
        return ts.tz_localize(None) if ts.tzinfo is not None else ts
    except (TypeError, ValueError, OverflowError):
        return pd.NaT


def _v23_utc(v, offset):
    ts = _v23_datetime(v)
    if pd.isna(ts) or offset is None:
        return pd.NaT
    return ts - pd.Timedelta(hours=float(offset))


def _v23_tipo(datos, dt):
    """Clasificación conservadora basada en inmersión y movimiento coherente.

    No exige profundidades mayores de 2-3 m ni ascenso final. Un período de
    atemperamiento/calibración operativa en superficie NO equivale a una
    calibración metrológica del sensor. No se descarta ningún HEX original.
    """
    p = pd.to_numeric(datos['presion_dbar'], errors='coerce').to_numpy(float)
    c = pd.to_numeric(datos['conductividad_S_m'], errors='coerce').to_numpy(float)
    ok = np.isfinite(p) & np.isfinite(c)
    if ok.sum() < 20 or not np.isfinite(dt) or dt <= 0:
        return 'REVISAR', 'Muestras instrumentales insuficientes'
    pp = p[ok]
    cc = c[ok]
    dur = len(pp) * float(dt)
    # Umbrales de C para agua de mar. Agua dulce/salobre muy diluida exige revisión.
    frac_agua = float(np.mean(cc > 0.5))
    frac_aire = float(np.mean(cc < 0.10))
    # Evidencia de prueba en aire requiere AMBAS señales; evita confundir
    # estaciones marinas someras (0-2 dbar) con prueba fuera del agua.
    if frac_aire >= 0.90 and np.nanpercentile(np.abs(pp), 90) < 1.0:
        return 'PRUEBA_AIRE_CANDIDATA', 'Presión casi atmosférica y conductividad compatible con aire'
    if frac_agua < 0.50:
        return 'REVISAR', 'Conductividad sin evidencia suficiente de inmersión marina'
    # Suavizado robusto: escala >= 1.5 s y a lo sumo 31 observaciones.
    window = min(31, max(3, int(round(2.0 / dt)) | 1))
    smooth = pd.Series(pp).rolling(window, center=True, min_periods=1).median().to_numpy()
    q05, q95 = np.nanpercentile(smooth, [5, 95])
    excursion = float(q95 - q05)
    # Ruido residual: diferencia respecto a curva medianamente suavizada.
    noise = 1.4826 * float(np.nanmedian(np.abs(pp - smooth)))
    threshold = max(0.75, 8.0 * noise)
    # Trayectoria de descenso que dura un tramo temporal (no un único pico).
    stride = max(1, int(round(2.0 / dt)))
    coarse = smooth[::stride]
    if len(coarse) < 4:
        coarse = smooth
    progress = float(np.nanmax(coarse - np.minimum.accumulate(coarse)))
    total_var = float(np.nansum(np.abs(np.diff(coarse))))
    directionality = progress / max(total_var, 1e-9)
    # Un verdadero perfil de 2-3 m suele mostrar una rampa coherente >=0.75 dbar.
    # Si las oscilaciones dominan (oleaje o fondeo), no declarar perfil automáticamente.
    if excursion >= threshold and progress >= threshold and directionality >= 0.20:
        return 'PERFIL_CANDIDATO', (f'Descenso coherente, excursión={excursion:.2f} dbar, '
                                   f'progreso={progress:.2f} dbar; compatible con perfil somero')
    # Estático: agua comprobada, duración suficiente, excursión pequeña O
    # oscilaciones sin descenso coherente. Variación fuerte queda en REVISAR.
    if dur >= 60 and excursion < max(1.0, 10.0 * noise):
        return 'ESTATICO_CANDIDATO', ('Registro sumergido estable sin descenso identificable; '
                                    f'excursión={excursion:.2f} dbar')
    return 'REVISAR', ('Variación ambigua: oleaje, operación somera, perfil parcial '
                      'o múltiples eventos; requiere gráfico de presión-tiempo')


def _v26_procesar_estatico(datos, cfg, inicio, archivo, salida, meta=None):
    """Serie temporal estática: sin downcast/Loop Edit/masa térmica ni corrección tau.

    Recorta tramo sumergido y estable; conserva tiempos originales y datos crudos
    convertidos, con indicadores QC, sin sustituir una calibración certificada.
    No inventa lat/lon ni calcula SA/CT/densidad TEOS-10 sin ellas.
    """
    meta = dict(meta or {})
    df = datos.copy().reset_index(drop=True)
    df['archivo_hex'] = archivo
    df['estacion'] = str(meta.get('estacion','SIN_ESTACION'))
    lat = pd.to_numeric(meta.get('latitud', np.nan), errors='coerce')
    lon = pd.to_numeric(meta.get('longitud', np.nan), errors='coerce')
    if np.isfinite(lat) and np.isfinite(lon) and abs(lat)<=90 and abs(lon)<=180:
        df['latitud'] = float(lat)
        df['longitud'] = float(lon)
    dt = float(cfg.navg / 4)
    if dt <= 0 or len(df) < 20:
        raise ValueError('Serie insuficiente para estático')
    p = pd.to_numeric(df['presion_dbar'], errors='coerce').to_numpy(float)
    c = pd.to_numeric(df['conductividad_S_m'], errors='coerce').to_numpy(float)
    t = pd.to_numeric(df['temperatura_C'], errors='coerce').to_numpy(float)
    # Umbral C marino: para aguas dulces, someras o salobres atípicas revisar manualmente.
    wet = np.isfinite(p) & np.isfinite(c) & np.isfinite(t) & (c > 0.5)
    if wet.sum() < max(20, round(30/dt)):
        raise ValueError('No existe suficiente evidencia de inmersión marina estable')
    pref = float(np.nanmedian(p[wet]))
    # 0.5 dbar admite oleaje pequeño; el MAD evita rigidez ante ruido instrumental.
    mad = float(np.nanmedian(np.abs(p[wet]-pref)))
    banda = max(0.5, 4 * 1.4826 * mad)
    estable = wet & (np.abs(p-pref) <= banda)
    # Evita mezclar tiempos discontinuos: usar el segmento estable continuo más largo.
    tramos = bloques(estable)
    a,b = max(tramos, key=lambda ab: ab[1]-ab[0])
    if (b-a)*dt < 30:
        raise ValueError('No se detectó intervalo estable continuo de al menos 30 s')
    df['qc_estatico'] = np.where(estable, 1, 3).astype('uint8')  # 1 pass, 3 suspect
    df.loc[:a-1,'qc_estatico'] = 4
    df.loc[b:,'qc_estatico'] = 4
    df['intervalo_estable'] = False
    df.loc[a:b-1,'intervalo_estable'] = True
    with np.errstate(invalid='ignore',divide='ignore'):
        df['salinidad_practica'] = gsw.SP_from_C(c*10,t,p)
    if np.isfinite(lat) and np.isfinite(lon) and abs(lat)<=90 and abs(lon)<=180:
        SA = gsw.SA_from_SP(df['salinidad_practica'].to_numpy(float),p,float(lon),float(lat))
        df['salinidad_absoluta_g_kg'] = SA
        df['temperatura_conservativa_C'] = gsw.CT_from_t(SA,t,p)
        df['sigma0_kg_m3'] = gsw.sigma0(SA,df['temperatura_conservativa_C'].to_numpy(float))
    if pd.notna(_v23_datetime(inicio)):
        df['fecha_hora_ctd'] = (_v23_datetime(inicio) + pd.to_timedelta(df['tiempo_s'],unit='s'))
    # Convertir SBE43 sin tau/histéresis dinámica: calibración del XMLCON sí se aplica.
    # El oxígeno SBE63 se transforma con sus coeficientes instrumentales.
    advertencias=[]
    for sensor in cfg.sensores[3:]:
        tag=f's{sensor.indice}'
        try:
            if sensor.tipo=='sbe43':
                ml=oxigeno43(df,sensor,dt,False,1)
            elif sensor.tipo=='sbe63':
                co=sbcoef.Oxygen63Coefficients(**{k:valor(sensor,k) for k in sbcoef.Oxygen63Coefficients.__annotations__})
                tc=sbcoef.Thermistor63Coefficients(**{k:valor(sensor,k) for k in sbcoef.Thermistor63Coefficients.__annotations__})
                ml=sbconv.convert_sbe63_oxygen(df[sensor.senal].to_numpy(),df.sbe63_temp_V.to_numpy(),p,df.salinidad_practica.to_numpy(),co,tc)
            else:
                continue
            df[f'oxigeno_{tag}_ml_L']=ml
            df[f'oxigeno_{tag}_mg_L']=sbconv.convert_oxygen_to_mg_per_l(ml)
        except Exception as exc:
            advertencias.append(f'{tag}: oxígeno no convertido ({exc})')
    usable = df.iloc[a:b].copy()
    columnas = ['presion_dbar','temperatura_C','conductividad_S_m','salinidad_practica']
    resumen = {'archivo_hex':archivo,'estacion':str(meta.get('estacion','SIN_ESTACION')),
               'latitud':lat,'longitud':lon,'fecha_bitacora':meta.get('fecha_bitacora',''),
               'hora_inmersion':meta.get('hora_inmersion',''),'inicio_hex':str(inicio),
               'inicio_estable_s':float(df.tiempo_s.iloc[a]),
               'fin_estable_s':float(df.tiempo_s.iloc[b-1]),
               'duracion_estable_min':round((b-a)*dt/60,3),
               'n_estable':int(b-a), 'banda_presion_dbar':float(banda),
               'estado':'ESTATICO_PROCESADO_PRELIMINAR',
               'observaciones':'; '.join(advertencias)}
    for col in columnas + [k for k in usable if k.startswith('oxigeno_')]:
        x=pd.to_numeric(usable[col],errors='coerce').to_numpy(float)
        resumen[col+'_media']=float(np.nanmean(x)) if np.isfinite(x).any() else None
        resumen[col+'_sd']=float(np.nanstd(x,ddof=1)) if np.isfinite(x).sum()>1 else None
    salida.mkdir(parents=True,exist_ok=True)
    sello=_v23_datetime(inicio)
    suf=sello.strftime('%Y%m%d_%H%M%S') if pd.notna(sello) else 'SIN_FECHA'
    nombre='ESTATICO_'+slug(meta.get('estacion','SIN_ESTACION'))+'_'+suf+'_'+hashlib.sha256(archivo.encode()).hexdigest()[:8]
    xlsx=salida/(nombre+'.xlsx')
    with pd.ExcelWriter(xlsx,engine='openpyxl') as writer:
        usable.to_excel(writer,sheet_name='serie_estable',index=False)
        df.to_excel(writer,sheet_name='serie_completa_QC',index=False)
        pd.DataFrame([resumen]).to_excel(writer,sheet_name='resumen',index=False)
    # PNG independientes y legibles: una variable por eje, tiempo real en minutos.
    cols=[c for c in columnas+['salinidad_absoluta_g_kg','temperatura_conservativa_C','sigma0_kg_m3'] if c in df]+[c for c in df if c.endswith('_ml_L')]
    fig,ax=plt.subplots(len(cols),1,figsize=(11,max(7,len(cols)*2.2)),sharex=True,squeeze=False)
    for axis,col in zip(ax[:,0],cols):
        axis.plot(df.tiempo_s/60,df[col],lw=.75,color='0.7',label='Registro completo')
        axis.plot(usable.tiempo_s/60,usable[col],lw=1.1,label='Intervalo estable')
        axis.set_ylabel(col.replace('_',' '),fontsize=9)
        axis.grid(alpha=.25)
    ax[0,0].legend(loc='best',fontsize=8)
    ax[-1,0].set_xlabel('Tiempo desde inicio del registro (min)')
    fig.suptitle(f"{meta.get('estacion','SIN_ESTACION')} | {archivo} | estático | sin correcciones de perfil",fontsize=12)
    fig.tight_layout(rect=(0,0,1,.97))
    fig.savefig(salida/(nombre+'.png'),dpi=300,bbox_inches='tight');plt.close(fig)
    return xlsx,resumen


def _v23_inventariar(carpeta, cfg, out, guardar_provisional, estaciones=None):
    registros=[]
    estaciones = estaciones or {}
    paths=sorted((p for p in carpeta.rglob('*') if p.is_file() and p.suffix.lower()=='.hex'),
                 key=lambda p: p.relative_to(carpeta).as_posix().lower())
    if not paths: raise ValueError('No hay archivos .hex en la carpeta.')
    diag=out/'00_DESCUBRIMIENTO';diag.mkdir(parents=True,exist_ok=True)
    for i,p in enumerate(paths,1):
        rel=p.relative_to(carpeta).as_posix()
        row=dict(archivo_hex=rel,sha256=sha256(p),inicio_hex='',fuente_hora='',
                 fin_hex='',duracion_min=np.nan,pmax_dbar=np.nan,delta_p_dbar=np.nan,
                 n_scans=0,tipo='REVISAR',motivo='',estado='INVENTARIADO')
        try:
            # El lector v2.2 verifica formato y serie de sensores contra XMLCON.
            head,raw=leer_hex(p,cfg)
            datos,_=convertir_base(raw,cfg)
            inicio,fuente=fecha_header(head)
            dt=cfg.navg/4
            row.update(inicio_hex=inicio,fuente_hora=fuente,
                       duracion_min=round(len(raw)*dt/60,3),n_scans=len(raw),
                       pmax_dbar=round(float(np.nanmax(datos.presion_dbar)),3),
                       delta_p_dbar=round(float(np.nanpercentile(datos.presion_dbar,95)-
                                                    np.nanpercentile(datos.presion_dbar,5)),3))
            tipo,motivo=_v23_tipo(datos,dt)
            row.update(tipo=tipo,motivo=motivo,estado='DECODIFICADO')
            ts=_v23_datetime(inicio)
            if not pd.isna(ts):row['fin_hex']=(ts+pd.Timedelta(seconds=(len(raw)-1)*dt)).isoformat()
            if tipo == 'ESTATICO_CANDIDATO':
                try:
                    xlsx, resumen_static = _v26_procesar_estatico(
                        datos,cfg,inicio,rel,out/'02_ESTATICOS_PROCESADOS',
                        meta=estaciones.get(rel, estaciones.get(p.name, {})))
                    row['producto_estatico']=str(xlsx)
                    row['qc_estatico']='PROCESADO_CON_ESTACION' if estaciones.get(rel, estaciones.get(p.name, {})) else 'PRELIMINAR'
                except Exception as exc_est:
                    row['qc_estatico']='REVISAR'
                    row['error_estatico']=f'{type(exc_est).__name__}: {exc_est}'
            if guardar_provisional and tipo == 'PERFIL_CANDIDATO':
                # Primera pasada: solo variables independientes de latitud/longitud.
                datos['salinidad_practica']=gsw.SP_from_C(
                    datos.conductividad_S_m.to_numpy(float)*10,
                    datos.temperatura_C.to_numpy(float),datos.presion_dbar.to_numpy(float))
                base=slug(p.stem)+'_'+row['sha256'][:10]
                datos.to_csv(diag/(base+'_calibrado_preliminar.csv'),index=False,encoding='utf-8-sig')
                fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True)
                axes[0].plot(datos.tiempo_s/60,datos.presion_dbar,lw=.8)
                axes[0].invert_yaxis();axes[0].set_ylabel('Presión (dbar)');axes[0].grid(alpha=.3)
                axes[1].plot(datos.tiempo_s/60,datos.temperatura_C,lw=.8)
                axes[1].set_ylabel('Temperatura (°C)');axes[1].set_xlabel('Tiempo (min)')
                axes[1].grid(alpha=.3);fig.suptitle(rel+' | preliminar, no QC final')
                fig.tight_layout();fig.savefig(diag/(base+'_diagnostico.png'),dpi=300)
                plt.close(fig)
        except Exception as exc:
            row.update(estado='ERROR',motivo=f'{type(exc).__name__}: {exc}')
        # Solo las copias de pruebas/estáticos reciben nombres cronológicos.
        if row['tipo'] in ('ESTATICO_CANDIDATO','PRUEBA_AIRE_CANDIDATA'):
            clase='ESTATICOS' if row['tipo']=='ESTATICO_CANDIDATO' else 'PRUEBAS'
            dst=out/'00_DESCUBRIMIENTO'/'NO_PROCESADOS'/clase
            dst.mkdir(parents=True,exist_ok=True)
            ts=_v23_datetime(row['inicio_hex'])
            stamp=ts.strftime('%Y%m%d_%H%M%S') if pd.notna(ts) else 'SIN_FECHA'
            nombre=f'{clase[:-1]}_{stamp}_{row["sha256"][:10]}.hex'
            destino=dst/nombre
            if not destino.exists():shutil.copy2(p,destino)
            row['copia_renombrada']=str(destino)
        registros.append(row)
        print(f'  {i:>3}/{len(paths)} {rel}: {row["tipo"]} | {row["inicio_hex"]} | {row["estado"]}')
    inv=pd.DataFrame(registros).sort_values(['inicio_hex','archivo_hex'],kind='stable')
    inv.to_csv(diag/'inventario_hex.csv',index=False,encoding='utf-8-sig')
    return inv.reset_index(drop=True)


def _v23_plantilla(inv,out):
    df=inv[['archivo_hex','tipo','inicio_hex','fuente_hora','motivo']].copy()
    df.insert(1,'estacion','');df.insert(2,'latitud','');df.insert(3,'longitud','')
    df.insert(4,'fecha','');df.insert(5,'hora','')
    path=out/'00_DESCUBRIMIENTO'/'PLANTILLA_ASIGNACION_ESTACIONES.xlsx'
    with pd.ExcelWriter(path,engine='openpyxl') as writer:
        df.to_excel(writer,sheet_name='asignacion',index=False)
        inv.to_excel(writer,sheet_name='inventario',index=False)
        ws=writer.sheets['asignacion'];ws.freeze_panes='B2';ws.auto_filter.ref=ws.dimensions
        for letter,width in {'A':34,'B':17,'C':14,'D':14,'E':15,'F':15,'G':20,'H':27,'I':35,'J':45}.items():
            ws.column_dimensions[letter].width=width
    return path


def _v23_conciliar(df, inv, out, carpeta):
    """Conciliación global 1:1 basada en fecha/hora, sin corte duro temporal.

    Fecha/hora del HEX es la hora de inicio del registro, NO necesariamente
    la inmersión. Se ordenan candidatos por proximidad al inicio del HEX.
    Propuestas lejanas o con relojes no contrastados SIEMPRE requieren revisión.
    """
    from scipy.optimize import linear_sum_assignment

    df = df.copy().reset_index(drop=True)
    if 'archivo_hex' not in df:
        df['archivo_hex'] = ''
    for c in ('fecha_bitacora', 'hora_inmersion'):
        if c not in df:
            df[c] = ''
    # Evitar tratar NaN como nombre de archivo (p. ej. Excel con celdas vacías).
    df['archivo_hex'] = df['archivo_hex'].fillna('').astype(str).str.strip().replace({'nan': '', 'None': ''})
    df['archivo_hex'] = df['archivo_hex'].str.replace(chr(92), '/', regex=False)
    known = set(inv.archivo_hex.astype(str))
    explicit = [x for x in df.archivo_hex if x]
    if len(explicit) != len(set(explicit)):
        raise ValueError('El Excel utiliza el mismo archivo HEX en varias estaciones.')
    for x in explicit:
        if x not in known:
            raise ValueError(f'Archivo HEX explícito no encontrado: {x}')
    explicit_set = set(explicit)

    # Mantener REVISAR como candidato por si es un perfil somero o atípico;
    # nunca asignar pruebas/estáticos detectados por defecto.
    viable = inv[(inv.estado == 'DECODIFICADO') &
                 (inv.tipo.isin(['PERFIL_CANDIDATO', 'REVISAR'])) &
                 (~inv.archivo_hex.isin(explicit_set))].copy()
    viable = viable.drop_duplicates(subset=['sha256']).reset_index(drop=True)
    paths = viable.archivo_hex.astype(str).tolist()
    start = []
    for _, z in viable.iterrows():
        t = _v23_datetime(z.inicio_hex)
        if pd.notna(t) and AUTO_CFG['hora_ctd_offset_h'] is not None:
            t = t - pd.Timedelta(hours=float(AUTO_CFG['hora_ctd_offset_h']))
        start.append(t)
    dates = []
    for _, r in df.iterrows():
        t = _v23_datetime(f"{r.fecha_bitacora} {r.hora_inmersion}")
        if pd.notna(t) and AUTO_CFG['hora_bitacora_offset_h'] is not None:
            t = t - pd.Timedelta(hours=float(AUTO_CFG['hora_bitacora_offset_h']))
        dates.append(t)
    # Comparar en hora civil si no está definida la zona de ambos relojes.
    # Una zona configurada sin la otra NO permite comparación UTC consistente.
    offset_incompleto = ((AUTO_CFG['hora_ctd_offset_h'] is None) !=
                         (AUTO_CFG['hora_bitacora_offset_h'] is None))
    if offset_incompleto:
        start = [_v23_datetime(v) for v in viable.inicio_hex]
        dates = [_v23_datetime(f"{r.fecha_bitacora} {r.hora_inmersion}") for _, r in df.iterrows()]

    rows = [i for i in range(len(df)) if not df.at[i, 'archivo_hex'] and pd.notna(dates[i])]
    columns = [j for j, t in enumerate(start) if pd.notna(t)]
    n, m = len(rows), len(columns)
    # Emparejar tantos pares como sea posible; después minimizar la suma de
    # diferencias temporales. Columnas ficticias únicamente para exceso de estaciones.
    # El algoritmo NO corta por tolerancia: umbrales son para confianza/revisión.
    cost = np.zeros((n, m+n), dtype=float)
    cost[:, :] = 1e12
    for a, i in enumerate(rows):
        for b, j in enumerate(columns):
            delta = abs((dates[i]-start[j]).total_seconds())/60.
            # Penalización pequeña para REVISAR; no puede dominar diferencias grandes.
            penal = 2.0 if viable.iloc[j].tipo == 'REVISAR' else 0.0
            cost[a,b] = delta + penal + (b+1)*1e-9
    # No penalizar por superar la tolerancia; solamente vacíos por falta de HEX.
    if n:
        cost[:,m:] = 1e9
    assigned = {}
    if n:
        rr,cc = linear_sum_assignment(cost)
        assigned = {rows[a]:columns[b] for a,b in zip(rr,cc) if b < m}

    suggested = []
    full_pairs = []
    for i, r in df.iterrows():
        if r.archivo_hex:
            suggested.append((r.archivo_hex, 'EXPLICITO', np.nan, '', '', ''))
            continue
        if pd.isna(dates[i]):
            suggested.append(('', 'SIN_FECHA_BITACORA', np.nan, '', '', ''))
            continue
        if i not in assigned:
            suggested.append(('', 'SIN_HEX_DISPONIBLE', np.nan, '', '', ''))
            continue
        j = assigned[i]
        z = viable.iloc[j]
        delta = (dates[i]-start[j]).total_seconds()/60.
        dist = abs(delta)
        ranking = sorted(((abs((dates[i]-start[k]).total_seconds())/60., viable.iloc[k].archivo_hex)
                          for k in columns), key=lambda x:(x[0],x[1]))
        otros = '; '.join(f'{path} ({minutes:.1f} min)' for minutes,path in ranking[:4]
                           if path != z.archivo_hex)
        lim = float(AUTO_CFG.get('tolerancia_min',15))
        if z.tipo != 'PERFIL_CANDIDATO':
            estado = 'REVISAR_TIPO'
        elif dist > lim:
            estado = 'REVISAR_DESFASE'
        elif offset_incompleto:
            estado = 'REVISAR_RELOJ'
        else:
            estado = 'SUGERIDO'
        suggested.append((z.archivo_hex, estado, round(delta,3), z.inicio_hex,
                          otros, z.tipo))
        full_pairs.append((str(r.estacion), z.archivo_hex, round(dist,3), estado))

    reporte = df.copy()
    reporte['archivo_hex'] = [v[0] for v in suggested]
    reporte['estado_asociacion'] = [v[1] for v in suggested]
    reporte['diferencia_inicio_min'] = [v[2] for v in suggested]
    reporte['inicio_hex'] = [v[3] for v in suggested]
    reporte['otros_candidatos'] = [v[4] for v in suggested]
    reporte['tipo_hex'] = [v[5] for v in suggested]
    reporte['reloj_ctd_offset_utc'] = str(AUTO_CFG['hora_ctd_offset_h'])
    reporte['reloj_bitacora_offset_utc'] = str(AUTO_CFG['hora_bitacora_offset_h'])
    reporte['criterio'] = 'Asignación global exclusiva minimizando distancia absoluta a inicio HEX'
    # Lista de archivos sin usar; se mantiene clasificación e información cronológica.
    usados = set(reporte.archivo_hex) - {''}
    sobrantes = inv[~inv.archivo_hex.isin(usados)].copy()
    destino = out/'00_DESCUBRIMIENTO'/'CONCILIACION_PROPUESTA.xlsx'
    destino.parent.mkdir(parents=True,exist_ok=True)
    with pd.ExcelWriter(destino,engine='openpyxl') as writer:
        reporte.to_excel(writer,sheet_name='asignaciones',index=False)
        sobrantes.to_excel(writer,sheet_name='hex_sin_asignar',index=False)
        inv.to_excel(writer,sheet_name='inventario',index=False)
        ws = writer.sheets['asignaciones']
        ws.freeze_panes = 'B2'; ws.auto_filter.ref = ws.dimensions
    print(f'Asignación global exclusiva: {len(full_pairs)} propuestas; '
          f'{len(reporte)-len(full_pairs)-len(explicit)} estaciones sin candidato temporal.')
    return reporte,destino

def preparar_entrada_v23():
    modo=AUTO_CFG.get('modo_entrada')
    if modo is None:
        modo=elegir('Entrada CTD: 1=tabla completa, 2=tabla sin HEX, 3=sin tabla',[
            'Excel/CSV completo, archivo_hex conocido','Excel/CSV con estaciones y horas, sin archivo_hex',
            'Solo HEX y XMLCON, inventario y primera pasada'],1)
    modo=int(modo)
    if modo not in (1,2,3):raise ValueError('modo_entrada debe ser 1, 2 o 3')
    tabla_path=None
    if modo in (1,2):
        tabla_path=_ruta_auto('tabla_estaciones','Selecciona tabla Excel/CSV de estaciones',extensiones=('.xlsx','.xlsm','.csv'))
    hex_dir=_ruta_auto('carpeta_hex','Selecciona carpeta de archivos .hex',carpeta=True)
    xml_path=_ruta_auto('xmlcon','Selecciona XMLCON del equipo',extensiones=('.xmlcon',))
    padre=_ruta_auto('carpeta_salida','Selecciona carpeta de salida',carpeta=True)
    raiz=padre / 'CTD_CORRIDAS_V26'
    raiz.mkdir(parents=True,exist_ok=True)
    marca=datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    out_pre=raiz / ('CTD_' + marca)
    out_pre.mkdir(parents=False, exist_ok=False)
    AUTO_CFG['_salida_corrida'] = out_pre
    print('\nNueva carpeta exclusiva de esta ejecución:',out_pre)
    cfg=leer_config(xml_path)
    print('\nExplorando HEX y verificando la configuración instrumental...')
    estaciones={}
    tabla_inicial=None
    if modo==1:
        tabla_inicial=leer_tabla(tabla_path)
        if 'archivo_hex' in tabla_inicial:
            for _, rr in tabla_inicial.iterrows():
                h=str(rr.get('archivo_hex','')).strip().replace('\\','/')
                if h and h.lower() not in ('nan','none'):
                    if h in estaciones:
                        raise ValueError(f'HEX duplicado en tabla: {h}')
                    estaciones[h]=rr.to_dict()
    inv=_v23_inventariar(hex_dir,cfg,out_pre,guardar_provisional=(modo==3),estaciones=estaciones)
    if modo==3:
        path=_v23_plantilla(inv,out_pre)
        print('\nPRIMERA PASADA TERMINADA (provisional; no es QC/publicación).')
        print('Inventario:',out_pre/'00_DESCUBRIMIENTO'/'inventario_hex.csv')
        print('Complete estación, latitud y longitud en:',path)
        print('Después vuelva a ejecutar en modo 1 con esa plantilla, o modo 2 con su bitácora.')
        return None
    df=tabla_inicial if tabla_inicial is not None else leer_tabla(tabla_path)
    if modo==2 or ('archivo_hex' not in df) or df['archivo_hex'].astype(str).str.strip().eq('').any():
        reporte,propuesta=_v23_conciliar(df,inv,out_pre,hex_dir)
        print('\nPropuesta de conciliación:',propuesta)
        print(reporte[['estacion','archivo_hex','estado_asociacion','diferencia_inicio_min']].to_string(index=False))
        problemas=reporte[reporte.estado_asociacion.str.startswith(('REVISAR','SIN_'))]
        if len(problemas):
            print('\nExisten propuestas que requieren validación (o estaciones sin HEX). Verifique la tabla y ejecútela como modo 1.')
            return None
        if any(reporte.estado_asociacion=='SUGERIDO'):
            if AUTO_CFG['hora_ctd_offset_h'] is None:
                print('AVISO: no se declaró desfase UTC del reloj del CTD. Verifique horas y zona.')
            if not si('¿Validó las asociaciones propuestas con bitácora y desea procesar?',False):
                print('Edite CONCILIACION_PROPUESTA.xlsx y vuelva a ejecutar en modo 1.')
                return None
        df=reporte[list(df.columns)].copy()
        tabla_path=out_pre/'00_DESCUBRIMIENTO'/'ESTACIONES_CONFIRMADAS.xlsx'
        df.to_excel(tabla_path,index=False)
    if 'archivo_hex' not in df or df.archivo_hex.astype(str).str.strip().eq('').any():
        raise ValueError('Faltan asociaciones HEX; revise la tabla.')
    clasificacion = inv.set_index('archivo_hex')['tipo'].to_dict()
    tipos=[]
    for _, r in df.iterrows():
        h=str(r['archivo_hex']).strip().replace('\\','/')
        auto=clasificacion.get(h,'DESCONOCIDO')
        tip=str(r.get('tipo','')).strip().upper()
        if tip.startswith('ESTATICO') and auto!='ESTATICO_CANDIDATO':
            print(f'AVISO: {h} es {auto} pero Excel indica ESTATICO. REVISAR clasificación.')
        tipos.append(auto)
    df=df.copy()
    df['_tipo_inventario']=tipos
    estaticos=df[df['_tipo_inventario']=='ESTATICO_CANDIDATO'].copy()
    pruebas=df[df['_tipo_inventario']=='PRUEBA_AIRE_CANDIDATA'].copy()
    perfiles=df[~df['_tipo_inventario'].isin(['ESTATICO_CANDIDATO','PRUEBA_AIRE_CANDIDATA'])].copy()
    print(f"\nSeparación automática: {len(perfiles)} perfiles/revisión, {len(estaticos)} estáticos, {len(pruebas)} pruebas en aire.")
    print('Los estáticos fueron procesados en 02_ESTATICOS_PROCESADOS con sus metadatos de Excel.')
    if len(estaticos):
        estaticos.drop(columns=['_tipo_inventario']).to_excel(out_pre/'00_DESCUBRIMIENTO'/'ESTATICOS_ASOCIADOS.xlsx',index=False)
    if len(pruebas):
        pruebas.drop(columns=['_tipo_inventario']).to_excel(out_pre/'00_DESCUBRIMIENTO'/'PRUEBAS_EXCLUIDAS.xlsx',index=False)
    if perfiles.empty:
        print('No hay perfiles verticales para enviar al procesador de downcast. Ejecución de estáticos terminada.')
        print('Resultados:',out_pre)
        return None
    df=perfiles.drop(columns=['_tipo_inventario'])
    df,rutas=vincular(df,hex_dir)
    # Guardar tabla exacta para integrarse con v2.2 (vincular ya validó).
    tabla_ok=out_pre/'00_DESCUBRIMIENTO'/'ESTACIONES_V23_VALIDADAS.xlsx'
    df.to_excel(tabla_ok,index=False)
    return tabla_ok,hex_dir,xml_path,padre

def main_auto():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white"
    })

    print(textwrap.dedent(f"""
    CTD CRUCERO AUTO v{VERSION}
    ------------------------------------------------------------
    Modo: inferencia automática + cola de excepciones.
    No pregunta límites, Loop Edit ni adelantos perfil por perfil.
    Calibraciones absolutas: SOLO XMLCON/certificados/referencias externas.
    Correcciones dinámicas: calculadas y documentadas con métricas.
    Producto final: Excel + PNG 300 dpi + logs/diagnósticos de trazabilidad.
    """))

    entradas=preparar_entrada_v23()
    if entradas is None:return None
    tabla_path,hex_dir,xml_path,padre=entradas
    df=leer_tabla(tabla_path)
    df,rutas=vincular(df,hex_dir)
    cfg=leer_config(xml_path)

    invent = pd.DataFrame(
        [[s.indice, s.etiqueta, s.serie, s.fecha_calibracion, s.tipo,
          s.motivo or "Ecuación identificada"] for s in cfg.sensores],
        columns=["Canal XML", "Sensor", "Serie", "Calibración", "Tipo", "Conversión"]
    )

    # Sesión persistente en modo incremental; carpeta con timestamp en modo lote.
    if AUTO_CFG["modo_incremental"]:
        out = AUTO_CFG.get("_salida_corrida", padre / AUTO_CFG["nombre_sesion"])
        out.mkdir(exist_ok=True)
        state, state_path, idx_nuevos, hashes = _preparar_incremental(out, rutas, df)
        if not idx_nuevos:
            print("\nNo hay HEX nuevos o modificados. La sesión ya está al día.")
            print(f"Resultados: {out}")
            return None
        print(f"\nModo incremental: {len(idx_nuevos)} HEX nuevos/cambiados de {len(rutas)} asociados.")
        old_base, old_resumen, old_decisiones = _cargar_base_previa(out)
        # Usa un adelanto C de campaña ya aprendido en ejecuciones anteriores.
        pcam = state.get("parametros_campana", {})
    else:
        out = AUTO_CFG["_salida_corrida"]
        out.mkdir(parents=True,exist_ok=True)
        state, state_path = {"archivos": {}, "parametros_campana": {}}, out / "estado_procesamiento.json"
        idx_nuevos = list(range(len(rutas)))
        hashes = {str(df.iloc[i].get("archivo_hex", rutas[i].name)): sha256(rutas[i]) for i in idx_nuevos}
        old_base, old_resumen, old_decisiones = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        pcam = {}

    diag = out / "diagnosticos_auto"
    diag.mkdir(exist_ok=True)

    print("\nDecodificando HEX nuevos y verificando integridad...")
    lances = []
    df_nuevo = df.iloc[idx_nuevos].reset_index(drop=True)
    rutas_nuevas = [rutas[i] for i in idx_nuevos]

    for j, ruta in enumerate(rutas_nuevas):
        header, raw = leer_hex(ruta, cfg)
        datos, vars = convertir_base(raw, cfg)
        hora, fuente = fecha_header(header)
        l = Lance(
            ruta, df_nuevo.iloc[j].to_dict(), cfg, header, raw, datos, vars,
            cfg.navg / 4, hora, fuente, cfg.avisos.copy()
        )
        lances.append(l)
        print(f"  {ruta.name}: {len(raw)} scans, dt={l.dt:g} s")

    opts = opciones_auto(lances)
    algo_ver = AUTO_CFG["algoritmo_dinamico_version"]
    # Solo reutiliza C si fue declarado estable por ESTA versión del algoritmo.
    prev_C_stable = bool(pcam.get("adelanto_C_estable", False) and pcam.get("algoritmo_version") == algo_ver)
    prev_O2_stable = bool(pcam.get("O2_campana_estable", False) and pcam.get("algoritmo_version") == algo_ver)
    if prev_C_stable and "adelanto_C_s" in pcam:
        opts["adelanto_C_s"] = float(pcam["adelanto_C_s"])
        opts["adelanto_C_fijado"] = True
        opts["C_confianza"] = float(pcam.get("C_confianza", 1.0))
        print(f"Adelanto C estable reutilizado: {opts['adelanto_C_s']:.3f} s")
    if prev_O2_stable and "O2_lag_campana_s" in pcam:
        opts["O2_campana_fijado"] = True
        opts["O2_sensor"] = str(pcam.get("O2_sensor", ""))
        opts["O2_tau_campana"] = bool(pcam.get("O2_tau_campana", False))
        opts["O2_lag_campana_s"] = float(pcam["O2_lag_campana_s"])
        opts["O2_confianza_campana"] = float(pcam.get("O2_confianza_campana", 1.0))
        print(f"O2 campaña reutilizado: Tau={opts['O2_tau_campana']}, lag={opts['O2_lag_campana_s']:.2f} s")

    decisiones_nuevas, diagC, diagO2, diagO2camp = preparar_auto(lances, opts, diag)

    # Persiste únicamente parámetros respaldados por evidencia suficiente.
    pc = state.setdefault("parametros_campana", {})
    pc["algoritmo_version"] = algo_ver
    pc["adelanto_C_s"] = float(opts["adelanto_C_s"])
    pc["C_confianza"] = float(opts.get("C_confianza", 0.0))
    pc["C_metodo"] = opts.get("C_metodo", "")
    pc["adelanto_C_estable"] = bool(
        opts.get("C_metodo") in ("minimo_robusto_campana", "parametro_campana_reutilizado") and
        opts.get("C_confianza", 0) >= 0.60
    )
    pc["O2_sensor"] = str(opts.get("O2_sensor", pcam.get("O2_sensor", "")))
    pc["O2_tau_campana"] = bool(opts.get("O2_tau_campana", pcam.get("O2_tau_campana", False)))
    pc["O2_lag_campana_s"] = float(opts.get("O2_lag_campana_s", pcam.get("O2_lag_campana_s", 0.0)))
    pc["O2_confianza_campana"] = float(opts.get("O2_confianza_campana", pcam.get("O2_confianza_campana", 0.0)))
    pc["O2_campana_estable"] = bool(opts.get("O2_campana_estable", prev_O2_stable))

    newly_stable = ((not prev_C_stable and pc["adelanto_C_estable"]) or
                    (not prev_O2_stable and pc["O2_campana_estable"]))
    if newly_stable and not old_base.empty and len(idx_nuevos) < len(rutas):
        # Los perfiles previos pudieron procesarse con parámetros provisionales.
        # Solo se fuerza otra corrida si esta ejecución NO reprocesó ya toda la campaña.
        state["forzar_reproceso_completo_proxima_ejecucion"] = True

    print("\nRESUMEN DE INFERENCIA AUTOMÁTICA")
    cols = [
        "estacion", "estado_CTD", "estado_O2", "estado_operativo",
        "pmax_dbar", "confianza_cast", "vmed_down_m_s", "vmin_loop_m_s", "retencion_loop",
        "align_C_s", "C_confianza", "O2_advance_effective_s", "O2_mejora_vs_cero_rel",
        "O2_penalizacion_optimo_rel", "O2_confianza", "motivo_revision"
    ]
    print(decisiones_nuevas[cols].to_string(index=False))

    revisar = decisiones_nuevas[decisiones_nuevas.estado == "REVISAR"]
    print(f"\nAUTO_OK: {(decisiones_nuevas.estado == 'AUTO_OK').sum()} | REVISAR: {len(revisar)}")
    if len(revisar):
        print("Solo estos perfiles requieren inspección humana; los demás continúan automáticamente.")

    json_dump(
        out / "plan_automatico_ultima_corrida.json",
        {"ajustes": opts, "lances": {l.ruta.name: l.plan for l in lances}}
    )

    _archivar_manifiesto_previo(out)
    base_nueva, resumen_nuevo, manifiesto = ejecutar(lances, opts, out, tabla_path, df_nuevo)

    # Consolidación persistente.
    base = _merge_por_clave(old_base, base_nueva, "id_perfil")
    resumen = _merge_por_clave(old_resumen, resumen_nuevo, "id")
    decisiones = _merge_por_clave(old_decisiones, decisiones_nuevas, "hex")

    base.to_csv(out / "base_ctd_procesada.csv", index=False, encoding="utf-8-sig")
    resumen.to_csv(out / "resumen_procesamiento.csv", index=False, encoding="utf-8-sig")
    decisiones.to_csv(out / "DECISIONES_AUTOMATICAS.csv", index=False, encoding="utf-8-sig")

    # Validación externa opcional + auditoría de estándares. Nunca modifica calibraciones.
    validacion, resumen_validacion, ref_status = validar_referencias_externas(
        base, AUTO_CFG.get("tabla_referencias"), out
    )
    audit = generar_auditoria_estandares(
        base, resumen, decisiones, cfg, validacion, resumen_validacion, ref_status, out
    )
    methods_path = generar_metodos_publicacion(
        base, decisiones, cfg, audit, resumen_validacion, out
    )
    xlsx = exportar_excel_crucero(
        base, resumen, manifiesto, decisiones, invent, out,
        audit=audit, validacion=validacion, resumen_validacion=resumen_validacion
    )

    # Actualiza estado solo después de que la corrida terminó correctamente.
    for i in idx_nuevos:
        rel = str(df.iloc[i].get("archivo_hex", rutas[i].name))
        state.setdefault("archivos", {})[rel] = hashes[rel]
    state["ultima_ejecucion"] = datetime.now(timezone.utc).isoformat()
    state["xmlcon_sha256"] = sha256(xml_path)
    state["tabla"] = str(tabla_path)
    state["software_version"] = VERSION
    state["algoritmo_dinamico_version"] = AUTO_CFG["algoritmo_dinamico_version"]
    # Si esta corrida ya fue el reproceso forzado y no acaba de estabilizar un
    # parámetro nuevo, la bandera puede retirarse.
    if state.get("forzar_reproceso_completo_proxima_ejecucion") and not newly_stable:
        state["forzar_reproceso_completo_proxima_ejecucion"] = False
    json_dump(state_path, state)

    print(f"\nExcel final consolidado: {xlsx}")
    print(f"Figuras: {out / 'figuras'} (PNG 300 dpi)")
    print(f"Cola de revisión: {out / 'DECISIONES_AUTOMATICAS.csv'}")
    print(f"Estado incremental: {state_path}")
    print(f"Auditoría de estándares: {out / 'AUDITORIA_ESTANDARES_CTD.csv'}")
    print(f"Métodos para publicación: {methods_path}")

    return base, resumen, manifiesto, decisiones


if __name__ == "__main__":
    try:
        resultado = main_auto()
        if resultado is not None:
            base_final, resumen_final, manifiesto_final, decisiones_automaticas = resultado
    except (KeyboardInterrupt, EOFError):
        print("\nEjecución cancelada. Los archivos originales no se modificaron.")
    except Exception as exc:
        print(f"\nNo se completó el procesamiento: {exc}")
        print("Detalle técnico para diagnóstico:")
        traceback.print_exc()
