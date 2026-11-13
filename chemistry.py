"""Transparent aqueous chemistry utilities for HydroGeoAI.
Input concentrations in mg/L; unqualified numeric values only. Not a substitute for lab QA/QC.
"""
import numpy as np
import pandas as pd

IONS = {
    'Ca': (40.078, 2), 'Mg': (24.305, 2), 'Na': (22.989769, 1),
    'K': (39.0983, 1), 'HCO3': (61.0168, -1), 'CO3': (60.0089, -2),
    'Cl': (35.453, -1), 'SO4': (96.06, -2), 'NO3': (62.0049, -1),
}
CATIONS = ['Ca', 'Mg', 'Na', 'K']
ANIONS = ['HCO3', 'CO3', 'Cl', 'SO4', 'NO3']
REQUIRED = ['Sample_ID','Ca','Mg','Na','K','HCO3','Cl','SO4']


def normalize_columns(frame):
    df = frame.copy()
    aliases = {
        'sample':'Sample_ID','sample_id':'Sample_ID','sample id':'Sample_ID',
        'ca2+':'Ca','mg2+':'Mg','na+':'Na','k+':'K','hco3-':'HCO3',
        'co3--':'CO3','cl-':'Cl','so4--':'SO4','no3-':'NO3',
        'sulfate':'SO4','sulphate':'SO4','bicarbonate':'HCO3','chloride':'Cl'
    }
    df.columns = [aliases.get(str(c).strip().lower(), str(c).strip()) for c in df.columns]
    if df.columns.duplicated().any():
        raise ValueError('Duplicate columns after name normalization. Rename duplicate ion columns.')
    return df


def classify_hardness(value):
    if pd.isna(value):
        return 'Unknown'
    if value <= 60:
        return 'Soft'
    if value <= 120:
        return 'Moderately hard'
    if value <= 180:
        return 'Hard'
    return 'Very hard'


def analyze(frame, tolerance=10):
    df = normalize_columns(frame)
    missing = [x for x in REQUIRED if x not in df.columns]
    if missing:
        raise ValueError('Missing required columns: ' + ', '.join(missing))
    if df.empty:
        raise ValueError('Upload at least one sample.')
    if len(df) > 10000:
        raise ValueError('This demo supports at most 10,000 samples per upload.')
    df['Sample_ID'] = df['Sample_ID'].astype('string').str.strip()
    if df['Sample_ID'].isna().any() or (df['Sample_ID'] == '').any():
        raise ValueError('Every row must have a non-empty Sample_ID.')
    if df['Sample_ID'].duplicated().any():
        raise ValueError('Sample_ID values must be unique.')

    for optional in ('CO3', 'NO3'):
        if optional not in df.columns:
            df[optional] = 0.0

    for ion in IONS:
        raw = df[ion]
        converted = pd.to_numeric(raw, errors='coerce')
        invalid = raw.notna() & converted.isna()
        if invalid.any():
            raise ValueError(f'{ion} contains nonnumeric entries (e.g. "<LOD"). Resolve censored data explicitly before analysis; never silently replace with zero.')
        if converted.isna().any():
            raise ValueError(f'{ion} contains missing values. Provide measured concentrations, or document an imputation strategy outside this demo.')
        if (~np.isfinite(converted.to_numpy(dtype=float))).any():
            raise ValueError(f'{ion} contains nonfinite concentrations.')
        if (converted < 0).any():
            raise ValueError(f'{ion} has negative concentrations.')
        df[ion] = converted.astype(float)
        df[f'{ion}_meqL'] = df[ion] * abs(IONS[ion][1]) / IONS[ion][0]

    cat = df[[f'{i}_meqL' for i in CATIONS]].sum(axis=1)
    ani = df[[f'{i}_meqL' for i in ANIONS]].sum(axis=1)
    df['Cations_meqL'] = cat
    df['Anions_meqL'] = ani
    total = cat + ani
    df['CBE_pct'] = np.where(total > 0, 100 * (cat - ani) / total, np.nan)
    df['QA_status'] = np.where(total == 0, 'INSUFFICIENT IONS', np.where(abs(df['CBE_pct']) <= tolerance, 'PASS', 'CHECK'))

    df['Hardness_as_CaCO3_mgL'] = 2.497 * df['Ca'] + 4.118 * df['Mg']
    df['Hardness_class'] = df['Hardness_as_CaCO3_mgL'].apply(classify_hardness)

    df['TDS_est_mgL'] = df[list(IONS)].sum(axis=1)
    df['NaCa_ratio'] = np.where((df['Na_meqL'] + df['Ca_meqL']) > 0, df['Na_meqL'] / (df['Na_meqL'] + df['Ca_meqL']), np.nan)
    df['ClHCO3_ratio'] = np.where((df['Cl_meqL'] + df['HCO3_meqL']) > 0, df['Cl_meqL'] / (df['Cl_meqL'] + df['HCO3_meqL']), np.nan)

    df['Ca_pct'] = np.where(cat > 0, 100 * df['Ca_meqL'] / cat, np.nan)
    df['Mg_pct'] = np.where(cat > 0, 100 * df['Mg_meqL'] / cat, np.nan)
    df['NaK_pct'] = np.where(cat > 0, 100 * (df['Na_meqL'] + df['K_meqL']) / cat, np.nan)
    df['HCO3CO3_pct'] = np.where(ani > 0, 100 * (df['HCO3_meqL'] + df['CO3_meqL']) / ani, np.nan)
    df['Cl_pct'] = np.where(ani > 0, 100 * df['Cl_meqL'] / ani, np.nan)
    df['SO4NO3_pct'] = np.where(ani > 0, 100 * (df['SO4_meqL'] + df['NO3_meqL']) / ani, np.nan)

    def facies(r):
        if r['QA_status'] == 'INSUFFICIENT IONS':
            return 'Unclassified'
        cats = {'Ca': r['Ca_pct'], 'Mg': r['Mg_pct'], 'Na+K': r['NaK_pct']}
        ans = {'HCO3+CO3': r['HCO3CO3_pct'], 'Cl': r['Cl_pct'], 'SO4+NO3': r['SO4NO3_pct']}
        c = max(cats, key=cats.get)
        a = max(ans, key=ans.get)
        return f'{c}–{a}' if cats[c] >= 50 and ans[a] >= 50 else f'Mixed ({c}–{a} dominant)'

    df['Facies'] = df.apply(facies, axis=1)
    return df
