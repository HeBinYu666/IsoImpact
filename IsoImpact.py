import os
import re
import argparse
import gzip

# =====================================================================
# 1. Command-line interface
# =====================================================================
REQUIRED_DOMAIN_COLUMNS = [
    "Protein_ID",
    "Domain_ID",
    "Domain_Name",
    "Genomic_Start",
    "Genomic_End",
]


# Fixed feature schema; order matches Supplementary Table S2 (966 features).
FEATURE_NAMES = [
    "pI", "GRAVY", "Instability_Index",
    "Extinction_Coefficient_Reduced", "Extinction_Coefficient_Oxidized",
    "SecStruct_Helix", "SecStruct_Turn", "SecStruct_Sheet",
]
FEATURE_NAMES += [f"AAC_{aa}" for aa in "ACDEFGHIKLMNPQRSTVWY"]
FEATURE_NAMES += ["Aliphatic_Index"]
for _prop in (
    "Hydrophobicity", "NormVDWVolume", "Polarity", "Polarizability",
    "Charge", "SecondaryStruct", "SolventAccess",
):
    FEATURE_NAMES += [f"CTD_C_{_prop}_G{group}" for group in range(1, 4)]
    FEATURE_NAMES += [f"CTD_T_{_prop}_{pair}" for pair in ("12", "13", "23")]
    for _group in range(1, 4):
        FEATURE_NAMES += [f"CTD_D_{_prop}_G{_group}_{quantile}p" for quantile in (1, 25, 50, 75, 100)]
FEATURE_NAMES += [f"PAAC{i}" for i in range(1, 31)]
FEATURE_NAMES += [f"APAAC{i}" for i in range(1, 41)]
for _method in ("MoreauBrotoAuto", "MoranAuto", "GearyAuto"):
    for _prop in (
        "Hydrophobicity", "AvFlexibility", "Polarizability", "FreeEnergy",
        "ResidueASA", "ResidueVol", "Steric", "Mutability",
    ):
        FEATURE_NAMES += [f"{_method}_{_prop}{lag}" for lag in range(1, 31)]
if len(FEATURE_NAMES) != 966 or len(set(FEATURE_NAMES)) != 966:
    raise RuntimeError("The fixed feature schema must contain 966 unique features.")


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

parser = argparse.ArgumentParser(
    description="IsoImpact: automated isoform functional impact analysis"
)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    "-i", dest="isoform", nargs="+",
    help="Compare at least two isoform IDs, for example: -i ENST1 ENST2"
)
group.add_argument(
    "-r", dest="reference_pair", nargs=2,
    metavar=("ALTERNATIVE1_ISOFORM", "ALTERNATIVE2_ISOFORM"),
    help="Compare two isoforms in the given order as alternative2 - alternative1"
)

parser.add_argument("-g", dest="gtf", required=True, help="Path to the matching GTF annotation file")
parser.add_argument("-f", dest="fasta", required=True, help="Path to the matching protein FASTA file")
parser.add_argument(
    "-d", dest="domain", required=True,
    help="Path to the matching domain-coordinate CSV file"
)
parser.add_argument(
    "-o", dest="outdir", default=".",
    help="Directory for output files"
)
parser.add_argument(
    "-z", dest="difference_view", action="store_true",
    help="Also write a difference-focused transcript-structure figure"
)

args = parser.parse_args()

GTF_FILE = args.gtf
FASTA_FILE = args.fasta
OUT_DIR = args.outdir

try:
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.gridspec import GridSpec
    import numpy as np
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError as exc:
    parser.error(
        "Missing Python dependency. Install the required packages with "
        "`pip install -r requirements.txt`. Original error: " + str(exc)
    )

# =====================================================================
# 2. Domain annotation input
# =====================================================================
DOMAIN_CSV = args.domain
print(f"[IsoImpact] Using domain-coordinate file: {DOMAIN_CSV}")

if not os.path.exists(GTF_FILE):
    parser.error(f"GTF file not found: {GTF_FILE}")

if not os.path.exists(FASTA_FILE):
    parser.error(f"FASTA file not found: {FASTA_FILE}")

if not os.path.exists(DOMAIN_CSV):
    parser.error(f"Domain-coordinate CSV not found: {DOMAIN_CSV}")

df_domains = pd.read_csv(DOMAIN_CSV)
missing_domain_cols = [c for c in REQUIRED_DOMAIN_COLUMNS if c not in df_domains.columns]
if missing_domain_cols:
    parser.error(
        "Domain-coordinate CSV is missing required column(s): "
        + ", ".join(missing_domain_cols)
    )

df_domains["Protein_ID"] = df_domains["Protein_ID"].astype(str)

# =====================================================================
# 3. Isoform comparison mode
# =====================================================================
if args.isoform:
    INPUT_TXS = args.isoform
    COMPARISON_MODE = "canonical"
    if len(INPUT_TXS) < 2:
        parser.error("-i requires at least two isoform IDs.")
else:
    INPUT_TXS = args.reference_pair
    COMPARISON_MODE = "reference_pair"
    if INPUT_TXS[0] == INPUT_TXS[1]:
        parser.error("-r requires two different isoform IDs.")

os.makedirs(OUT_DIR, exist_ok=True)

# =====================================================================
# 4. Protein sequence feature extraction
# =====================================================================
def get_high_dim_propy3_features(seq):
    if not seq or len(seq) < 35:
        return {}
    feats = {}
    try:
        from propy.PyPro import GetProDes
        des = GetProDes(seq)
    except Exception:
        return feats
    for method_name, method in (
        ("PAAC", des.GetPAAC),
        ("APAAC", des.GetAPAAC),
        ("MoreauBrotoAuto", des.GetMoreauBrotoAuto),
        ("MoranAuto", des.GetMoranAuto),
        ("GearyAuto", des.GetGearyAuto),
    ):
        try:
            feats.update(method())
        except Exception:
            # The fixed schema below records any unavailable descriptors as NA.
            pass
    return feats


def calculate_ctd(seq):
    if not seq or len(seq) == 0: return {}
    properties = {
        'Hydrophobicity': [{'R','K','E','D','Q','N'}, {'G','A','S','T','P','H','Y'}, {'C','V','L','I','M','F','W'}],
        'NormVDWVolume': [{'G','A','S','T','P','D','C'}, {'N','V','E','Q','I','L'}, {'M','H','K','F','R','Y','W'}],
        'Polarity': [{'L','I','F','W','C','M','V','Y'}, {'P','A','T','G','S'}, {'H','Q','R','K','N','E','D'}],
        'Polarizability': [{'G','A','S','D','T'}, {'C','P','N','V','E','Q','I','L'}, {'K','M','H','F','R','Y','W'}],
        'Charge': [{'K','R'}, {'A','N','C','Q','G','I','L','M','F','P','S','T','W','Y','V','H'}, {'D','E'}],
        'SecondaryStruct': [{'E','A','L','M','Q','K','R','H'}, {'V','I','Y','C','W','F','T'}, {'G','N','P','S','D'}],
        'SolventAccess': [{'A','L','F','C','G','I','V','W'}, {'R','K','Q','E','N','D'}, {'M','S','P','T','H','Y'}]
    }
    ctd = {}
    seq_len = len(seq)
    for prop_name, groups in properties.items():
        counts = [sum(1 for aa in seq if aa in g) for g in groups]
        for i in range(3): ctd[f"CTD_C_{prop_name}_G{i+1}"] = round(counts[i] / seq_len, 4)
        t_counts = [0, 0, 0]
        for i in range(seq_len - 1):
            g1 = next((j for j, g in enumerate(groups) if seq[i] in g), -1)
            g2 = next((j for j, g in enumerate(groups) if seq[i+1] in g), -1)
            if g1 != -1 and g2 != -1 and g1 != g2:
                if (g1==0 and g2==1) or (g1==1 and g2==0): t_counts[0] += 1
                elif (g1==0 and g2==2) or (g1==2 and g2==0): t_counts[1] += 1
                elif (g1==1 and g2==2) or (g1==2 and g2==1): t_counts[2] += 1
        ctd[f"CTD_T_{prop_name}_12"] = round(t_counts[0] / (seq_len-1), 4) if seq_len>1 else 0
        ctd[f"CTD_T_{prop_name}_13"] = round(t_counts[1] / (seq_len-1), 4) if seq_len>1 else 0
        ctd[f"CTD_T_{prop_name}_23"] = round(t_counts[2] / (seq_len-1), 4) if seq_len>1 else 0
        for i, g in enumerate(groups):
            pos_list = [idx+1 for idx, aa in enumerate(seq) if aa in g]
            for percent in [1, 25, 50, 75, 100]:
                feat_name = f"CTD_D_{prop_name}_G{i+1}_{percent}p"
                if len(pos_list) == 0: ctd[feat_name] = 0
                else:
                    target_idx = max(0, int((percent/100) * len(pos_list)) - 1)
                    if percent == 1: target_idx = 0
                    if percent == 100: target_idx = len(pos_list) - 1
                    ctd[feat_name] = round(pos_list[target_idx] / seq_len, 4)
    return ctd


def get_all_features(seq):
    # Always return the complete fixed schema. Unavailable values remain NA.
    feats = {name: float("nan") for name in FEATURE_NAMES}
    feats["MW"] = float("nan")
    if not seq:
        return feats
    seq = seq.replace("X", "").replace("U", "").replace("Z", "").replace("B", "")
    if not seq:
        return feats
    try:
        pa = ProteinAnalysis(seq)
        feats["MW"] = round(pa.molecular_weight(), 2)
        feats["pI"] = round(pa.isoelectric_point(), 4)
        feats["GRAVY"] = round(pa.gravy(), 4)
        feats["Instability_Index"] = round(pa.instability_index(), 4)
        extinction = pa.molar_extinction_coefficient()
        feats["Extinction_Coefficient_Reduced"] = extinction[0]
        feats["Extinction_Coefficient_Oxidized"] = extinction[1]
        sec = pa.secondary_structure_fraction()
        feats["SecStruct_Helix"], feats["SecStruct_Turn"], feats["SecStruct_Sheet"] = (
            round(sec[0], 4), round(sec[1], 4), round(sec[2], 4)
        )
        aac = pa.amino_acids_percent
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            if aa in aac:
                feats[f"AAC_{aa}"] = round(aac[aa], 4)
        required_aas = [aac.get(aa) for aa in "AVIL"]
        if all(value is not None for value in required_aas):
            mole_pct_A, mole_pct_V, mole_pct_I, mole_pct_L = [value * 100 for value in required_aas]
            feats["Aliphatic_Index"] = round(
                mole_pct_A + 2.9 * mole_pct_V + 3.9 * (mole_pct_I + mole_pct_L), 4
            )
    except Exception:
        # Keep any successfully initialized values and leave failed features as NA.
        pass
    try:
        for key, value in calculate_ctd(seq).items():
            if key in feats:
                feats[key] = value
    except Exception:
        pass
    for key, value in get_high_dim_propy3_features(seq).items():
        if key in feats:
            feats[key] = value
    return feats


def signed_normalized_difference(compared_value, reference_value):
    try:
        compared_value = float(compared_value)
        reference_value = float(reference_value)
    except (TypeError, ValueError, OverflowError):
        return float('nan')
    # Relative change is measured against the canonical/reference isoform.
    # A zero reference value has no finite relative-change denominator.
    if not np.isfinite(compared_value) or not np.isfinite(reference_value) or reference_value == 0:
        return float('nan')
    relative_change = (compared_value - reference_value) / reference_value
    return relative_change if np.isfinite(relative_change) else float('nan')

# =====================================================================
# 5. Parse the GTF annotation
# =====================================================================
tx_info = {tx: {'exons': [], 'cds': [], 'strand': '+', 'is_canonical': False, 
                'protein_id': '', 'gene_id': '', 'biotype': 'unknown',
                'span_start': float('inf'), 'span_end': 0} for tx in INPUT_TXS}

if os.path.exists(GTF_FILE):
    with open_text(GTF_FILE) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts) < 9: continue
            info_str = parts[8]
            for tx in INPUT_TXS:
                if tx in info_str:
                    tx_info[tx]['strand'] = parts[6]
                    if 'Ensembl_canonical' in info_str or 'MANE_Select' in info_str: tx_info[tx]['is_canonical'] = True
                    match_prot = re.search(r'protein_id "([^"]+)"', info_str)
                    if match_prot: tx_info[tx]['protein_id'] = match_prot.group(1)
                    match_gene = re.search(r'gene_id "([^"]+)"', info_str)
                    if match_gene: tx_info[tx]['gene_id'] = match_gene.group(1)
                    match_biotype = re.search(r'transcript_biotype "([^"]+)"', info_str)
                    if not match_biotype: match_biotype = re.search(r'transcript_type "([^"]+)"', info_str)
                    if match_biotype: tx_info[tx]['biotype'] = match_biotype.group(1)
                    
                    start, end = int(parts[3]), int(parts[4])
                    if parts[2] == 'exon':
                        tx_info[tx]['span_start'] = min(tx_info[tx]['span_start'], start)
                        tx_info[tx]['span_end'] = max(tx_info[tx]['span_end'], end)
                        tx_info[tx]['exons'].append({'start': start, 'end': end})
                    elif parts[2] == 'CDS':
                        tx_info[tx]['cds'].append({'start': start, 'end': end})

missing_txs = [tx for tx in INPUT_TXS if tx_info[tx]['span_start'] == float('inf')]
if missing_txs:
    parser.error(
        "The following isoform ID(s) were not found in the GTF exon records: "
        + ", ".join(missing_txs)
    )

missing_proteins = [tx for tx in INPUT_TXS if not tx_info[tx]['protein_id']]
if missing_proteins:
    parser.error(
        "The following isoform ID(s) do not have protein_id annotations in the GTF: "
        + ", ".join(missing_proteins)
    )

if COMPARISON_MODE == "reference_pair":
    gene_ids = [tx_info[tx]['gene_id'] for tx in INPUT_TXS]
    if any(not gene_id for gene_id in gene_ids) or len(set(gene_ids)) != 1:
        parser.error("-r requires two isoforms from the same gene.")
    can_tx = INPUT_TXS[0]
    alt_txs = [INPUT_TXS[1]]
    reference_label = "alternative1"
    compared_label = "alternative2"
    reference_feature_prefix = "alternative1"
    compared_feature_prefix = "alternative2"
else:
    can_tx = next((tx for tx in INPUT_TXS if tx_info[tx]['is_canonical']), None)
    if not can_tx:
        can_tx = max(INPUT_TXS, key=lambda x: tx_info[x]['span_end'] - tx_info[x]['span_start'])
    alt_txs = [tx for tx in INPUT_TXS if tx != can_tx]
    reference_label = "Canonical"
    compared_label = "Alternative"
    reference_feature_prefix = "Can"
    compared_feature_prefix = "Alt"

# =====================================================================
# 6. Extract protein FASTA sequences and build the feature matrix
# =====================================================================
target_prots = [tx_info[tx]['protein_id'] for tx in INPUT_TXS if tx_info[tx]['protein_id']]
sequences = {pid: "" for pid in target_prots}
current_id = None
if os.path.exists(FASTA_FILE):
    with open_text(FASTA_FILE) as f:
        for line in f:
            if line.startswith('>'):
                current_id = next((pid for pid in target_prots if pid in line), None)
            elif current_id: sequences[current_id] += line.strip()


missing_sequences = [pid for pid, seq in sequences.items() if not seq]
if missing_sequences:
    parser.error(
        "The following protein sequence(s) were not found in the FASTA file: "
        + ", ".join(missing_sequences)
    )

can_prot_id = tx_info[can_tx]['protein_id']
feat_c = get_all_features(sequences.get(can_prot_id, ""))

d_c_list = df_domains[df_domains['Protein_ID'].str.contains(can_prot_id, na=False, regex=False)]['Domain_Name'].tolist() if can_prot_id and not df_domains.empty else []
sp_c = tx_info[can_tx]['span_end'] - tx_info[can_tx]['span_start'] + 1

rows = []
normalized_rows = []
for alt_tx in alt_txs:
    alt_prot_id = tx_info[alt_tx]['protein_id']
    feat_a = get_all_features(sequences.get(alt_prot_id, ""))
    
    d_a_list = df_domains[df_domains['Protein_ID'].str.contains(alt_prot_id, na=False, regex=False)]['Domain_Name'].tolist() if alt_prot_id and not df_domains.empty else []
    d_shared = sorted(set(d_c_list) & set(d_a_list))
    d_lost = sorted(set(d_c_list) - set(d_a_list))
    d_gained = sorted(set(d_a_list) - set(d_c_list))
    sp_a = tx_info[alt_tx]['span_end'] - tx_info[alt_tx]['span_start'] + 1
    mw_c = feat_c.get("MW", np.nan)
    mw_a = feat_a.get("MW", np.nan)
    mw_delta = round(mw_a - mw_c, 2) if np.isfinite(mw_a) and np.isfinite(mw_c) else np.nan
    row = {
        "Gene_ID": tx_info[can_tx]['gene_id'],
        f"{reference_label}_Transcript_ID": can_tx,
        f"{compared_label}_Transcript_ID": alt_tx,
        f"{reference_label}_Protein_ID": can_prot_id,
        f"{compared_label}_Protein_ID": alt_prot_id,
        f"Coding_Potential_{reference_label}": tx_info[can_tx]['biotype'],
        f"Coding_Potential_{compared_label}": tx_info[alt_tx]['biotype'],
        f"Protein_Length_{reference_label}": len(sequences.get(can_prot_id, "")),
        f"Protein_Length_{compared_label}": len(sequences.get(alt_prot_id, "")),
        "Protein_Length_Difference": len(sequences.get(alt_prot_id, "")) - len(sequences.get(can_prot_id, "")),
        f"Molecular_Weight_{reference_label}_Da": mw_c,
        f"Molecular_Weight_{compared_label}_Da": mw_a,
        "Molecular_Weight_Difference_Da": mw_delta,
        f"Exon_Count_{reference_label}": len(tx_info[can_tx]['exons']),
        f"Exon_Count_{compared_label}": len(tx_info[alt_tx]['exons']),
        "Exon_Count_Difference": len(tx_info[alt_tx]['exons']) - len(tx_info[can_tx]['exons']),
        f"Genomic_Span_{reference_label}_bp": sp_c if sp_c != float('inf') else 0,
        f"Genomic_Span_{compared_label}_bp": sp_a if sp_a != float('inf') else 0,
        "Genomic_Span_Difference_bp": (sp_a - sp_c) if sp_c != float('inf') else 0,
        f"Domains_List_{reference_label}": ";".join(d_c_list) if d_c_list else "None",
        f"Domains_List_{compared_label}": ";".join(d_a_list) if d_a_list else "None",
        "Domains_List_Shared": ";".join(d_shared) if d_shared else "None",
        "Domains_List_Lost": ";".join(d_lost) if d_lost else "None",
        "Domains_List_Gained": ";".join(d_gained) if d_gained else "None",
        "Domain_Count_Lost": len(d_lost),
        "Domain_Count_Gained": len(d_gained),
        "Domain_Count_Total_Change": len(d_lost) + len(d_gained)
    }
    normalized_values = {"MW": signed_normalized_difference(mw_a, mw_c)}
    for k in FEATURE_NAMES:
        reference_value = feat_c.get(k, np.nan)
        compared_value = feat_a.get(k, np.nan)
        row[f"{reference_feature_prefix}_{k}"] = reference_value
        row[f"{compared_feature_prefix}_{k}"] = compared_value
        if np.isfinite(reference_value) and np.isfinite(compared_value):
            row[f"Delta_{k}"] = round(compared_value - reference_value, 4)
        else:
            row[f"Delta_{k}"] = np.nan
        normalized_values[k] = signed_normalized_difference(compared_value, reference_value)
    rows.append(row)
    normalized_rows.append(normalized_values)

df_all_features = pd.DataFrame(rows)
if df_all_features.empty:
    print("Error: no isoform comparison rows were generated. Please check isoform IDs, GTF, FASTA, and domain files.")
    exit(1)

normalized_deltas = pd.DataFrame(
    [
        {
            "Normalized_Delta_MW": values["MW"],
            **{f"Normalized_Delta_{name}": values[name] for name in FEATURE_NAMES},
        }
        for values in normalized_rows
    ]
)
normalized_delta_cols = list(normalized_deltas.columns)
if len(normalized_delta_cols) != 967:
    parser.error("The in-memory normalized feature panel must contain 967 values.")

# Normalized_Delta_* values are retained for ranking and plotting only; they are
# intentionally not written to the feature CSV.
csv_filename = os.path.join(OUT_DIR, "IsoImpact_features.csv")
df_all_features.to_csv(csv_filename, index=False, na_rep='NA')
print(f"[IsoImpact] Feature matrix written to {csv_filename} ({len(rows[0])} columns).")

# =====================================================================
# 7. Dynamic visualization layout
# =====================================================================
print("[IsoImpact] Rendering domain map and ranked feature-difference plot...")

exons_dict = {tx: sorted(tx_info[tx]['exons'], key=lambda x: x['start']) for tx in INPUT_TXS}
cds_dict = {tx: sorted(tx_info[tx]['cds'], key=lambda x: x['start']) for tx in INPUT_TXS}

boundaries = set()
for ex_list in exons_dict.values():
    for e in ex_list: boundaries.update([e['start'], e['end']])
for c_list in cds_dict.values():
    for c in c_list: boundaries.update([c['start'], c['end']])
boundaries = sorted(list(boundaries))

segments = []
if boundaries:
    for i in range(len(boundaries)-1):
        s, e = boundaries[i], boundaries[i+1]
        is_cds = any(c['start'] <= s and e <= c['end'] for c_list in cds_dict.values() for c in c_list)
        is_exon = any(ex['start'] <= s and e <= ex['end'] for ex_list in exons_dict.values() for ex in ex_list)
        segments.append({'start': s, 'end': e, 'is_cds': is_cds, 'is_exon': is_exon})

MAX_INTRON_LEN = 200
MAX_UTR_LEN = 150
mapping = {boundaries[0]: 0} if boundaries else {}
current_vis_pos = 0

for seg in segments:
    length = seg['end'] - seg['start']
    if not seg['is_exon']: vis_len = min(length, MAX_INTRON_LEN)
    elif not seg['is_cds']: vis_len = min(length, MAX_UTR_LEN)
    else: vis_len = length         
    current_vis_pos += vis_len
    mapping[seg['end']] = current_vis_pos

def map_pos(pos):
    if not boundaries: return 0
    if pos <= boundaries[0]: return mapping[boundaries[0]]
    if pos >= boundaries[-1]: return mapping[boundaries[-1]]
    for i in range(len(boundaries)-1):
        if boundaries[i] <= pos <= boundaries[i+1]:
            s, e = boundaries[i], boundaries[i+1]
            if s == e: return mapping[s]
            ratio = (pos - s) / (e - s)
            return mapping[s] + ratio * (mapping[e] - mapping[s])
    return 0

# Prepare plot data and domain colors.
all_txs_to_plot = [can_tx] + alt_txs
current_used_domains = []
for tx in all_txs_to_plot:
    prot_id = tx_info[tx]['protein_id']
    if prot_id and not df_domains.empty:
        tx_doms = df_domains[df_domains['Protein_ID'].str.contains(prot_id, na=False, regex=False)]
        current_used_domains.extend(tx_doms['Domain_Name'].dropna().tolist())

unique_domains = sorted(list(set(current_used_domains)))
DOMAIN_PALETTE = [
    '#F4A460',
    '#5DCAA5',
    '#9B8DC4',
    '#E8735A',
    '#87CEEB',
    '#FFB347',
    '#77DD77',
    '#DDA0DD',
    '#48D1CC',
    '#F08080',
    '#AEC6CF',
    '#C8A96A',
    '#6FB8A3',
    '#E0A080',
    '#A8C8A0',
    '#C4A0C4',
    '#80C4C4',
    '#D4B080',
    '#90B8D0',
    '#C8A090',
]
dynamic_domain_colors = {dom: DOMAIN_PALETTE[i % len(DOMAIN_PALETTE)] 
                         for i, dom in enumerate(unique_domains)}

# Size A independently of the B grid so panel labels and the domain legend
# cannot squeeze the transcript tracks or create large gaps between panels.
alternative_count = max(len(alt_txs), 1)
panel_columns = min(3, alternative_count)
panel_rows = int(np.ceil(alternative_count / panel_columns))
legacy_height = max(12, 4 + len(INPUT_TXS) * 2.0)
track_height_ratio = len(INPUT_TXS) * 0.8
# Original A-axis height, measured with the legacy 20-inch-wide layout.
a_height = (legacy_height - 1.62223) * track_height_ratio / (track_height_ratio + 3.0)
b_height = 3.2
b_bottom_margin = 1.35
b_row_pitch = b_height + 1.8
a_bottom = b_bottom_margin + b_height + (panel_rows - 1) * b_row_pitch + 0.8
fig_height = a_bottom + a_height + 0.5
fig = plt.figure(figsize=(20, fig_height), dpi=300, constrained_layout=False)

# ==================== [A] Genomic structure and domain map ====================
ax_a = fig.add_axes([2.9 / 20, a_bottom / fig_height, 13.1 / 20, a_height / fig_height])
ax_a.set_title('A', loc='left', fontsize=22, fontweight='bold', pad=10)

H_utr, H_cds, H_domain = 0.20, 0.40, 0.55

# Use one shared left label anchor to keep isoform labels aligned.
global_min_x = float('inf')
for tx in all_txs_to_plot:
    if exons_dict[tx]: global_min_x = min(global_min_x, map_pos(exons_dict[tx][0]['start']))
text_align_x = global_min_x - 400 if global_min_x != float('inf') else -400

y_positions = np.arange(len(all_txs_to_plot), 0, -1) * 2.0 
track_labels = [f"{reference_label}:\n{can_tx}"] + [f"{compared_label}:\n{tx}" for tx in alt_txs]

for tx, y, label in zip(all_txs_to_plot, y_positions, track_labels):
    exons = exons_dict[tx]
    cds_list = cds_dict[tx]
    if not exons: continue
    
    vis_gene_start, vis_gene_end = map_pos(exons[0]['start']), map_pos(exons[-1]['end'])
    ax_a.plot([vis_gene_start, vis_gene_end], [y, y], color='dimgray', linewidth=1.5, zorder=1)
    
    ax_a.text(text_align_x, y, label, ha='right', va='center', fontweight='bold', color='#222222', fontsize=20)

    for i in range(len(exons) - 1):
        intron_start = map_pos(exons[i]['end'])
        intron_end = map_pos(exons[i+1]['start'])
        intron_len = intron_end - intron_start
        if intron_len >= 40: 
            if intron_len < 120: ax_a.plot(intron_start + intron_len/2, y, marker='>', markersize=6, color='#888888', linestyle='None', zorder=2) 
            else:
                ax_a.plot(intron_start + intron_len/3, y, marker='>', markersize=6, color='#888888', linestyle='None', zorder=2) 
                ax_a.plot(intron_start + 2*intron_len/3, y, marker='>', markersize=6, color='#888888', linestyle='None', zorder=2) 

    for exon in exons:
        v_start, v_end = map_pos(exon['start']), map_pos(exon['end'])
        ax_a.add_patch(patches.Rectangle((v_start, y - H_utr/2), v_end - v_start, H_utr, linewidth=1.0, edgecolor='gray', facecolor='#e0e0e0', zorder=3)) 

    for cds in cds_list:
        v_start, v_end = map_pos(cds['start']), map_pos(cds['end'])
        ax_a.add_patch(patches.Rectangle((v_start, y - H_cds/2), v_end - v_start, H_cds, linewidth=1.2, edgecolor='#2C5F8A' , facecolor='#5B8DB8', zorder=4)) 
        
    if not df_domains.empty and tx_info[tx]['protein_id']:
        tx_domains = df_domains[df_domains['Protein_ID'].str.contains(tx_info[tx]['protein_id'], na=False, regex=False)]
        for _, row_domain in tx_domains.iterrows():
            d_start, d_end = row_domain['Genomic_Start'], row_domain['Genomic_End'] 
            color = dynamic_domain_colors.get(row_domain['Domain_Name'], "orange") 
            for cds in cds_list:
                overlap_start = max(cds['start'], d_start) 
                overlap_end = min(cds['end'], d_end) 
                if overlap_start < overlap_end:
                    v_start, v_end = map_pos(overlap_start), map_pos(overlap_end) 
                    ax_a.add_patch(patches.Rectangle((v_start, y - H_domain/2), v_end - v_start, H_domain, linewidth=1.0, edgecolor='black', facecolor=color, alpha=1.0, zorder=5)) 

ax_a.set_ylim(0.5, len(all_txs_to_plot) * 2.0 + 1.0) 
ax_a.axis('off')

legend_elements = [patches.Patch(facecolor='#e0e0e0', edgecolor='#999999', label='UTR'), 
                   patches.Patch(facecolor='#5B8DB8', edgecolor='#2C5F8A', label='CDS')]
for d_name, color in dynamic_domain_colors.items(): 
    legend_elements.append(patches.Patch(facecolor=color, edgecolor='dimgray', label=d_name)) 
ax_a.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=False, 
            prop={'size': 20, 'weight': 'bold'}, title="Regions & Domains", title_fontproperties={'size': 22, 'weight': 'bold'})

# ==================== [B] Top-15 relative feature differences ====================
# Selection uses absolute relative change; the selected features are then
# ordered by signed value to match the paper's Figure 2B panels.
def abbreviated_feature_name(column_name, max_len=8):
    name = column_name.replace('Normalized_Delta_', '')
    name = re.sub(r'\s+', '_', name.strip())
    match = re.fullmatch(
        r'(MoreauBrotoAuto|MoranAuto|GearyAuto)_(.+?)(\d+)', name
    )
    if match:
        method, property_name, lag = match.groups()
        method_short = {'MoreauBrotoAuto': 'MB', 'MoranAuto': 'MO', 'GearyAuto': 'GE'}
        property_short = {
            'Hydrophobicity': 'H', 'AvFlexibility': 'F', 'Polarizability': 'P',
            'FreeEnergy': 'E', 'ResidueASA': 'A', 'ResidueVol': 'V',
            'Steric': 'S', 'Mutability': 'M',
        }
        if property_name in property_short:
            return f"{method_short[method]}.{property_short[property_name]}{lag}"
    replacements = [
        ('Extinction_Coefficient_Oxidized', 'ExtCoef_Ox'),
        ('Extinction_Coefficient_Reduced', 'ExtCoef_Red'),
        ('Aliphatic_Index', 'Aliphatic'),
        ('Instability_Index', 'Instability'),
        ('MoreauBrotoAuto', 'MB'), ('MoranAuto', 'Moran'), ('GearyAuto', 'Geary'),
        ('SecStruct_Helix', 'Helix'), ('SecStruct_Turn', 'Turn'),
        ('SecStruct_Sheet', 'Sheet'), ('Hydrophobicity', 'Hydrophob.'),
        ('Polarizability', 'Polariz.'), ('SecondaryStruct', 'SecStruct'),
        ('SolventAccess', 'SolvAcc'), ('NormVDWVolume', 'VDWVol'),
        ('Composition', 'Comp.'), ('Transition', 'Trans.'), ('Distribution', 'Dist.'),
    ]
    for old, new in replacements:
        name = name.replace(old, new)
    name = re.sub(r'CTD_([CTD])_', r'CTD\1_', name)
    return name if len(name) <= max_len else name[:max_len - 3] + '...'

if not df_all_features.empty:
    plot_limit = 50.0
    plot_records = []
    for i, (_, normalized_row) in enumerate(normalized_deltas.iterrows()):
        metadata_row = df_all_features.iloc[i]
        alt_id = metadata_row[f'{compared_label}_Transcript_ID']
        row_delta = pd.to_numeric(normalized_row, errors='coerce')
        row_delta = row_delta[np.isfinite(row_delta)]
        selected_columns = (
            row_delta.abs()
            .sort_values(ascending=False, kind='mergesort')
            .head(15)
            .index
            .tolist()
        )
        selected_columns = sorted(
            selected_columns,
            key=lambda column: float(row_delta[column]),
            reverse=True,
        )
        plot_records.append(
            (
                alt_id,
                [abbreviated_feature_name(column) for column in selected_columns],
                row_delta.loc[selected_columns].to_numpy(dtype=float),
            )
        )

    for panel_index in range(panel_rows * panel_columns):
        panel_row = panel_index // panel_columns
        panel_column = panel_index % panel_columns
        if panel_index >= len(plot_records):
            continue
        b_left = 2.9 if panel_columns == 1 else 1.5
        b_total_width = 13.1 if panel_columns == 1 else 18.1
        b_gap = 0.6
        b_width = (b_total_width - b_gap * (panel_columns - 1)) / panel_columns
        b_bottom = b_bottom_margin + (panel_rows - panel_row - 1) * b_row_pitch
        ax_b = fig.add_axes([
            (b_left + panel_column * (b_width + b_gap)) / 20,
            b_bottom / fig_height, b_width / 20, b_height / fig_height,
        ])

        alt_id, feature_names, y_vals = plot_records[panel_index]
        x_indices = np.arange(len(y_vals))
        plot_values = np.clip(y_vals, -plot_limit, plot_limit)
        style_color = '#2F6F73'
        ax_b.vlines(
            x_indices, 0, plot_values, color=style_color,
            alpha=0.78, linewidth=8.4, zorder=3,
        )
        for mask, marker, boundary in (
            (y_vals > plot_limit, '^', plot_limit),
            (y_vals < -plot_limit, 'v', -plot_limit),
        ):
            if mask.any():
                ax_b.scatter(
                    x_indices[mask], np.full(mask.sum(), boundary),
                    marker=marker, s=44, color=style_color,
                    clip_on=False, zorder=4,
                )
        ax_b.axhline(0, color='black', linewidth=1.2, zorder=2)
        ax_b.set_ylim(-plot_limit, plot_limit)
        ax_b.set_yticks(np.arange(-50, 51, 25))
        ax_b.set_xlim(-0.6, max(len(y_vals) - 0.4, 0.6))
        ax_b.set_xticks(x_indices)
        ax_b.set_xticklabels(
            feature_names, rotation=65, ha='right', rotation_mode='anchor',
            fontsize=18, fontweight='bold',
        )
        ax_b.tick_params(axis='y', labelsize=18, labelleft=(panel_column == 0))
        for tick_label in ax_b.get_yticklabels():
            tick_label.set_fontweight('bold')
        ax_b.grid(axis='y', color='#D9D9D9', linewidth=0.6, alpha=0.85)
        ax_b.spines['top'].set_visible(False)
        ax_b.spines['right'].set_visible(False)
        ax_b.set_title(
            str(alt_id), fontsize=18, fontweight='bold', pad=8,
        )
        if panel_column == 0:
            ax_b.set_ylabel('Relative feature change', fontsize=20, fontweight='bold', labelpad=8)
        else:
            ax_b.set_ylabel('')
        if panel_index == 0:
            fig.text(
                b_left / 20, (b_bottom + b_height + 0.45) / fig_height,
                'B', fontsize=22, fontweight='bold', ha='left', va='bottom',
            )

output_filename = os.path.join(OUT_DIR, "IsoImpact_figure.png")
fig.savefig(output_filename, dpi=300, bbox_inches='tight', facecolor='white')

if args.difference_view:
    def interval_overlaps(items, start, end):
        return any(item['start'] < end and item['end'] > start for item in items)

    difference_segments = []
    for segment_start, segment_end in zip(boundaries[:-1], boundaries[1:]):
        reference_has_cds = interval_overlaps(cds_dict[can_tx], segment_start, segment_end)
        compared_cds_states = [
            interval_overlaps(cds_dict[tx], segment_start, segment_end)
            for tx in alt_txs
        ]
        if any(state != reference_has_cds for state in compared_cds_states):
            difference_segments.append((map_pos(segment_start), map_pos(segment_end)))

    if not difference_segments:
        for segment_start, segment_end in zip(boundaries[:-1], boundaries[1:]):
            reference_has_exon = interval_overlaps(exons_dict[can_tx], segment_start, segment_end)
            compared_exon_states = [
                interval_overlaps(exons_dict[tx], segment_start, segment_end)
                for tx in alt_txs
            ]
            if any(state != reference_has_exon for state in compared_exon_states):

                difference_segments.append((map_pos(segment_start), map_pos(segment_end)))

    full_visible_start = map_pos(boundaries[0])
    full_visible_end = map_pos(boundaries[-1])
    merged_difference_groups = []
    for segment_start, segment_end in difference_segments:
        if not merged_difference_groups or segment_start - merged_difference_groups[-1][1] > MAX_INTRON_LEN * 1.25:
            merged_difference_groups.append([segment_start, segment_end, segment_end - segment_start])
        else:
            merged_difference_groups[-1][1] = segment_end
            merged_difference_groups[-1][2] += segment_end - segment_start

    if merged_difference_groups:
        focus_group = max(
            merged_difference_groups,
            key=lambda group: (group[2], group[1] - group[0])
        )
        focus_start, focus_end = focus_group[0], focus_group[1]
    else:
        focus_start, focus_end = full_visible_start, full_visible_end

    focus_width = max(focus_end - focus_start, 1)
    full_visible_width = max(full_visible_end - full_visible_start, 1)
    focus_padding = max(focus_width * 0.15, min(500, full_visible_width * 0.05))
    view_start = max(full_visible_start, focus_start - focus_padding)
    view_end = min(full_visible_end, focus_end + focus_padding)
    if view_end <= view_start:
        view_start, view_end = full_visible_start, full_visible_end

    overview_height = max(3.5, 1.2 + len(all_txs_to_plot) * 1.0)
    difference_height = max(4.5, 1.5 + len(all_txs_to_plot) * 1.25)
    combined_height = overview_height + difference_height + 1.5
    difference_figure = plt.figure(
        figsize=(20, combined_height), dpi=300, constrained_layout=False
    )
    combined_grid = GridSpec(
        2, 1, figure=difference_figure,
        height_ratios=[overview_height, difference_height], hspace=0.48
    )
    overview_axis = difference_figure.add_subplot(combined_grid[0, 0])
    difference_axis = difference_figure.add_subplot(combined_grid[1, 0])
    difference_figure.subplots_adjust(left=0.18, right=0.98, top=0.93, bottom=0.12)

    def draw_transcript_tracks(target_axis, x_start, x_end, label_fontsize, highlighted_segments=None):
        if highlighted_segments:
            for segment_start, segment_end in highlighted_segments:
                highlighted_start = max(segment_start, x_start)
                highlighted_end = min(segment_end, x_end)
                if highlighted_end > highlighted_start:
                    target_axis.axvspan(
                        highlighted_start, highlighted_end,
                        color='#F4A3A3', alpha=0.22, zorder=0
                    )

        target_y_positions = np.arange(len(all_txs_to_plot), 0, -1) * 2.0
        label_x = x_start - max((x_end - x_start) * 0.025, 20)
        for tx, y, label in zip(all_txs_to_plot, target_y_positions, track_labels):
            exons = exons_dict[tx]
            cds_list = cds_dict[tx]
            gene_start = max(map_pos(exons[0]['start']), x_start)
            gene_end = min(map_pos(exons[-1]['end']), x_end)
            if gene_end > gene_start:
                target_axis.plot(
                    [gene_start, gene_end], [y, y],
                    color='dimgray', linewidth=1.5, zorder=1
                )
            target_axis.text(
                label_x, y, label, ha='right', va='center',
                fontweight='bold', color='#222222',
                fontsize=label_fontsize, clip_on=False
            )

            for exon in exons:
                exon_start, exon_end = map_pos(exon['start']), map_pos(exon['end'])
                clipped_start = max(exon_start, x_start)
                clipped_end = min(exon_end, x_end)
                if clipped_end > clipped_start:
                    target_axis.add_patch(
                        patches.Rectangle(
                            (clipped_start, y - H_utr / 2),
                            clipped_end - clipped_start, H_utr,
                            linewidth=0.8, edgecolor='gray',
                            facecolor='#e0e0e0', zorder=3
                        )
                    )

            for cds in cds_list:
                cds_start, cds_end = map_pos(cds['start']), map_pos(cds['end'])
                clipped_start = max(cds_start, x_start)
                clipped_end = min(cds_end, x_end)
                if clipped_end > clipped_start:
                    target_axis.add_patch(
                        patches.Rectangle(
                            (clipped_start, y - H_cds / 2),
                            clipped_end - clipped_start, H_cds,
                            linewidth=0.8, edgecolor='#2C5F8A',
                            facecolor='#5B8DB8', zorder=4
                        )
                    )

            if not df_domains.empty and tx_info[tx]['protein_id']:
                tx_domains = df_domains[
                    df_domains['Protein_ID'].str.contains(
                        tx_info[tx]['protein_id'], na=False, regex=False
                    )
                ]
                for _, row_domain in tx_domains.iterrows():
                    domain_start = row_domain['Genomic_Start']
                    domain_end = row_domain['Genomic_End']
                    domain_color = dynamic_domain_colors.get(row_domain['Domain_Name'], 'orange')
                    for cds in cds_list:
                        overlap_start = max(cds['start'], domain_start)
                        overlap_end = min(cds['end'], domain_end)
                        if overlap_start >= overlap_end:
                            continue
                        visible_start = max(map_pos(overlap_start), x_start)
                        visible_end = min(map_pos(overlap_end), x_end)
                        if visible_end > visible_start:
                            target_axis.add_patch(
                                patches.Rectangle(
                                    (visible_start, y - H_domain / 2),
                                    visible_end - visible_start, H_domain,
                                    linewidth=0.35, edgecolor='#333333',
                                    facecolor=domain_color, alpha=0.9, zorder=5
                                )
                            )

        target_axis.set_xlim(x_start, x_end)
        target_axis.set_ylim(0.5, len(all_txs_to_plot) * 2.0 + 1.0)
        target_axis.axis('off')

    overview_axis.set_title(
        'A  Full transcript overview and zoom location',
        loc='left', fontsize=16, fontweight='bold', pad=8
    )
    draw_transcript_tracks(
        overview_axis, full_visible_start, full_visible_end, label_fontsize=10
    )
    overview_y_bottom = 0.62
    overview_y_top = len(all_txs_to_plot) * 2.0 + 0.88
    overview_axis.add_patch(
        patches.Rectangle(
            (view_start, overview_y_bottom),
            view_end - view_start, overview_y_top - overview_y_bottom,
            linewidth=2.2, edgecolor='#D62728', facecolor='none',
            linestyle='--', zorder=8
        )
    )
    overview_axis.text(
        (view_start + view_end) / 2, overview_y_top,
        'Zoomed region', color='#D62728', fontsize=10,
        fontweight='bold', ha='center', va='bottom'
    )

    difference_axis.set_title(
        'B  Difference-focused transcript structure',
        loc='left', fontsize=16, fontweight='bold', pad=8
    )
    draw_transcript_tracks(
        difference_axis, view_start, view_end, label_fontsize=12,
        highlighted_segments=difference_segments
    )

    difference_y_top = len(all_txs_to_plot) * 2.0 + 0.88
    for connector_x in (view_start, view_end):
        connector = patches.ConnectionPatch(
            xyA=(connector_x, overview_y_bottom),
            xyB=(connector_x, difference_y_top),
            coordsA=overview_axis.transData,
            coordsB=difference_axis.transData,
            color='#D62728', linewidth=1.2, linestyle='--',
            alpha=0.75, zorder=7, clip_on=False
        )
        difference_figure.add_artist(connector)

    overview_legend = [
        patches.Patch(facecolor='#e0e0e0', edgecolor='#999999', label='UTR'),
        patches.Patch(facecolor='#5B8DB8', edgecolor='#2C5F8A', label='CDS'),
        patches.Patch(facecolor='none', edgecolor='#D62728', linestyle='--', label='Zoomed region')
    ]
    for domain_name, domain_color in dynamic_domain_colors.items():
        overview_legend.append(
            patches.Patch(facecolor=domain_color, edgecolor='dimgray', label=domain_name)
        )
    overview_axis.legend(
        handles=overview_legend,
        loc='center left', bbox_to_anchor=(1.005, 0.5),
        frameon=False, prop={'size': 9, 'weight': 'bold'}
    )

    difference_legend = [
        patches.Patch(facecolor='#F4A3A3', edgecolor='none', alpha=0.35, label='Different CDS segment'),
        patches.Patch(facecolor='#e0e0e0', edgecolor='#999999', label='UTR'),
        patches.Patch(facecolor='#5B8DB8', edgecolor='#2C5F8A', label='CDS')
    ]
    focused_domain_names = set()
    for tx in all_txs_to_plot:
        if df_domains.empty or not tx_info[tx]['protein_id']:
            continue
        tx_domains = df_domains[
            df_domains['Protein_ID'].str.contains(
                tx_info[tx]['protein_id'], na=False, regex=False
            )
        ]
        for _, row_domain in tx_domains.iterrows():
            if map_pos(row_domain['Genomic_End']) > view_start and map_pos(row_domain['Genomic_Start']) < view_end:
                focused_domain_names.add(row_domain['Domain_Name'])
    for domain_name in sorted(focused_domain_names):
        domain_color = dynamic_domain_colors.get(domain_name, 'orange')
        difference_legend.append(
            patches.Patch(facecolor=domain_color, edgecolor='dimgray', label=domain_name)
        )
    difference_axis.legend(
        handles=difference_legend,
        loc='upper center', bbox_to_anchor=(0.5, -0.12),
        ncol=min(6, len(difference_legend)), frameon=False,
        prop={'size': 9, 'weight': 'bold'}
    )
    difference_filename = os.path.join(OUT_DIR, "IsoImpact_overview_zoom.png")
    difference_figure.savefig(
        difference_filename, dpi=300, bbox_inches='tight', facecolor='white'
    )
    plt.close(difference_figure)
    print(f"[IsoImpact] Overview-and-zoom figure written to {difference_filename}")

plt.close(fig)

print("[IsoImpact] Analysis complete.")
print(f"[IsoImpact] Figure written to {output_filename}")


