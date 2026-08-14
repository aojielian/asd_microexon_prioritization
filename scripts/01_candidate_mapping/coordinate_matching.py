import os
#!/usr/bin/env python3
"""Analysis: Coordinate-based matching of CHyMErA human events to Parikshak ASD universe.

Strategy:
1. Parse VastDB hg38 coordinates from existing mapping file
2. Convert hg38 -> hg19 using pyliftover
3. Parse Parikshak SE event coordinates (hg19) from event IDs
4. Match by gene + chromosome + coordinate proximity (±3bp tolerance)
5. Report all matches with delta_PSI, p-value, FDR
"""
import pandas as pd
import numpy as np
import re, sys
from pathlib import Path
from pyliftover import LiftOver

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", "."))
RESOURCE = PROJECT_ROOT / "09_resource_schema"
MAPPING = PROJECT_ROOT / "10_event_mapping"

SEED = 42
np.random.seed(SEED)

print("=" * 70)
print("Analysis: Coordinate-based Parikshak matching")
print("=" * 70)

# ============================================================
# 1. Load the mouse-human mapping (already done in 01_core_mapping.py)
# ============================================================
mapping = pd.read_csv(MAPPING / "05_mouse_human_event_mapping/01_mouse_human_all_conservation_links.tsv", sep='\t')
print(f"\n1. Loaded {len(mapping)} mouse-human conservation links")

# Parse hg38 coordinates from human_coord column (format: chr10:60082150-60082176)
coord_parsed = mapping['human_coord'].str.extract(r'(chr[\dXY]+):(\d+)-(\d+)')
mapping['vastdb_chr_hg38'] = coord_parsed[0]
mapping['vastdb_start_hg38'] = coord_parsed[1].astype(int)
mapping['vastdb_end_hg38'] = coord_parsed[2].astype(int)
mapping['exon_length_hg38'] = mapping['vastdb_end_hg38'] - mapping['vastdb_start_hg38']

print(f"   Parsed hg38 coordinates for {mapping['vastdb_chr_hg38'].notna().sum()}/36 events")
print(f"   Exon length range: {mapping['exon_length_hg38'].min()}-{mapping['exon_length_hg38'].max()} bp")

# ============================================================
# 2. Convert VastDB hg38 -> hg19 using pyliftover
# ============================================================
print("\n2. Converting VastDB hg38 coordinates to hg19 via pyliftover...")
lo_hg38_to_hg19 = LiftOver('hg38', 'hg19')

hg19_data = []
for _, row in mapping.iterrows():
    chrom = row['vastdb_chr_hg38']
    start = int(row['vastdb_start_hg38'])
    end = int(row['vastdb_end_hg38'])

    result_start = lo_hg38_to_hg19.convert_coordinate(chrom, start)
    result_end = lo_hg38_to_hg19.convert_coordinate(chrom, end)

    if result_start and result_end:
        hg19_chr = result_start[0][0]
        hg19_start = result_start[0][1]
        hg19_end = result_end[0][1]
        # Check strand flip
        strand_start = result_start[0][2]
        status = 'SUCCESS'
    else:
        hg19_chr = chrom
        hg19_start = None
        hg19_end = None
        strand_start = None
        status = 'UNMAPPED'

    hg19_data.append({
        'HsaEX_ID': row['HsaEX_ID'],
        'MmuEX_ID': row['MmuEX_ID'],
        'vastdb_chr_hg38': chrom,
        'vastdb_start_hg38': start,
        'vastdb_end_hg38': end,
        'lifted_chr_hg19': hg19_chr,
        'lifted_start_hg19': hg19_start,
        'lifted_end_hg19': hg19_end,
        'lifted_strand': strand_start,
        'liftover_status': status
    })

hg19_df = pd.DataFrame(hg19_data)
n_success = (hg19_df['liftover_status'] == 'SUCCESS').sum()
print(f"   LiftOver success: {n_success}/{len(hg19_df)}")
print(f"   LiftOver unmapped: {(hg19_df['liftover_status'] == 'UNMAPPED').sum()}")

# Save liftover results
hg19_df.to_csv(MAPPING / "04_parikshak_coordinate_harmonization/02_Parikshak_SE_hg38_liftover.tsv",
               sep='\t', index=False)

# ============================================================
# 3. Load Parikshak events and parse SE coordinates
# ============================================================
print("\n3. Loading Parikshak event universe...")
pk = pd.read_csv(RESOURCE / "06_parikshak/07_Parikshak_full_event_universe.tsv", sep='\t')
print(f"   Total Parikshak events: {len(pk)}")

# Parse SE exon coordinates from event IDs
def parse_se_coords(event_id):
    """Parse SE event coordinates from Parikshak event ID format.
    Format: chr10-str[-]-eS-61905725-eE-61905779-ueS-...-deS-...-deE-...
    """
    if not isinstance(event_id, str):
        return None, None, None, None
    es_match = re.search(r'eS-(\d+)-eE-(\d+)', event_id)
    chr_match = re.search(r'(chr[\dXY]+)-str\[([+-])\]', event_id)
    if es_match and chr_match:
        chrom = chr_match.group(1)
        strand = chr_match.group(2)
        es = int(es_match.group(1))
        ee = int(es_match.group(2))
        return chrom, strand, es, ee
    return None, None, None, None

def parse_a3ss_a5ss_coords(event_id):
    """Parse A3SS/A5SS coordinates from ENSG|chr:start-end-... format."""
    if not isinstance(event_id, str):
        return None, None, None, None
    match = re.search(r'\|(\d+|chr[\dXY]+):(\d+)-(\d+)', event_id)
    if match:
        chrom = match.group(1)
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        start = int(match.group(2))
        end = int(match.group(3))
        return chrom, None, start, end
    return None, None, None, None

# Parse all SE events
print("   Parsing SE event coordinates...")
se_events = pk[pk['event_type_standardized'] == 'SE'].copy()
se_coords_list = [parse_se_coords(eid) for eid in se_events['original_event_id']]
se_events['parsed_chr'] = [c[0] for c in se_coords_list]
se_events['parsed_strand'] = [c[1] for c in se_coords_list]
se_events['parsed_start'] = [c[2] for c in se_coords_list]
se_events['parsed_end'] = [c[3] for c in se_coords_list]
se_parsed = se_events[se_events['parsed_start'].notna()].copy()
print(f"   SE events with parsed coordinates: {len(se_parsed)}/{len(se_events)}")

# Parse A3SS/A5SS events
print("   Parsing A3SS/A5SS event coordinates...")
alt_events = pk[pk['event_type_standardized'].isin(['A3SS', 'A5SS'])].copy()
alt_coords_list = [parse_a3ss_a5ss_coords(eid) for eid in alt_events['original_event_id']]
alt_events['parsed_chr'] = [c[0] for c in alt_coords_list]
alt_events['parsed_strand'] = [c[1] for c in alt_coords_list]
alt_events['parsed_start'] = [c[2] for c in alt_coords_list]
alt_events['parsed_end'] = [c[3] for c in alt_coords_list]
alt_parsed = alt_events[alt_events['parsed_start'].notna()].copy()
print(f"   A3SS/A5SS events with parsed coordinates: {len(alt_parsed)}/{len(alt_events)}")

# ============================================================
# 4. Match CHyMErA human events to Parikshak by coordinate
# ============================================================
print("\n4. Matching CHyMErA human events to Parikshak...")

# Merge mapping with hg19 liftover coordinates
mapping_hg19 = mapping.merge(hg19_df[['HsaEX_ID', 'lifted_chr_hg19', 'lifted_start_hg19', 'lifted_end_hg19', 'lifted_strand', 'liftover_status']],
                              on='HsaEX_ID', how='left')

TOLERANCE_BP = 3  # Allow ±3bp for coordinate convention differences

all_matches = []
match_status = {}

for _, mrow in mapping_hg19.iterrows():
    mmuex = mrow['MmuEX_ID']
    hsaex = mrow['HsaEX_ID']
    gene = mrow['human_gene']
    lifted_chr = mrow['lifted_chr_hg19']
    lifted_start = mrow['lifted_start_hg19']
    lifted_end = mrow['lifted_end_hg19']
    vastdb_chr_hg38 = mrow['vastdb_chr_hg38']
    vastdb_start_hg38 = mrow['vastdb_start_hg38']
    vastdb_end_hg38 = mrow['vastdb_end_hg38']

    if pd.isna(lifted_start) or pd.isna(lifted_end):
        match_status[mmuex] = 'LIFTOVER_ERROR'
        continue

    lifted_start = int(lifted_start)
    lifted_end = int(lifted_end)

    # Search in SE events for this gene
    gene_se = se_parsed[se_parsed['gene_symbol_original'].str.upper() == gene.upper()]
    gene_alt = alt_parsed[alt_parsed['gene_symbol_original'].str.upper() == gene.upper()]

    found_match = False

    # Try SE matching (exact exon coordinates ±tolerance)
    for _, pk_row in gene_se.iterrows():
        pk_start = int(pk_row['parsed_start'])
        pk_end = int(pk_row['parsed_end'])
        pk_chr = pk_row['parsed_chr']

        if pk_chr != lifted_chr:
            continue

        dist_start = abs(pk_start - lifted_start)
        dist_end = abs(pk_end - lifted_end)

        if dist_start <= TOLERANCE_BP and dist_end <= TOLERANCE_BP:
            match_level = 'MATCH_B_EXACT_EXON_COORDINATE_AND_STRAND' if (dist_start <= 1 and dist_end <= 1) else f'MATCH_D_EXON_BOUNDARY_TOLERANCE_{max(dist_start,dist_end)}BP'
            all_matches.append({
                'MmuEX_ID': mmuex,
                'HsaEX_ID': hsaex,
                'gene': gene,
                'mouse_gene': mrow['mouse_gene'],
                'match_type': 'SE',
                'match_level': match_level,
                'Parikshak_event_id': pk_row['original_event_id'],
                'Parikshak_region': pk_row['region'],
                'Parikshak_event_type': pk_row['event_type_standardized'],
                'Parikshak_chr': pk_chr,
                'Parikshak_start_hg19': pk_start,
                'Parikshak_end_hg19': pk_end,
                'VastDB_chr_hg38': vastdb_chr_hg38,
                'VastDB_start_hg38': vastdb_start_hg38,
                'VastDB_end_hg38': vastdb_end_hg38,
                'lifted_start_hg19': lifted_start,
                'lifted_end_hg19': lifted_end,
                'coord_distance_start': dist_start,
                'coord_distance_end': dist_end,
                'strand': pk_row['parsed_strand'],
                'delta_psi': pk_row['delta_psi'],
                'p_value': pk_row['p_value'],
                'fdr': pk_row['fdr'],
                'mapping_method': 'VastDB_hg38_pyliftover_hg19_SE_coordinate_match',
                'exon_length_parikshak': pk_end - pk_start,
                'exon_length_vastdb': vastdb_end_hg38 - vastdb_start_hg38,
            })
            found_match = True

    # Try A3SS/A5SS matching (more tolerant - alternative splice sites)
    for _, pk_row in gene_alt.iterrows():
        pk_start = int(pk_row['parsed_start'])
        pk_end = int(pk_row['parsed_end'])
        pk_chr = pk_row['parsed_chr']

        if pk_chr != lifted_chr:
            continue

        dist_start = abs(pk_start - lifted_start)
        dist_end = abs(pk_end - lifted_end)

        # For A3SS/A5SS, one boundary should match (the shared exon boundary)
        # Use 5bp tolerance for the matching boundary
        boundary_match = (dist_start <= 5) or (dist_end <= 5)
        if boundary_match and max(dist_start, dist_end) <= 500:
            all_matches.append({
                'MmuEX_ID': mmuex,
                'HsaEX_ID': hsaex,
                'gene': gene,
                'mouse_gene': mrow['mouse_gene'],
                'match_type': 'A3SS_A5SS',
                'match_level': 'MATCH_E_CONSERVED_LOCAL_STRUCTURE',
                'Parikshak_event_id': pk_row['original_event_id'],
                'Parikshak_region': pk_row['region'],
                'Parikshak_event_type': pk_row['event_type_standardized'],
                'Parikshak_chr': pk_chr,
                'Parikshak_start_hg19': pk_start,
                'Parikshak_end_hg19': pk_end,
                'VastDB_chr_hg38': vastdb_chr_hg38,
                'VastDB_start_hg38': vastdb_start_hg38,
                'VastDB_end_hg38': vastdb_end_hg38,
                'lifted_start_hg19': lifted_start,
                'lifted_end_hg19': lifted_end,
                'coord_distance_start': dist_start,
                'coord_distance_end': dist_end,
                'strand': pk_row.get('strand_original', ''),
                'delta_psi': pk_row['delta_psi'],
                'p_value': pk_row['p_value'],
                'fdr': pk_row['fdr'],
                'mapping_method': 'VastDB_hg38_pyliftover_hg19_A3SS_A5SS_boundary_match',
                'exon_length_parikshak': pk_end - pk_start,
                'exon_length_vastdb': vastdb_end_hg38 - vastdb_start_hg38,
            })
            found_match = True

    if not found_match:
        # Check if gene exists in Parikshak at all
        gene_in_pk = pk['gene_symbol_original'].str.upper().eq(gene.upper()).any()
        match_status[mmuex] = f'NO_MATCH_{gene}_{"IN_PARIKSHAK" if gene_in_pk else "NOT_IN_PARIKSHAK"}'
    else:
        match_status[mmuex] = 'MATCHED'

# ============================================================
# 5. Results summary
# ============================================================
matches_df = pd.DataFrame(all_matches)
print(f"\n5. RESULTS:")
print(f"   Total match records: {len(matches_df)}")

if len(matches_df) > 0:
    n_matched_events = matches_df['MmuEX_ID'].nunique()
    print(f"   Unique CHyMErA events matched: {n_matched_events}/36")

    ctx = matches_df[matches_df['Parikshak_region'] == 'Cortex']
    cb = matches_df[matches_df['Parikshak_region'] == 'Cerebellum']
    print(f"   Cortex match records: {len(ctx)} ({ctx['MmuEX_ID'].nunique()} events)")
    print(f"   Cerebellum match records: {len(cb)} ({cb['MmuEX_ID'].nunique()} events)")

    # Strict matches (MATCH_B)
    strict = matches_df[matches_df['match_level'].str.startswith('MATCH_B')]
    print(f"   Strict (MATCH_B) records: {len(strict)} ({strict['MmuEX_ID'].nunique()} events)")

    # Match levels
    print(f"\n   Match level breakdown:")
    for level, count in matches_df['match_level'].value_counts().items():
        print(f"     {level}: {count}")

    # Show matched events
    print(f"\n   Matched events (deduplicated by event+region):")
    shown = matches_df.drop_duplicates(['MmuEX_ID', 'Parikshak_region'])
    for _, row in shown.sort_values('p_value').iterrows():
        sig = '**' if row['p_value'] < 0.01 else ('*' if row['p_value'] < 0.05 else '')
        print(f"     {row['gene']:10s} {row['MmuEX_ID']} | {row['Parikshak_region']:12s} | "
              f"dPSI={row['delta_psi']:+.4f} p={row['p_value']:.4f} FDR={row['fdr']:.3f} "
              f"| {row['match_level'][:30]} {sig}")

# Unmatched events
unmatched = [k for k, v in match_status.items() if v != 'MATCHED']
print(f"\n   Unmatched events ({len(unmatched)}):")
for ev in unmatched:
    gene = mapping[mapping['MmuEX_ID'] == ev]['human_gene'].values[0]
    print(f"     {ev} ({gene}): {match_status[ev]}")

# ============================================================
# 6. Save all outputs
# ============================================================
print("\n6. Saving outputs...")

# All matches
matches_df.to_csv(MAPPING / "06_asd_event_overlap/00_CHyMErA_human_Parikshak_all_matches.tsv",
                  sep='\t', index=False)

if len(matches_df) > 0:
    # Strict matches
    strict = matches_df[matches_df['match_level'].str.startswith('MATCH_B')]
    ctx_strict = strict[strict['Parikshak_region'] == 'Cortex']
    cb_strict = strict[strict['Parikshak_region'] == 'Cerebellum']
    ctx_strict.to_csv(MAPPING / "06_asd_event_overlap/01_CHyMErA_CTX_strict_event_overlap.tsv", sep='\t', index=False)
    cb_strict.to_csv(MAPPING / "06_asd_event_overlap/02_CHyMErA_CB_strict_event_overlap.tsv", sep='\t', index=False)

    # Tolerant matches
    ctx_all = matches_df[matches_df['Parikshak_region'] == 'Cortex']
    cb_all = matches_df[matches_df['Parikshak_region'] == 'Cerebellum']
    ctx_all.to_csv(MAPPING / "06_asd_event_overlap/03_CHyMErA_CTX_tolerant_overlap.tsv", sep='\t', index=False)
    cb_all.to_csv(MAPPING / "06_asd_event_overlap/04_CHyMErA_CB_tolerant_overlap.tsv", sep='\t', index=False)

    # Multiplicity check
    mult = matches_df.groupby(['MmuEX_ID', 'HsaEX_ID', 'gene']).agg(
        n_parikshak_records=('Parikshak_event_id', 'count'),
        n_regions=('Parikshak_region', 'nunique'),
        regions=('Parikshak_region', lambda x: ','.join(sorted(x.unique()))),
        match_types=('match_type', lambda x: ','.join(sorted(x.unique())))
    ).reset_index()
    mult.to_csv(MAPPING / "06_asd_event_overlap/05_event_overlap_multiplicity_check.tsv", sep='\t', index=False)
else:
    for f in ['01_CHyMErA_CTX_strict_event_overlap.tsv', '02_CHyMErA_CB_strict_event_overlap.tsv',
              '03_CHyMErA_CTX_tolerant_overlap.tsv', '04_CHyMErA_CB_tolerant_overlap.tsv',
              '05_event_overlap_multiplicity_check.tsv']:
        pd.DataFrame().to_csv(MAPPING / "06_asd_event_overlap" / f, sep='\t', index=False)

# Errors
errors = []
for ev in unmatched:
    gene = mapping[mapping['MmuEX_ID'] == ev]['human_gene'].values[0]
    hsaex = mapping[mapping['MmuEX_ID'] == ev]['HsaEX_ID'].values[0]
    lo_status = hg19_df[hg19_df['HsaEX_ID'] == hsaex]['liftover_status'].values
    errors.append({
        'MmuEX_ID': ev,
        'HsaEX_ID': hsaex,
        'gene': gene,
        'error_reason': match_status[ev],
        'liftover_status': lo_status[0] if len(lo_status) > 0 else 'NOT_FOUND',
        'gene_in_parikshak': 'YES' if 'IN_PARIKSHAK' in match_status[ev] else 'NO'
    })
pd.DataFrame(errors).to_csv(MAPPING / "06_asd_event_overlap/06_event_overlap_errors.tsv", sep='\t', index=False)

# Summary
n_ctx_matched = ctx['MmuEX_ID'].nunique() if len(matches_df) > 0 else 0
n_cb_matched = cb['MmuEX_ID'].nunique() if len(matches_df) > 0 else 0
n_both = len(set(ctx['MmuEX_ID'].unique()) & set(cb['MmuEX_ID'].unique())) if len(matches_df) > 0 else 0

summary_data = [
    {'metric': 'N_CHYMERA_STRICT_HUMAN_EVENTS', 'value': 36},
    {'metric': 'N_LIFTOVER_SUCCESS', 'value': n_success},
    {'metric': 'N_CTX_DIRECT_MATCHED_EVENTS', 'value': n_ctx_matched},
    {'metric': 'N_CB_DIRECT_MATCHED_EVENTS', 'value': n_cb_matched},
    {'metric': 'N_MATCHED_IN_BOTH_REGIONS', 'value': n_both},
    {'metric': 'N_CTX_RECORDS', 'value': len(ctx) if len(matches_df) > 0 else 0},
    {'metric': 'N_CB_RECORDS', 'value': len(cb) if len(matches_df) > 0 else 0},
    {'metric': 'N_CTX_NOMINAL_P_LT_0_05', 'value': int((ctx['p_value'] < 0.05).sum()) if len(ctx) > 0 else 0},
    {'metric': 'N_CTX_NOMINAL_P_LT_0_01', 'value': int((ctx['p_value'] < 0.01).sum()) if len(ctx) > 0 else 0},
    {'metric': 'N_CTX_FDR_LT_0_05', 'value': int((ctx['fdr'] < 0.05).sum()) if len(ctx) > 0 else 0},
    {'metric': 'N_CB_NOMINAL_P_LT_0_05', 'value': int((cb['p_value'] < 0.05).sum()) if len(cb) > 0 else 0},
    {'metric': 'N_CB_FDR_LT_0_05', 'value': int((cb['fdr'] < 0.05).sum()) if len(cb) > 0 else 0},
    {'metric': 'COORDINATE_TOLERANCE_BP', 'value': TOLERANCE_BP},
    {'metric': 'MATCHING_METHOD', 'value': 'VastDB_hg38_pyliftover_hg19_coordinate_proximity'},
    {'metric': 'PARIKSHAK_SE_PARSED', 'value': len(se_parsed)},
    {'metric': 'PARIKSHAK_A3SS_A5SS_PARSED', 'value': len(alt_parsed)},
]
pd.DataFrame(summary_data).to_csv(MAPPING / "06_asd_event_overlap/07_event_overlap_summary.tsv", sep='\t', index=False)

# Phase status
n_strict_ctx = strict['MmuEX_ID'].nunique() if len(matches_df) > 0 and len(strict) > 0 else 0
if n_strict_ctx >= 5:
    phase = 'OK'
elif n_ctx_matched >= 3:
    phase = 'OK_TOLERANT'
else:
    phase = 'HOLD'

check_df = pd.DataFrame([
    {'check_item': 'PARIKSHAK_OVERLAP_STATUS', 'status': phase,
     'evidence': f'{n_ctx_matched} CTX events matched ({n_strict_ctx} strict), {n_cb_matched} CB events matched'}
])
check_df.to_csv(MAPPING / "06_asd_event_overlap/08_event_overlap_check.tsv", sep='\t', index=False)

print(f"\n   Phase status: {phase}")
print(f"   All outputs saved to {MAPPING / '06_asd_event_overlap/'}")
print("\n" + "=" * 70)
print("COORDINATE MATCHING COMPLETE")
print("=" * 70)
