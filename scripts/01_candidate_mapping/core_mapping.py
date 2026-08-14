import os
#!/usr/bin/env python3
"""Analysis Core: Mouse->Human event mapping via VastDB EVENT_CONSERVATION, then Parikshak overlap."""
import gzip, csv, os, sys, hashlib
from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime, timezone

import numpy as np
np.random.seed(42)

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", "."))
MAPPING = PROJECT_ROOT / "10_event_mapping"
RESOURCE = PROJECT_ROOT / "09_resource_schema"
VASTDB = PROJECT_ROOT / "05_vastdb"
RANDOM_SEED = 42

print("=" * 70)
print("ANALYSIS: EVENT-LEVEL CROSS-SPECIES MAPPING & ASD FUNCTIONAL BRIDGE")
print(f"Started: {datetime.now(timezone.utc).isoformat()}")
print(f"Random seed: {RANDOM_SEED}")
print("=" * 70)

# ============================================================
# STEP 1: Load CHyMErA target events
# ============================================================
print("\n[STEP 1] Loading CHyMErA target microexon events...")
chymera_master = RESOURCE / "03_chymera/01_CHyMErA_target_microexon_master.tsv"
chymera_events = []
with open(chymera_master) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['target_class'] == 'MICROEXON_DELETION':
            chymera_events.append(row)

mmuex_ids = [r['chymera_target_id'] for r in chymera_events]
print(f"  CHyMErA microexon targeting events: {len(mmuex_ids)}")
print(f"  MmuEX IDs: {sorted(mmuex_ids)}")

# ============================================================
# STEP 2: Get mouse event coordinates from VastDB mm10 EVENT_INFO
# ============================================================
print("\n[STEP 2] Looking up CHyMErA events in VastDB mm10 EVENT_INFO...")
mm10_event_info = {}
with gzip.open(VASTDB / "mm10/EVENT_INFO-mm10.tab.gz", 'rt', encoding='latin-1') as f:
    header = f.readline().strip().split('\t')
    for line in f:
        parts = line.strip().split('\t')
        event_id = parts[1]  # EVENT column
        if event_id in mmuex_ids:
            mm10_event_info[event_id] = {
                'GENE': parts[0],
                'EVENT': parts[1],
                'COORD_o': parts[2],
                'LE_o': parts[3],
                'FULL_CO': parts[4],
                'COMPLEX': parts[5],
                'REF_CO': parts[6],
            }

print(f"  Found {len(mm10_event_info)}/{len(mmuex_ids)} CHyMErA events in VastDB mm10")
for eid in sorted(mmuex_ids):
    if eid in mm10_event_info:
        info = mm10_event_info[eid]
        print(f"    {eid}: {info['GENE']}, {info['COORD_o']}, type={info['COMPLEX']}")
    else:
        print(f"    {eid}: NOT FOUND IN VastDB mm10!")

# ============================================================
# STEP 3: Query EVENT_CONSERVATION for mouse->human mappings
# ============================================================
print("\n[STEP 3] Querying EVENT_CONSERVATION for mouse->human mappings...")
# EVENT_CONSERVATION: EventID, Ass1, Ass2, Chro, Start, End, Type, ConservedID
# We need to find rows where EventID is one of our MmuEX and Ass2=hg38
# OR where ConservedID is one of our MmuEX (reverse direction)

conservation_hits = defaultdict(list)  # MmuEX -> list of HsaEX
cons_path = VASTDB / "general/EVENT_CONSERVATION.tab.gz"

print("  Scanning EVENT_CONSERVATION (21M rows)...")
with gzip.open(cons_path, 'rt', encoding='latin-1') as f:
    f.readline()  # header
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 8:
            continue
        event_id = parts[0]
        ass1 = parts[1]
        ass2 = parts[2]
        conserved_id = parts[7]

        # Direction 1: MmuEX is source, hg38 is target
        if event_id in mmuex_ids and ass2 == 'hg38':
            conservation_hits[event_id].append({
                'HsaEX': conserved_id,
                'chr': parts[3],
                'start': parts[4],
                'end': parts[5],
                'type': parts[6],
                'direction': 'MmuEX_is_source',
            })
        # Direction 2: HsaEX is source, MmuEX is conserved target
        elif conserved_id in mmuex_ids and ass1 == 'hg38':
            conservation_hits[conserved_id].append({
                'HsaEX': event_id,
                'chr': parts[3],
                'start': parts[4],
                'end': parts[5],
                'type': parts[6],
                'direction': 'MmuEX_is_target',
            })

print(f"\n  Conservation results:")
n_with_human = 0
total_links = 0
for mmuex in sorted(mmuex_ids):
    hits = conservation_hits.get(mmuex, [])
    # Deduplicate by HsaEX
    unique_hsa = set(h['HsaEX'] for h in hits)
    if unique_hsa:
        n_with_human += 1
        total_links += len(unique_hsa)
        print(f"    {mmuex} -> {len(unique_hsa)} human event(s): {sorted(unique_hsa)[:5]}")
    else:
        print(f"    {mmuex} -> NO HUMAN CONSERVATION")

print(f"\n  Summary: {n_with_human}/{len(mmuex_ids)} CHyMErA events have human orthologs")
print(f"  Total unique links: {total_links}")

# ============================================================
# STEP 4: Get human event info from VastDB hg38 EVENT_INFO
# ============================================================
print("\n[STEP 4] Looking up human events in VastDB hg38 EVENT_INFO...")
all_hsa_ids = set()
for hits in conservation_hits.values():
    for h in hits:
        all_hsa_ids.add(h['HsaEX'])

print(f"  Need to look up {len(all_hsa_ids)} unique HsaEX IDs in hg38 EVENT_INFO")

hg38_event_info = {}
with gzip.open(VASTDB / "hg38/EVENT_INFO-hg38.tab.gz", 'rt', encoding='latin-1') as f:
    f.readline()  # header
    for line in f:
        parts = line.strip().split('\t')
        event_id = parts[1]
        if event_id in all_hsa_ids:
            hg38_event_info[event_id] = {
                'GENE': parts[0],
                'EVENT': parts[1],
                'COORD_o': parts[2],
                'LE_o': parts[3],
                'FULL_CO': parts[4],
                'COMPLEX': parts[5],
                'REF_CO': parts[6],
            }

print(f"  Found {len(hg38_event_info)}/{len(all_hsa_ids)} human events in VastDB hg38")

# ============================================================
# STEP 5: Load Parikshak full event universe
# ============================================================
print("\n[STEP 5] Loading Parikshak event universe...")
parikshak_path = RESOURCE / "06_parikshak/07_Parikshak_full_event_universe.tsv"
parikshak_events = []
with open(parikshak_path) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        parikshak_events.append(row)
print(f"  Loaded {len(parikshak_events)} Parikshak events")

# Build lookup by event ID and by gene
parikshak_by_id = defaultdict(list)
parikshak_by_gene = defaultdict(list)
for pe in parikshak_events:
    parikshak_by_id[pe['original_event_id']].append(pe)
    if pe['gene_symbol_original']:
        parikshak_by_gene[pe['gene_symbol_original']].append(pe)

# Also check if any Parikshak A3SS/A5SS events contain HsaEX IDs
# Format: ENSG...|chr:coords|HsaEX...
parikshak_hsa_lookup = {}
for pe in parikshak_events:
    eid = pe['original_event_id']
    if 'Hsa' in eid:
        # Extract HsaEX ID
        parts = eid.split('|')
        for p in parts:
            if p.startswith('Hsa'):
                parikshak_hsa_lookup[p] = pe
                break

print(f"  Parikshak events with embedded HsaEX IDs: {len(parikshak_hsa_lookup)}")

# ============================================================
# STEP 6: Match CHyMErA human events to Parikshak
# ============================================================
print("\n[STEP 6] Matching CHyMErA human orthologs to Parikshak events...")

matches = []
for mmuex in sorted(mmuex_ids):
    hits = conservation_hits.get(mmuex, [])
    unique_hsa = {}
    for h in hits:
        hsa = h['HsaEX']
        if hsa not in unique_hsa:
            unique_hsa[hsa] = h

    mouse_info = mm10_event_info.get(mmuex, {})
    mouse_gene = mouse_info.get('GENE', '')

    for hsa_id, hit_info in unique_hsa.items():
        human_info = hg38_event_info.get(hsa_id, {})
        human_gene = human_info.get('GENE', '')

        # Check if this HsaEX is directly in Parikshak
        match_method = 'NO_MATCH'
        parikshak_record = None

        # Method A: Direct HsaEX ID match
        if hsa_id in parikshak_hsa_lookup:
            match_method = 'MATCH_A_EXACT_HsaEX_ID'
            parikshak_record = parikshak_hsa_lookup[hsa_id]
        else:
            # Method B: Search by gene + coordinate overlap
            # Get human exon coordinates
            human_coord = human_info.get('COORD_o', '')
            if human_gene:
                gene_events = parikshak_by_gene.get(human_gene, [])
                for pe in gene_events:
                    # Check coordinate overlap
                    pe_id = pe['original_event_id']
                    if hsa_id in pe_id:
                        match_method = 'MATCH_A_EXACT_HsaEX_ID'
                        parikshak_record = pe
                        break
                    # For SE events, check if coordinates overlap
                    if pe['event_type_original'] == 'SE' and human_coord:
                        # Parse human coord: chr12:120210930-120211106
                        try:
                            h_parts = human_coord.split(':')
                            h_chr = h_parts[0]
                            h_range = h_parts[1]
                            h_start, h_end = h_range.split('-')
                            h_start, h_end = int(h_start), int(h_end)

                            # Parse Parikshak SE coord: chr12-str[-]-eS-123425353-eE-123425530
                            if 'eS-' in pe_id and 'eE-' in pe_id:
                                p_parts = pe_id.split('-')
                                for i, p in enumerate(p_parts):
                                    if p == 'eS' and i+1 < len(p_parts):
                                        p_start = int(p_parts[i+1])
                                    if p == 'eE' and i+1 < len(p_parts):
                                        p_end = int(p_parts[i+1])
                                p_chr = p_parts[0]
                                if p_chr == h_chr and abs(p_start - h_start) <= 3 and abs(p_end - h_end) <= 3:
                                    match_method = 'MATCH_C_EXON_BOUNDARY_TOLERANCE_3BP'
                                    parikshak_record = pe
                                    break
                        except:
                            ok

        matches.append({
            'MmuEX_ID': mmuex,
            'mouse_gene': mouse_gene,
            'mouse_coord': mouse_info.get('COORD_o', ''),
            'mouse_type': mouse_info.get('COMPLEX', ''),
            'HsaEX_ID': hsa_id,
            'human_gene': human_gene,
            'human_coord': human_info.get('COORD_o', ''),
            'human_type': human_info.get('COMPLEX', ''),
            'conservation_type': hit_info['type'],
            'match_to_parikshak': match_method,
            'parikshak_region': parikshak_record['region'] if parikshak_record else '',
            'parikshak_delta_psi': parikshak_record['delta_psi'] if parikshak_record else '',
            'parikshak_p_value': parikshak_record['p_value'] if parikshak_record else '',
            'parikshak_fdr': parikshak_record['fdr'] if parikshak_record else '',
            'parikshak_event_id': parikshak_record['original_event_id'][:60] if parikshak_record else '',
        })

print(f"\n  Total mapping links: {len(matches)}")
match_counts = Counter(m['match_to_parikshak'] for m in matches)
for mc, cnt in match_counts.most_common():
    print(f"    {mc}: {cnt}")

# ============================================================
# STEP 7: Write output files
# ============================================================
print("\n[STEP 7] Writing output files...")

# 05_mouse_human_event_mapping/01_mouse_human_all_conservation_links.tsv
out_dir = MAPPING / "05_mouse_human_event_mapping"
with open(out_dir / "01_mouse_human_all_conservation_links.tsv", 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=matches[0].keys(), delimiter='\t')
    w.writeheader()
    w.writerows(matches)

# 02_mouse_human_event_map_strict.tsv (LEVEL_1 only: direct conservation)
strict = [m for m in matches if m['conservation_type'] == 'ALTERNATIVE']
with open(out_dir / "02_mouse_human_event_map_strict.tsv", 'w', newline='') as f:
    if strict:
        w = csv.DictWriter(f, fieldnames=strict[0].keys(), delimiter='\t')
        w.writeheader()
        w.writerows(strict)
    else:
        f.write('\t'.join(matches[0].keys()) + '\n')

# Summary stats
n_level1 = sum(1 for m in matches if m['conservation_type'] == 'ALTERNATIVE')
n_total_links = len(matches)
n_mmuex_with_human = len(set(m['MmuEX_ID'] for m in matches))
n_mmuex_with_parikshak = len(set(m['MmuEX_ID'] for m in matches if m['match_to_parikshak'] != 'NO_MATCH'))
n_parikshak_matched = sum(1 for m in matches if m['match_to_parikshak'] != 'NO_MATCH')

# 06_mouse_human_mapping_QC.tsv
with open(out_dir / "06_mouse_human_mapping_QC.tsv", 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['metric', 'value'])
    w.writerow(['N_CHYMERA_EVENTS', len(mmuex_ids)])
    w.writerow(['N_IN_VASTDB_MM10', len(mm10_event_info)])
    w.writerow(['N_WITH_ANY_HUMAN_EVENT', n_mmuex_with_human])
    w.writerow(['N_LEVEL1_ALTERNATIVE', n_level1])
    w.writerow(['N_TOTAL_CONSERVATION_LINKS', n_total_links])
    w.writerow(['N_ONE_TO_ONE', sum(1 for mmuex in mmuex_ids if len(set(h['HsaEX'] for h in conservation_hits.get(mmuex, []))) == 1)])
    w.writerow(['N_ONE_TO_MANY', sum(1 for mmuex in mmuex_ids if len(set(h['HsaEX'] for h in conservation_hits.get(mmuex, []))) > 1)])
    w.writerow(['N_UNMAPPED', len(mmuex_ids) - n_mmuex_with_human])
    w.writerow(['N_WITH_PARIKSHAK_MATCH', n_mmuex_with_parikshak])
    w.writerow(['N_PARIKSHAK_MATCHED_RECORDS', n_parikshak_matched])

# 07_mouse_human_mapping_check.tsv
with open(out_dir / "07_mouse_human_mapping_check.tsv", 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['check', 'status', 'evidence'])
    w.writerow(['chymera_events_in_vastdb', 'OK' if len(mm10_event_info) == len(mmuex_ids) else 'PARTIAL', f'{len(mm10_event_info)}/{len(mmuex_ids)}'])
    w.writerow(['conservation_mapping_available', 'OK' if n_mmuex_with_human > 0 else 'ERROR', f'{n_mmuex_with_human}/{len(mmuex_ids)} have human orthologs'])
    w.writerow(['parikshak_overlap', 'OK' if n_parikshak_matched > 0 else 'ERROR', f'{n_parikshak_matched} records matched'])
    w.writerow(['mapping_check', 'OK' if n_mmuex_with_human >= 10 else 'HOLD', f'{n_mmuex_with_human} events mappable'])

# ============================================================
# STEP 8: Build ASD evidence table for matched events
# ============================================================
print("\n[STEP 8] Building ASD evidence summary for matched events...")
asd_evidence = [m for m in matches if m['match_to_parikshak'] != 'NO_MATCH']
print(f"  Events with Parikshak ASD evidence: {len(asd_evidence)}")
for m in asd_evidence:
    print(f"    {m['MmuEX_ID']} ({m['mouse_gene']}) -> {m['HsaEX_ID']} ({m['human_gene']}): "
          f"ΔPSI={m['parikshak_delta_psi']}, p={m['parikshak_p_value']}, FDR={m['parikshak_fdr']} [{m['parikshak_region']}]")

# 06_asd_event_overlap
overlap_dir = MAPPING / "06_asd_event_overlap"
with open(overlap_dir / "00_CHyMErA_human_Parikshak_all_matches.tsv", 'w', newline='') as f:
    if asd_evidence:
        w = csv.DictWriter(f, fieldnames=asd_evidence[0].keys(), delimiter='\t')
        w.writeheader()
        w.writerows(asd_evidence)
    else:
        f.write('\t'.join(matches[0].keys()) + '\n')

# Summary
ctx_matches = [m for m in asd_evidence if m['parikshak_region'] == 'Cortex']
cb_matches = [m for m in asd_evidence if m['parikshak_region'] == 'Cerebellum']

with open(overlap_dir / "07_event_overlap_summary.tsv", 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['metric', 'value'])
    w.writerow(['N_CHYMERA_STRICT_HUMAN_EVENTS', n_mmuex_with_human])
    w.writerow(['N_CTX_DIRECT_MATCHED_EVENTS', len(set(m['MmuEX_ID'] for m in ctx_matches))])
    w.writerow(['N_CB_DIRECT_MATCHED_EVENTS', len(set(m['MmuEX_ID'] for m in cb_matches))])
    w.writerow(['N_CTX_RECORDS', len(ctx_matches)])
    w.writerow(['N_CB_RECORDS', len(cb_matches)])
    # Count significant
    ctx_sig = [m for m in ctx_matches if m['parikshak_p_value'] and m['parikshak_p_value'] != '']
    ctx_p05 = 0
    ctx_p01 = 0
    ctx_fdr05 = 0
    for m in ctx_matches:
        try:
            p = float(m['parikshak_p_value'])
            if p < 0.05: ctx_p05 += 1
            if p < 0.01: ctx_p01 += 1
        except:
            ok
        try:
            fdr = float(m['parikshak_fdr'])
            if fdr < 0.05: ctx_fdr05 += 1
        except:
            ok
    w.writerow(['N_CTX_NOMINAL_P_LT_0_05', ctx_p05])
    w.writerow(['N_CTX_NOMINAL_P_LT_0_01', ctx_p01])
    w.writerow(['N_CTX_FDR_LT_0_05', ctx_fdr05])

# ============================================================
# STEP 9: Background universe and set-level test
# ============================================================
print("\n[STEP 9] Building background and running set-level tests...")

# Background: all Parikshak events that have VastDB HsaEX IDs (from A3SS/A5SS format)
# and are NOT in our target set
target_hsa_ids = set(m['HsaEX_ID'] for m in matches)
background_events = []
for hsa_id, pe in parikshak_hsa_lookup.items():
    if hsa_id not in target_hsa_ids:
        background_events.append({'HsaEX_ID': hsa_id, **pe})

print(f"  Background events with HsaEX IDs: {len(background_events)}")

# For set-level test: compare |deltaPSI| of target vs background
target_dpsi = []
for m in ctx_matches:
    try:
        target_dpsi.append(abs(float(m['parikshak_delta_psi'])))
    except:
        ok

bg_dpsi_ctx = []
for bg in background_events:
    if bg.get('region') == 'Cortex':
        try:
            bg_dpsi_ctx.append(abs(float(bg['delta_psi'])))
        except:
            ok

print(f"  Target |ΔPSI| values (CTX): n={len(target_dpsi)}")
print(f"  Background |ΔPSI| values (CTX): n={len(bg_dpsi_ctx)}")

if target_dpsi and bg_dpsi_ctx:
    from scipy import stats
    # Wilcoxon rank-sum test
    stat, p_val = stats.mannwhitneyu(target_dpsi, bg_dpsi_ctx, alternative='greater')
    target_median = np.median(target_dpsi)
    bg_median = np.median(bg_dpsi_ctx)

    # Permutation test
    n_perm = 10000
    combined = np.array(target_dpsi + bg_dpsi_ctx)
    n_target = len(target_dpsi)
    observed_diff = np.mean(target_dpsi) - np.mean(bg_dpsi_ctx)
    perm_diffs = []
    rng = np.random.default_rng(RANDOM_SEED)
    for _ in range(n_perm):
        perm = rng.permutation(combined)
        perm_diff = perm[:n_target].mean() - perm[n_target:].mean()
        perm_diffs.append(perm_diff)
    perm_p = (sum(1 for d in perm_diffs if d >= observed_diff) + 1) / (n_perm + 1)

    print(f"\n  === SET-LEVEL TEST RESULTS (CTX) ===")
    print(f"  Target median |ΔPSI|: {target_median:.4f}")
    print(f"  Background median |ΔPSI|: {bg_median:.4f}")
    print(f"  Mean difference: {observed_diff:.4f}")
    print(f"  Mann-Whitney U p-value (greater): {p_val:.4f}")
    print(f"  Permutation p-value ({n_perm} perms): {perm_p:.4f}")
else:
    p_val = 'N/A'
    perm_p = 'N/A'
    observed_diff = 'N/A'
    target_median = 'N/A'
    bg_median = 'N/A'

# Write set-level results
set_dir = MAPPING / "08_set_level_tests"
with open(set_dir / "01_CTX_set_level_continuous_tests.tsv", 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['test', 'statistic', 'effect_size', 'p_value', 'n_target', 'n_background', 'notes'])
    w.writerow(['mannwhitney_abs_delta_psi', f'{stat if target_dpsi and bg_dpsi_ctx else "N/A"}', f'{observed_diff}', f'{p_val}', len(target_dpsi), len(bg_dpsi_ctx), 'one-sided greater'])
    w.writerow(['permutation_abs_delta_psi', f'{observed_diff}', f'{observed_diff}', f'{perm_p}', len(target_dpsi), len(bg_dpsi_ctx), f'{RANDOM_SEED} seed, {n_perm if target_dpsi else 0} permutations'])

with open(set_dir / "10_set_level_effect_summary.tsv", 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['metric', 'value'])
    w.writerow(['target_median_abs_delta_psi', target_median])
    w.writerow(['background_median_abs_delta_psi', bg_median])
    w.writerow(['mean_difference', observed_diff])
    w.writerow(['mannwhitney_p', p_val])
    w.writerow(['permutation_p', perm_p])
    w.writerow(['n_target_matched_CTX', len(ctx_matches)])
    w.writerow(['n_background_CTX', len(bg_dpsi_ctx)])

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS CORE MAPPING SUMMARY")
print("=" * 70)
print(f"  CHyMErA events: {len(mmuex_ids)}")
print(f"  In VastDB mm10: {len(mm10_event_info)}")
print(f"  With human conservation: {n_mmuex_with_human}")
print(f"  Matched to Parikshak: {n_mmuex_with_parikshak}")
print(f"  CTX matches: {len(ctx_matches)} records")
print(f"  CB matches: {len(cb_matches)} records")
print(f"  Background eligible: {len(background_events)}")
if target_dpsi and bg_dpsi_ctx:
    print(f"  Set-level p (permutation): {perm_p}")
print("=" * 70)
