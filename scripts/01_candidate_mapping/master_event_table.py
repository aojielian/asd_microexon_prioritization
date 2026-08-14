#!/usr/bin/env python3
"""Phase D: 19-event master evidence table (TSV+XLSX) + logic check.
Built from Final models authoritative master table + Primary mixed models (adds SE, 95% CI).
"""
import os, csv
ROOT=os.environ.get("PROJECT_ROOT", ".")
OUT=os.path.join(ROOT,"25_master_evidence")
MET=os.path.join(OUT,"06_master_event_table")
os.makedirs(MET,exist_ok=True)

def rd(rel):
    with open(os.path.join(ROOT,rel)) as f: return list(csv.DictReader(f,delimiter="\t"))

MT = rd("24_event_annotation_finalization/03_master_table_logic_repair/MASTER_19_EVENT_EVIDENCE_TABLE_FINAL.tsv")
MODELS = {r["HsaEX_ID"]:r for r in rd("21_coordinate_inference/05_mixed_model_inference/01_Satterthwaite_models.tsv")}

cols=["HsaEX_ID","MmuEX_ID","gene",
 "source_chr_hg19","source_start_hg19","source_end_hg19","chr_hg38","start_hg38","end_hg38",
 "liftover_status","roundtrip_status","GENCODE_v33_local_structure_status",
 "CHyMErA_perturbation_status","CHyMErA_direction","CHyMErA_direction_concordant","CHyMErA_bridge_classification","CHyMErA_reason_direction_not_determined",
 "Parikshak_delta_PSI","Parikshak_P","Parikshak_FDR","Parikshak_direction",
 "developmental_dynamic_status","developmental_trajectory","developmental_timing_tier","PSI_range","prenatal_mean_PSI","postnatal_mean_PSI","monotonicity_rho","monotonicity_p",
 "network_module_or_pathway",
 "GSE30573_mapping_status","GSE30573_direction","GSE30573_direction_concordant","GSE30573_support_level",
 "PsychENCODE_beta","PsychENCODE_SE","PsychENCODE_CI95_lower","PsychENCODE_CI95_upper","PsychENCODE_direction","direction_concordant",
 "P_KR","BH_FDR_KR","P_LRT","BH_FDR_LRT","PsychENCODE_primary_significance","PsychENCODE_sensitivity_significance",
 "inclusion_transcript_count","exclusion_transcript_count","transcript_usage_definition",
 "final_evidence_tier","positive_evidence_summary","negative_evidence_summary","manuscript_claim_level"]

# Publication-facing vocabulary applied at the render layer only
# (upstream prespecified enums are internal; numeric values are never altered).
VOCAB = {
    "NOT_DIRECTION_TESTABLE": "NOT_DIRECTION_ELIGIBLE",
    "DIRECTION_TESTABLE": "DIRECTION_ELIGIBLE",
    "UNMAPPED_NOT_TESTABLE": "UNMAPPED_NOT_ELIGIBLE",
    "NOT_EVALUABLE_NO_ASD_DIRECTION": "NO_ASD_EFFECT_DIRECTION",
    "ASD_DIRECTION_UNRESOLVED_P_GT_0.05": "ASD_DIRECTION_UNRESOLVED",
}

def vocab(v):
    return VOCAB.get(v, v)

rows=[]
for r in MT:
    m=MODELS[r["HsaEX_ID"]]
    beta=float(m["beta_ASD"]); se=float(m["SE"])
    nr=dict(r)
    nr["PsychENCODE_SE"]="%.6f"%se
    nr["PsychENCODE_CI95_lower"]="%.6f"%(beta-1.96*se)
    nr["PsychENCODE_CI95_upper"]="%.6f"%(beta+1.96*se)
    nr["PsychENCODE_beta"]="%.6f"%beta
    rows.append(nr)

with open(os.path.join(MET,"MASTER_19_EVENT_EVIDENCE_TABLE.tsv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter="\t",extrasaction="ignore")
    w.writeheader()
    for nr in rows:
        w.writerow({k: vocab(v) for k, v in nr.items()})

# XLSX
try:
    import openpyxl
    wb=openpyxl.Workbook(); ws=wb.active; ws.title="MASTER_19_EVENTS"
    ws.append(cols)
    for nr in rows:
        ws.append([vocab(nr.get(c, "")) for c in cols])
    wb.save(os.path.join(MET,"MASTER_19_EVENT_EVIDENCE_TABLE.xlsx"))
    xlsx_ok="YES"
except Exception as e:
    xlsx_ok="NO:"+str(e)

# LOGIC CHECK (4 rules)
check=[]
def allowed_claim_ok(nr):
    tier=nr["final_evidence_tier"]; claim=nr["manuscript_claim_level"]
    disc_pe = nr["direction_concordant"]=="FALSE"
    fdrkr=float(nr["BH_FDR_KR"])
    gse=nr["GSE30573_mapping_status"]
    issues=[]
    # Rule1: discordant event must NOT claim directionally replicated/supported at individual level
    if disc_pe and claim in ("PRIMARY_SIGNIFICANT_EVENT","SUPPORTED_SENSITIVITY_EVENT"):
        issues.append("R1_discordant_event_with_individual_positive_claim")
    # Rule2: KR FDR<0.05 event must NOT say No_event_level_significance
    if fdrkr<0.05 and "No_event_level_significance" in nr["negative_evidence_summary"]:
        issues.append("R2_KR_FDR005_but_No_event_level_significance")
    # Rule3: GSE30573 unmapped event must NOT be LIMITED_N3
    if gse=="UNMAPPED_NOT_TESTABLE" and "LIMITED_N3" in nr["GSE30573_support_level"]:
        issues.append("R3_unmapped_but_LIMITED_N3")
    # Rule4: Tier D event must NOT have individual positive claim
    if tier.startswith("TIER_D") and claim in ("PRIMARY_SIGNIFICANT_EVENT","SUPPORTED_SENSITIVITY_EVENT","SET_MEMBER_DIRECTIONALLY_SUPPORTED"):
        # SET_MEMBER allowed? Tier D = discordant/no signal -> should be NO_INDIVIDUAL_CLAIM
        issues.append("R4_TierD_with_positive_claim")
    return issues

n_ok=0
for nr in rows:
    iss=allowed_claim_ok(nr)
    ok = len(iss)==0
    n_ok += 1 if ok else 0
    check.append([nr["HsaEX_ID"],nr["gene"],nr["final_evidence_tier"],nr["direction_concordant"],
                  nr["BH_FDR_KR"],nr["GSE30573_mapping_status"],nr["manuscript_claim_level"],
                  "OK" if ok else "ERROR",";".join(iss) if iss else "all_4_rules_ok"])

with open(os.path.join(MET,"MASTER_EVENT_TABLE_LOGIC_CHECK.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["HsaEX_ID","gene","final_evidence_tier","direction_concordant","BH_FDR_KR","GSE30573_mapping_status","manuscript_claim_level","logic_status","issues"])
    for a in check: w.writerow(a)

# phase
nf=sum(1 for a in check if a[7]=="ERROR")
with open(os.path.join(MET,"MASTER_TABLE_CHECK.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["key","value"])
    for k,v in [("N_EVENTS",len(rows)),("N_LOGIC_OK",n_ok),("N_LOGIC_ERROR",nf),
                ("XLSX_GENERATED",xlsx_ok),
                ("MASTER_TABLE_STATUS","OK" if nf==0 and len(rows)==19 else "HOLD")]:
        w.writerow([k,v])
print("EVENTS",len(rows),"LOGIC_OK",n_ok,"ERROR",nf,"XLSX",xlsx_ok)
for a in check:
    if a[7]=="ERROR": print("LOGIC ERROR:",a)
