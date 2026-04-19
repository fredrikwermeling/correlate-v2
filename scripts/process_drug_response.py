#!/usr/bin/env python3
"""
Build web_data/drug_response.json from DepMap PRISM Repurposing Secondary AUC.

Input:
  correlate_v2/REPURPOSINGAUCMatrix.csv                       (AUC matrix; rows=cell, cols=DPC-id)
  correlate_v2/REPURPOSINGLog2ViabilityCollapsedConditions.csv (compound metadata)

Output:
  web_data/drug_response.json  — ~200 KB, per cell line sensitivity for ~150 curated drugs.

Structure:
{
  "panelSize": 150,
  "compounds": [
    { "id": "DPC-004744", "name": "OLAPARIB", "target": "PARP1/2", "moa": "PARP inhibitor",
      "indication": "BRCA-mut breast/ovarian/prostate/pancreatic",
      "auc": { "ACH-000001": 0.73, "ACH-000004": 0.42, ... },
      "mean": 0.81, "sd": 0.11 },
    ...
  ]
}

AUC interpretation: 0 = complete kill across dose range, 1 = no effect.
Low AUC = sensitive. High AUC = resistant.
"""

import csv
import json
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, '..')
AUC_CSV = os.path.join(REPO, 'REPURPOSINGAUCMatrix.csv')
CONDITIONS_CSV = os.path.join(REPO, 'REPURPOSINGLog2ViabilityCollapsedConditions.csv')
OUTPUT = os.path.join(REPO, 'web_data', 'drug_response.json')


# Curated panel. Target, MOA, indication are editorial — intended for wiki display.
# Drug names must match how they appear in REPURPOSING* (case-insensitive).
# Missing from the panel → just omitted with a warning.
CURATED_PANEL = [
    # ----- Kinase: RAS-MAPK -----
    ('VEMURAFENIB',   'BRAF V600',      'BRAF inhibitor',              'BRAF-mut melanoma, hairy-cell leukemia'),
    ('DABRAFENIB',    'BRAF V600',      'BRAF inhibitor',              'BRAF-mut melanoma, NSCLC, thyroid'),
    ('TRAMETINIB',    'MEK1/2',         'MEK inhibitor',               'BRAF-mut melanoma/NSCLC (combo); NF1'),
    ('COBIMETINIB',   'MEK1/2',         'MEK inhibitor',               'BRAF-mut melanoma (combo with vemurafenib)'),
    ('BINIMETINIB',   'MEK1/2',         'MEK inhibitor',               'NRAS-mut melanoma (combo)'),
    ('SELUMETINIB',   'MEK1/2',         'MEK inhibitor',               'NF1-associated plexiform neurofibroma'),
    # ----- EGFR / HER -----
    ('ERLOTINIB',     'EGFR',           'EGFR TKI',                    'EGFR-mut NSCLC; pancreatic'),
    ('GEFITINIB',     'EGFR',           'EGFR TKI',                    'EGFR-mut NSCLC'),
    ('OSIMERTINIB',   'EGFR T790M',     '3rd-gen EGFR TKI',            'EGFR-mut NSCLC (T790M+)'),
    ('AFATINIB',      'EGFR/HER2',      'pan-HER TKI',                 'EGFR-mut NSCLC'),
    ('LAPATINIB',     'HER2/EGFR',      'HER2 TKI',                    'HER2+ breast'),
    ('NERATINIB',     'HER2/EGFR',      'pan-HER TKI',                 'HER2+ breast'),
    ('TUCATINIB',     'HER2',           'HER2 TKI',                    'HER2+ breast'),
    # ----- ALK / ROS1 -----
    ('CRIZOTINIB',    'ALK/ROS1/MET',   'ALK/ROS1 inhibitor',          'ALK/ROS1-rearranged NSCLC'),
    ('CERITINIB',     'ALK',            'ALK inhibitor',               'ALK-rearranged NSCLC'),
    ('ALECTINIB',     'ALK',            'ALK inhibitor',               'ALK-rearranged NSCLC'),
    ('BRIGATINIB',    'ALK',            'ALK inhibitor',               'ALK-rearranged NSCLC'),
    ('LORLATINIB',    'ALK/ROS1',       'ALK/ROS1 inhibitor',          'ALK/ROS1-rearranged NSCLC'),
    ('ENTRECTINIB',   'TRK/ALK/ROS1',   'TRK/ALK/ROS1 inhibitor',      'NTRK-fusion tumors, ROS1+ NSCLC'),
    # ----- MET / RET / FGFR -----
    ('CAPMATINIB',    'MET',            'MET inhibitor',               'MET exon-14 NSCLC'),
    ('TEPOTINIB',     'MET',            'MET inhibitor',               'MET exon-14 NSCLC'),
    ('CABOZANTINIB',  'MET/VEGFR/RET',  'multi-kinase inhibitor',      'RCC, MTC, HCC'),
    ('SELPERCATINIB', 'RET',            'RET inhibitor',               'RET-fusion/mut cancers'),
    ('PRALSETINIB',   'RET',            'RET inhibitor',               'RET-fusion NSCLC, MTC'),
    ('ERDAFITINIB',   'FGFR1-4',        'FGFR inhibitor',              'FGFR-mut urothelial'),
    # ----- NTRK -----
    ('LAROTRECTINIB', 'TRKA/B/C',       'TRK inhibitor',               'NTRK-fusion solid tumors'),
    # ----- BCR-ABL / KIT / CSF1R -----
    ('IMATINIB',      'BCR-ABL/KIT',    'BCR-ABL/KIT TKI',             'CML, Ph+ ALL, GIST'),
    ('DASATINIB',     'BCR-ABL/SRC',    '2nd-gen BCR-ABL TKI',         'CML, Ph+ ALL'),
    ('NILOTINIB',     'BCR-ABL',        '2nd-gen BCR-ABL TKI',         'CML'),
    ('PONATINIB',     'BCR-ABL T315I',  '3rd-gen BCR-ABL TKI',         'T315I CML, Ph+ ALL'),
    ('BOSUTINIB',     'BCR-ABL',        'BCR-ABL TKI',                 'CML'),
    ('RIPRETINIB',    'KIT/PDGFRA',     'KIT/PDGFRA inhibitor',        'GIST (4th-line)'),
    # ----- VEGFR multi-kinase -----
    ('SORAFENIB',     'multi-kinase',   'multi-kinase inhibitor',      'HCC, RCC, thyroid'),
    ('SUNITINIB',     'VEGFR/PDGFR',    'multi-kinase inhibitor',      'RCC, GIST, pNET'),
    ('PAZOPANIB',     'VEGFR/PDGFR',    'multi-kinase inhibitor',      'RCC, soft-tissue sarcoma'),
    ('AXITINIB',      'VEGFR',          'VEGFR inhibitor',             'RCC'),
    ('LENVATINIB',    'VEGFR/FGFR',     'multi-kinase inhibitor',      'HCC, thyroid, RCC (combo)'),
    ('REGORAFENIB',   'multi-kinase',   'multi-kinase inhibitor',      'CRC, HCC, GIST'),
    # ----- BTK / BCL-2 / BCR pathway (hematologic) -----
    ('IBRUTINIB',     'BTK',            'BTK inhibitor',               'CLL, MCL, WM'),
    ('ACALABRUTINIB', 'BTK',            '2nd-gen BTK inhibitor',       'CLL, MCL'),
    ('ZANUBRUTINIB',  'BTK',            '2nd-gen BTK inhibitor',       'MCL, CLL, WM'),
    ('VENETOCLAX',    'BCL-2',          'BCL-2 inhibitor',             'CLL, AML (combo)'),
    ('NAVITOCLAX',    'BCL-2/BCL-xL',   'BCL-2/BCL-xL inhibitor',      'investigational, hematologic'),
    ('IDELALISIB',    'PI3K-delta',     'PI3K-δ inhibitor',            'CLL, FL'),
    ('DUVELISIB',     'PI3K-delta/gamma','PI3K-δ/γ inhibitor',         'CLL, FL'),
    ('COPANLISIB',    'PI3K pan',       'pan-PI3K inhibitor',          'FL (relapsed)'),
    # ----- PI3K / AKT / mTOR (solid) -----
    ('ALPELISIB',     'PI3Kalpha',      'PI3K-α inhibitor',            'PIK3CA-mut HR+ HER2- breast'),
    ('IPATASERTIB',   'AKT1/2/3',       'AKT inhibitor',               'triple-neg breast (investigational)'),
    ('EVEROLIMUS',    'mTORC1',         'mTOR inhibitor',              'HR+ breast, RCC, pNET'),
    ('TEMSIROLIMUS',  'mTORC1',         'mTOR inhibitor',              'RCC, MCL'),
    # ----- CDK4/6 (cell cycle) -----
    ('PALBOCICLIB',   'CDK4/6',         'CDK4/6 inhibitor',            'HR+ HER2- breast'),
    ('RIBOCICLIB',    'CDK4/6',         'CDK4/6 inhibitor',            'HR+ HER2- breast'),
    ('ABEMACICLIB',   'CDK4/6',         'CDK4/6 inhibitor',            'HR+ HER2- breast'),
    # ----- PARP (HR-deficient) -----
    ('OLAPARIB',      'PARP1/2',        'PARP inhibitor',              'BRCA-mut breast/ovary/prostate/pancreas'),
    ('NIRAPARIB',     'PARP1/2',        'PARP inhibitor',              'ovarian maintenance'),
    ('RUCAPARIB',     'PARP1/2',        'PARP inhibitor',              'BRCA-mut ovarian/prostate'),
    ('TALAZOPARIB',   'PARP1/2',        'PARP inhibitor',              'BRCA-mut breast'),
    # ----- JAK / STAT -----
    ('RUXOLITINIB',   'JAK1/2',         'JAK inhibitor',               'myelofibrosis, PV, GVHD'),
    ('TOFACITINIB',   'JAK1/3',         'JAK inhibitor',               'RA, UC (approved); cancer research'),
    # ----- IDH / TET / DNMT -----
    ('IVOSIDENIB',    'IDH1 mutant',    'IDH1-R132 inhibitor',         'IDH1-mut AML, cholangio'),
    ('ENASIDENIB',    'IDH2 mutant',    'IDH2-R140/R172 inhibitor',    'IDH2-mut AML'),
    ('AZACITIDINE',   'DNMT',           'hypomethylating agent',       'MDS, AML'),
    ('DECITABINE',    'DNMT',           'hypomethylating agent',       'MDS, AML'),
    # ----- FLT3 / KIT (AML) -----
    ('MIDOSTAURIN',   'FLT3/KIT',       'multi-kinase inhibitor',      'FLT3-mut AML, mastocytosis'),
    ('GILTERITINIB',  'FLT3',           'FLT3 inhibitor',              'FLT3-mut relapsed AML'),
    ('QUIZARTINIB',   'FLT3',           'FLT3 inhibitor',              'FLT3-ITD AML'),
    # ----- Chemo: DNA damage -----
    ('CISPLATIN',     'DNA crosslink',  'platinum',                    'many solid tumors'),
    ('CARBOPLATIN',   'DNA crosslink',  'platinum',                    'ovarian, NSCLC, many'),
    ('OXALIPLATIN',   'DNA crosslink',  'platinum',                    'CRC, gastric, pancreatic'),
    ('DOXORUBICIN',   'TOP2A',          'anthracycline',               'many hematologic/solid'),
    ('EPIRUBICIN',    'TOP2A',          'anthracycline',               'breast, gastric'),
    ('DAUNORUBICIN',  'TOP2A',          'anthracycline',               'AML, ALL'),
    ('ETOPOSIDE',     'TOP2A',          'topo-II poison',              'SCLC, testicular, lymphoma'),
    ('IRINOTECAN',    'TOP1',           'topo-I poison',               'CRC, SCLC, pancreatic'),
    ('TOPOTECAN',     'TOP1',           'topo-I poison',               'SCLC, ovarian'),
    ('TEMOZOLOMIDE',  'DNA alkyl',      'alkylating agent',            'GBM, melanoma'),
    ('DACARBAZINE',   'DNA alkyl',      'alkylating agent',            'melanoma, HL, sarcoma'),
    ('CYCLOPHOSPHAMIDE','DNA crosslink','alkylating agent',            'many heme/solid'),
    ('IFOSFAMIDE',    'DNA crosslink',  'alkylating agent',            'sarcoma, testicular'),
    ('BENDAMUSTINE',  'DNA alkyl',      'alkylating agent',            'CLL, NHL'),
    ('CHLORAMBUCIL',  'DNA alkyl',      'alkylating agent',            'CLL'),
    ('CARMUSTINE',    'DNA alkyl',      'nitrosourea',                 'GBM, HL'),
    ('LOMUSTINE',     'DNA alkyl',      'nitrosourea',                 'brain tumors'),
    # ----- Chemo: antimetabolites -----
    ('GEMCITABINE',   'dNTP pool',      'nucleoside analog',           'pancreatic, NSCLC, bladder'),
    ('CYTARABINE',    'dNTP / pol',     'nucleoside analog',           'AML, ALL, lymphoma'),
    ('FLUDARABINE',   'nucleoside',     'purine analog',               'CLL, conditioning'),
    ('METHOTREXATE',  'DHFR',           'antifolate',                  'ALL, lymphoma, autoimmune'),
    ('PEMETREXED',    'TYMS/DHFR',      'antifolate',                  'NSCLC, mesothelioma'),
    ('FLUOROURACIL',  'TYMS',           'pyrimidine analog',           'CRC, breast, many'),
    ('CAPECITABINE',  'TYMS (prodrug)', 'pyrimidine analog',           'CRC, breast'),
    # ----- Chemo: microtubule -----
    ('PACLITAXEL',    'tubulin',        'taxane',                      'breast, ovary, NSCLC, many'),
    ('DOCETAXEL',     'tubulin',        'taxane',                      'breast, NSCLC, prostate'),
    ('VINCRISTINE',   'tubulin',        'vinca alkaloid',              'ALL, lymphoma'),
    ('VINBLASTINE',   'tubulin',        'vinca alkaloid',              'HL, testicular'),
    ('VINORELBINE',   'tubulin',        'vinca alkaloid',              'NSCLC, breast'),
    ('ERIBULIN',      'tubulin',        'halichondrin analog',         'breast, liposarcoma'),
    # ----- Proteasome / IMiDs -----
    ('BORTEZOMIB',    'proteasome',     'proteasome inhibitor',        'MM, MCL'),
    ('CARFILZOMIB',   'proteasome',     'proteasome inhibitor',        'MM'),
    ('IXAZOMIB',      'proteasome',     'proteasome inhibitor',        'MM'),
    ('LENALIDOMIDE',  'CRBN',           'IMiD (cereblon degrader)',    'MM, MDS-del5q, lymphoma'),
    ('POMALIDOMIDE',  'CRBN',           'IMiD',                        'MM'),
    ('THALIDOMIDE',   'CRBN',           'IMiD',                        'MM, ENL'),
    # ----- Epigenetic (non-IDH / DNMT) -----
    ('VORINOSTAT',    'HDAC pan',       'HDAC inhibitor',              'CTCL'),
    ('PANOBINOSTAT',  'HDAC pan',       'HDAC inhibitor',              'MM (combo)'),
    ('ROMIDEPSIN',    'HDAC1/2',        'HDAC inhibitor',              'CTCL, PTCL'),
    ('BELINOSTAT',    'HDAC pan',       'HDAC inhibitor',              'PTCL'),
    ('TAZEMETOSTAT',  'EZH2',           'EZH2 inhibitor',              'epithelioid sarcoma, FL (EZH2-mut)'),
    # ----- Hormonal -----
    ('TAMOXIFEN',     'ER',             'SERM',                        'ER+ breast'),
    ('FULVESTRANT',   'ER',             'SERD',                        'ER+ breast'),
    ('ENZALUTAMIDE',  'AR',             'AR antagonist',               'prostate (CRPC)'),
    ('ABIRATERONE',   'CYP17A1',        'androgen biosynth inhibitor', 'prostate (CRPC)'),
    ('APALUTAMIDE',   'AR',             'AR antagonist',               'prostate'),
    ('DARolutamide',  'AR',             'AR antagonist',               'prostate'),
    ('LETROZOLE',     'aromatase',      'aromatase inhibitor',         'ER+ breast (postmenopause)'),
    ('ANASTROZOLE',   'aromatase',      'aromatase inhibitor',         'ER+ breast (postmenopause)'),
    ('EXEMESTANE',    'aromatase',      'steroidal aromatase inhib.',  'ER+ breast'),
    # ----- Hedgehog -----
    ('VISMODEGIB',    'SMO',            'SMO/Hedgehog inhibitor',      'basal cell carcinoma'),
    ('SONIDEGIB',     'SMO',            'SMO/Hedgehog inhibitor',      'basal cell carcinoma'),
    # ----- Apoptosis / MDM2 -----
    ('IDASANUTLIN',   'MDM2',           'MDM2 inhibitor',              'TP53-WT (investigational)'),
    # ----- Tool compounds / preclinical -----
    ('JQ1',           'BRD4',           'BET bromodomain inhibitor',   'preclinical, MYC-driven'),
    ('ONALESPIB',     'HSP90',          'HSP90 inhibitor',             'preclinical'),
    ('TANESPIMYCIN',  'HSP90',          'HSP90 inhibitor',             'preclinical'),
    # ----- Misc targeted -----
    ('AFLIBERCEPT',   'VEGF-A',         'VEGF trap',                   'CRC'),
    ('ENCORAFENIB',   'BRAF V600',      'BRAF inhibitor',              'BRAF-mut melanoma/CRC (combo)'),
    ('DUSIGITUMAB',   'IGF-1R',         'anti-IGF1R antibody',         'investigational'),
    ('AG-270',        'MAT2A',          'MAT2A inhibitor',             'MTAP-del cancers (preclinical)'),
    ('MK-6482',       'HIF-2alpha',     'HIF-2α inhibitor',            'VHL-associated RCC'),
    ('SIROLIMUS',     'mTORC1',         'mTOR inhibitor',              'transplant, LAM'),
]


def norm(name):
    """Normalize drug name for matching."""
    return re.sub(r'[^A-Z0-9]', '', (name or '').upper())


def find_id(name, compound_map):
    """Try exact, then fuzzy match. Returns (compound_id, actual_name) or (None, None)."""
    key = norm(name)
    if key in compound_map:
        return compound_map[key]
    # Fuzzy: any compound that starts with the key or vice versa
    for norm_name, (dpc, actual) in compound_map.items():
        if norm_name == key or norm_name.startswith(key) or key.startswith(norm_name):
            return (dpc, actual)
    return (None, None)


def main():
    # Load compound ID → name
    compound_map = {}  # norm_name → (compound_id, actual_name)
    with open(CONDITIONS_CSV) as f:
        for row in csv.DictReader(f):
            n = norm(row['CompoundName'])
            if n and n not in compound_map:
                compound_map[n] = (row['CompoundID'], row['CompoundName'])
    print(f"Loaded {len(compound_map)} unique compounds from conditions")

    # Load AUC matrix
    print(f"Reading AUC matrix from {AUC_CSV}...")
    with open(AUC_CSV) as f:
        reader = csv.reader(f)
        header = next(reader)  # first col = '', rest = DPC-ids
        dpc_ids = header[1:]
        dpc_to_idx = {d: i for i, d in enumerate(dpc_ids)}
        auc_by_cl = {}  # cl_id → list of floats aligned to dpc_ids
        for row in reader:
            cl = row[0]
            vals = []
            for v in row[1:]:
                vals.append(float(v) if v else float('nan'))
            auc_by_cl[cl] = vals
    print(f"  {len(auc_by_cl)} cell lines × {len(dpc_ids)} compounds")

    # Build compound-centric output
    compounds_out = []
    missed = []
    for (name, target, moa, indication) in CURATED_PANEL:
        dpc, actual = find_id(name, compound_map)
        if not dpc:
            missed.append(name)
            continue
        col_idx = dpc_to_idx.get(dpc)
        if col_idx is None:
            missed.append(f"{name} (id {dpc} not in AUC matrix)")
            continue
        auc = {}
        vals_for_stats = []
        for cl, row in auc_by_cl.items():
            v = row[col_idx]
            if not math.isnan(v):
                auc[cl] = round(v, 4)
                vals_for_stats.append(v)
        if len(vals_for_stats) < 10:
            missed.append(f"{name} (only {len(vals_for_stats)} cell lines with data)")
            continue
        mean = sum(vals_for_stats) / len(vals_for_stats)
        sd = math.sqrt(sum((v - mean) ** 2 for v in vals_for_stats) / (len(vals_for_stats) - 1)) if len(vals_for_stats) > 1 else 0
        compounds_out.append({
            'id': dpc,
            'name': actual.upper(),
            'target': target,
            'moa': moa,
            'indication': indication,
            'auc': auc,
            'mean': round(mean, 4),
            'sd': round(sd, 4),
            'nCells': len(vals_for_stats),
        })

    print(f"\nPanel: {len(compounds_out)} compounds resolved, {len(missed)} missed/skipped")
    for m in missed:
        print(f"  ✗ {m}")

    out = {
        'panelSize': len(compounds_out),
        'dataSource': 'DepMap PRISM Repurposing Secondary 25Q2 (AUC; 0 = complete kill, 1 = no effect)',
        'compounds': compounds_out,
    }
    with open(OUTPUT, 'w') as f:
        json.dump(out, f)
    size_mb = os.path.getsize(OUTPUT) / (1024 * 1024)
    print(f"\nWrote {OUTPUT} ({size_mb:.2f} MB)")


if __name__ == '__main__':
    main()
