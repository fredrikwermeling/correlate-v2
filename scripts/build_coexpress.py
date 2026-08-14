#!/usr/bin/env python3
"""
Derive the CoExpress build from the V1 build.

CoExpress is V1 (so: no PCA/UMAP, no version badge or changelog, no mRNA
basis toggle) with EXPRESSION as the primary matrix. The app loads
metadata.json + geneEffects.bin.gz as its main matrix; in CoExpress those
hold the log2(TPM+1) expression data, so every feature — correlations,
network, clusters, mutation analysis, wiki, Cell Line Browser — runs on
expression with no code fork.

Run scripts/build_v1.py first (this reads V1's index.html / app.js), then:
    python3 scripts/build_coexpress.py

Data files are rebuilt separately by build_coexpress_data.py, only when the
DepMap release changes.
"""

import os
import re
import shutil

V2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V1 = "/Users/fredrikwermeling/Documents/correlate app feb 2026 (färdig)/correlation app"
CO = "/Users/fredrikwermeling/Documents/coexpress-app"

ASSETS = ["network_example.svg", "scatter_example.svg", "tsc_pathway.svg",
          "favicon-32x32.png", "favicon.ico",
          # UI images live in web_data but are build assets, not release data.
          "web_data/clb_banner.png"]

# Order matters: the specific phrases must be replaced before the generic
# "gene effect" sweep, which would otherwise consume them. The spaced forms
# are safe because identifiers are camelCase (geneEffects, colorByGeneEffect).
REPLACEMENTS = [
    # branding
    ('<title>Correlate</title>', '<title>CoExpress</title>'),
    ('src="web_data/correlate_logo.png" alt="Correlate - a Gene Effect Analysis Tool"',
     'src="web_data/coexpress_logo.png" alt="CoExpress - an Expression Analysis Tool"'),
    ('<p>Correlate | Part of the', '<p>CoExpress | Part of the'),
    # ---- AI-export sign conventions -------------------------------------
    # These sentences explain what the SIGN of a number means, and the generic
    # sweep turns them into falsehoods: "gene effect is NEGATIVE for a
    # dependency" becomes "expression is NEGATIVE for a dependency", which is
    # untrue of log₂(TPM+1) and would have an assistant reading every polarity
    # in the file backwards. They must be rewritten for expression, not swept.
    ('Polarity, and note the sign convention because it is easy to get backwards: gene effect is NEGATIVE for a dependency, so two genes both needed by the same lines move down together and their r is POSITIVE. Positive r therefore means partner and focal gene are co-essential, required by the same lines, which is the same-complex / same-pathway pattern: BRAF against MAP2K1, MAP2K2 and MAPK1 sits near r = +0.6 to +0.8. Negative r means the opposite, lines that depend on the partner do NOT depend on the focal gene, which is the buffering pattern seen between paralogues, where one covers for the other and only the line that has lost one needs the other. Neither sign is a mechanism on its own: check the pair against the geneEffect matrix, which carries every gene named here.',
     'Polarity: this is mRNA abundance on a log₂(TPM+1) scale, which is never negative, and HIGHER means more transcript. Positive r therefore means the two genes are CO-EXPRESSED, switched on and off together across the panel, which is the shared-transcriptional-programme pattern and is typical of genes in one complex, one pathway or one lineage programme. Negative r means they are reciprocal, high in different cell lines, which is typical of two competing differentiation states. Note what this is NOT: it says nothing about whether either gene is REQUIRED, because no knockout was performed anywhere in this dataset. Co-expression is not co-essentiality, and a strong correlation here is evidence of shared regulation, not of shared function. Check the pair against the expression matrix, which carries every gene named here.'),
    ('so a NEGATIVE delta_ge means the gene is more essential in the altered lines',
     'so a NEGATIVE delta_ge means the gene is expressed LOWER in the altered lines and a positive one means higher'),
    ('so a NEGATIVE delta means the gene is more essential in the MSI-high lines',
     'so a NEGATIVE delta means the gene is expressed LOWER in the MSI-high lines and a positive one means higher'),
    ("'CRISPR gene effect (Chronos): which genes the same cell lines depend on. NOT expression'",
     "'mRNA expression (log2 TPM+1): which genes are switched on together across the panel'"),
    ("POSITIVE r = the two are co-essential, needed by the same lines, the same-complex / same-pathway pattern, because gene effect is negative for a dependency and partners move down together; NEGATIVE r = the buffering pattern between paralogues, where the lines needing one do not need the other",
     "POSITIVE r = the two are CO-EXPRESSED, switched on and off together across the panel, which is the shared-programme pattern typical of one complex, one pathway or one lineage; NEGATIVE r = reciprocal, high in different lines, typical of competing differentiation states. This says nothing about either gene being REQUIRED: nothing was knocked out anywhere in this dataset"),
    ("'mRNA expression' : 'CRISPR gene effect'", "'mRNA expression' : 'mRNA expression'"),
    # semantics that contain the generic phrases
    ('how much the cell lines depend on that gene',
     'how strongly the cell lines express the gene (log₂ TPM+1)'),
    ('Negative means knocking the gene out slows growth.', 'Higher means more mRNA.'),
    ('How much the cell lines depend on', 'How strongly the cell lines express'),
    ('CRISPRGeneEffect (Chronos)', 'OmicsExpressionTPMLogp1 (log₂-TPM+1)'),
    ('CRISPRGeneEffect.csv <span style="color:#999;">(gene effects)</span>',
     'OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv <span style="color:#999;">(expression, the primary matrix)</span>'),
    # generic terminology
    ('Gene Effect', 'Expression'),
    ('Gene effect', 'Expression'),
    ('gene effect', 'expression'),
    ('gene-effect', 'expression'),
    # bare "GE" only through enumerated labels, never on its own
    (' (GE)', ''),
    ('GE mean', 'Expr mean'),
    ('GE SD', 'Expr SD'),
    ('Mean GE', 'Mean Expr'),
    ('GE=', 'Expr='),
    ('Δ GE', 'Δ Expr'),
    ('ΔGE', 'ΔExpr'),
    ('Color by GE', 'Color by expression'),
    ('colored by GE', 'colored by expression'),
    ('GE for gene(s)', 'Expression for gene(s)'),
    ('GE + expression', 'expression'),
    ('GE label', 'expression label'),
    ('GE score', 'expression value'),
]

# In CoExpress the primary matrix IS expression, so an axis choice of "GE" and
# one of "Expr" would read the same numbers. Collapse those duplicated options.
HTML_FIXES = [
    ("""<option value="ge">GE</option>
                                    <option value="expr">Expr</option>
                                    <option value="cn">CN</option>""",
     """<option value="ge">Expr</option>
                                    <option value="cn">CN</option>"""),
    ("""<option value="ge">GE</option>
                                <option value="expr">Expression</option>""",
     """<option value="ge">Expression</option>"""),
    # Cohort-specific counts: CoExpress runs on the larger expression cohort,
    # so the virus-transformation tallies written for the CRISPR panel are wrong here.
    ('53 lines in this panel are confirmed transformed by EBV, HPV16/18, HBV, KSHV or HTLV-1',
     '76 lines in this panel are confirmed transformed by EBV, HPV16/18, HBV, KSHV, HTLV-1 or adenovirus 5'),
    ('53 of the 1,208 lines carry one', '76 of the 1,719 lines carry one'),
    ('title="Switch between expression (CRISPR knockout) and expression (log2-TPM+1) as the y-axis data"',
     'title="Expression (log2-TPM+1) as the y-axis data"'),
    # the mutation-analysis Measure toggle would offer expression twice
    ("""<span style="display:inline-flex; align-items:center; gap:4px;">
                            <span style="font-size:10px; color:#9ca3af;">Measure:</span>""",
     """<span style="display:none;">
                            <span style="font-size:10px; color:#9ca3af;">Measure:</span>"""),
]

JS_FIXES = [
    (": 'Expr' : 'GE';", ": 'Expr' : 'Expr';"),
    ("alert('Please enter genes for GE/Expr axes.');", "alert('Please enter genes for both axes.');"),
]


def main():
    for name in ("index.html", "app.js"):
        src = os.path.join(V1, name)
        if not os.path.exists(src):
            raise SystemExit(f"{src} missing — run scripts/build_v1.py first")
        text = open(src).read()
        for old, new in REPLACEMENTS:
            text = text.replace(old, new)
        for old, new in (HTML_FIXES if name == "index.html" else JS_FIXES):
            text = text.replace(old, new)
        open(os.path.join(CO, name), "w").write(text)
        print(f"  wrote {name}")

    for a in ASSETS:
        s = os.path.join(V1, a)
        if os.path.exists(s):
            shutil.copy2(s, os.path.join(CO, a))

    left = sum(open(os.path.join(CO, n)).read().count("Gene Effect")
               + open(os.path.join(CO, n)).read().count("gene effect")
               for n in ("index.html", "app.js"))
    print(f"\nCoExpress written to {CO}")
    print(f"remaining 'gene effect' mentions: {left} (expect 0)")
    print("Data files: run build_coexpress_data.py when the DepMap release changes.")


if __name__ == "__main__":
    main()
