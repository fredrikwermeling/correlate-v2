"""
Make every gene in the Green Listed CRISPR libraries answerable in the gene box
-> web_data/library_genes.json.

Libraries (greenlistedv2/libraries): human Brunello, Gattinara, GeCKO v2,
Jacquere, VBC, MinLibCas9, TKOv3, Yusa v1; mouse Brie, GeCKO v2, Gouda,
Julianna, VBC, mTKO, Yusa v2.

  alias        {ALIAS: DEPMAP_GENE}. Names the app's own tables (synonyms.json,
               orthologs.json) cannot resolve but Green Listed's human / mouse
               synonym tables can, kept only when they land on exactly one gene
               with CRISPR data (an alias shared by two genes is left out rather
               than guessed). Mostly old mouse names (RIKEN clone IDs, Gm genes
               since renamed).
  notScreened  {SYMBOL: code} for library genes that resolve to no gene with
               CRISPR data, so the app can say why instead of "not found":
                 m  microRNA
                 o  olfactory receptor
                 g  mouse gene with no human ortholog
                 n  gene not screened by DepMap

Multi-gene library entries ("EBP|nan", "Becn1|Cntd1") are split by the app at
run time; their parts are covered here like any other symbol.

Run: python3 scripts/build_library_genes.py   (idempotent; needs openpyxl)
"""
import csv
import json
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(ROOT, "web_data")
LIB = os.path.expanduser("~/Documents/greenlistedv2/libraries")

TSV = {"Brunello (human).txt": 1, "Gattinara (human).txt": 1, "GeCKO v2 (human) A+B.txt": 0,
       "Jacquere (human).txt": 1, "VBC (human).txt": 0,
       "Brie (mouse).txt": 1, "Gouda (mouse).txt": 1, "Julianna (mouse).txt": 1,
       "VBC (mouse).txt": 0, "GeCKO v2 (mouse) A+B.txt": 0}
XLSX = {"minlibcas9_raw.xlsx": 3, "tkov3_raw.xlsx": 0, "yusa_human_v1_raw.xlsx": 1,
        "mtko_raw.xlsx": 2, "yusa_mouse_v2_raw.xlsx": 1}


def is_mouse(name):
    return "mouse" in name or name.startswith(("mtko", "yusa_mouse"))


def library_symbols():
    import openpyxl
    human, mouse = set(), set()
    for fn, col in TSV.items():
        with open(os.path.join(LIB, fn), newline="") as fh:
            r = csv.reader(fh, delimiter="\t")
            next(r)
            for row in r:
                if len(row) > col and row[col].strip():
                    (mouse if is_mouse(fn) else human).add(row[col].strip())
    for fn, col in XLSX.items():
        wb = openpyxl.load_workbook(os.path.join(LIB, fn), read_only=True)
        ws = wb[wb.sheetnames[0]]
        for i, row in enumerate(ws.iter_rows(values_only=True)):
            if i and len(row) > col and row[col] and str(row[col]).strip():
                (mouse if is_mouse(fn) else human).add(str(row[col]).strip())
    split = lambda s: {p for x in s for p in x.split("|") if p and p.lower() != "nan"}
    return split(human), split(mouse)


def gl_table(fn):
    d = defaultdict(set)
    with open(os.path.join(LIB, fn), newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[r["Gene Synonym"].strip().upper()].add(r["Gene name"].strip().upper())
    return d


def main():
    ge = {g.upper() for g in json.load(open(os.path.join(WEB, "metadata.json")))["genes"]}
    syn = json.load(open(os.path.join(WEB, "synonyms.json")))
    m2h = json.load(open(os.path.join(WEB, "orthologs.json")))["mouseToHuman"]
    syn_d = lambda u: (syn.get(u) or {}).get("d", "").upper()

    def own(u, mouse):
        """What the app already resolves without this file."""
        if u in ge:
            return u
        if mouse and m2h.get(u, "").upper() in ge:
            return m2h[u].upper()
        return syn_d(u) if syn_d(u) in ge else None

    def via(name):
        for cand in (name, m2h.get(name, "").upper(), syn_d(name)):
            if cand in ge:
                return cand
        return None

    human, mouse = library_symbols()
    tables = {"human": gl_table("human synonym.txt"), "mouse": gl_table("mouse synonym.txt")}

    alias = {}
    for species, table in tables.items():
        for a, names in table.items():
            if own(a, species == "mouse"):
                continue
            targets = {t for t in (via(n) for n in names) if t}
            # The gene box splits on whitespace, so an alias with a space in
            # it can never be typed as one name.
            if len(targets) == 1 and not any(ch.isspace() for ch in a):
                alias.setdefault(a, next(iter(targets)))

    not_screened = {}
    for species, symbols in (("human", human), ("mouse", mouse)):
        for s in symbols:
            u = s.upper()
            if own(u, species == "mouse") or u in alias:
                continue
            if u.startswith(("HSA-MIR", "MMU-MIR", "MIR")):
                code = "m"
            elif u.startswith("OLFR") or (u.startswith("OR") and u[2:3].isdigit()):
                code = "o"
            elif species == "mouse":
                code = "g"
            else:
                code = "n"
            not_screened.setdefault(u, code)

    out = {"_doc": __doc__.strip().split("\n\n")[0], "alias": alias, "notScreened": not_screened}
    path = os.path.join(WEB, "library_genes.json")
    json.dump(out, open(path, "w"), separators=(",", ":"))
    print(f"library symbols: {len(human)} human, {len(mouse)} mouse")
    print(f"alias: {len(alias)}   notScreened: {len(not_screened)}   ({os.path.getsize(path)/1e3:.0f} KB)")
    for sp, symbols in (("human", human), ("mouse", mouse)):
        ok = sum(1 for s in symbols if own(s.upper(), sp == "mouse") or s.upper() in alias)
        print(f"  {sp}: {ok/len(symbols)*100:.1f}% of library genes resolve to a gene with CRISPR data")


if __name__ == "__main__":
    main()
