#!/usr/bin/env python3
"""
Build one TSV mapping every reviewed (Swiss-Prot) human UniProt entry to its
protein names, accession, gene symbol and gene name.

Official sources:
  * UniProt REST API (UniProtKB, organism_id 9606, reviewed:true)
      https://rest.uniprot.org/uniprotkb/stream
  * HGNC complete set (approved gene symbols and full gene names)
      https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

Usage (Python 3.8+, standard library only):
    python swissprot_hgnc_map.py                      # downloads both sources
    python swissprot_hgnc_map.py -o my_map.tsv
    python swissprot_hgnc_map.py --uniprot-file u.tsv --hgnc-file hgnc.txt   # offline

Output columns:
    UniProt_Accession   primary accession (e.g. P04637)
    Entry_Name          UniProt entry name (e.g. P53_HUMAN)
    Recommended_Protein_Name  UniProt recommended name
    All_Protein_Names   full UniProt "Protein names" field (recommended + alternative,
                        EC numbers, [Cleaved into: ...] etc.)
    Gene_Symbol         UniProt primary gene name (= HGNC approved symbol for most entries)
    Gene_Name           HGNC approved full gene name (e.g. "tumor protein p53")
    Gene_Synonyms       other gene names UniProt lists for the entry
    HGNC_ID             HGNC identifier(s) cross-referenced by UniProt
    
"""
import argparse
import csv
import gzip
import io
import sys
import urllib.parse
import urllib.request
from collections import defaultdict

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_QUERY = "(organism_id:9606) AND (reviewed:true)"
# Requested in this order; the response is parsed by position, not header text.
UNIPROT_FIELDS = ["accession", "id", "protein_name", "gene_primary", "gene_names", "xref_hgnc"]
HGNC_URL = ("https://storage.googleapis.com/public-download-files/"
            "hgnc/tsv/tsv/hgnc_complete_set.txt")

csv.field_size_limit(sys.maxsize)


def fetch_text(url, params=None):
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    print(f"Downloading {url}", file=sys.stderr)
    req = urllib.request.Request(url, headers={"Accept-Encoding": "gzip"})
    with urllib.request.urlopen(req, timeout=600) as resp:
        data = resp.read()
    # UniProt's compressed=true payload may also be gzipped for transport,
    # so unwrap every gzip layer present.
    while data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("utf-8")


def read_text(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        return fh.read()


def recommended_name(full):
    """Text before the first top-level ' (' or ' [' group in UniProt's Protein names field."""
    depth = 0
    for i, ch in enumerate(full):
        if ch in "([":
            if depth == 0 and i > 0 and full[i - 1] == " ":
                return full[:i - 1].strip()
            depth += 1
        elif ch in ")]":
            depth = max(depth - 1, 0)
    return full.strip()


def split_ids(value):
    return [v.strip() for v in value.replace(",", ";").replace("|", ";").split(";") if v.strip()]


def load_hgnc(text):
    """Return {hgnc_id: (symbol, name)} and {uniprot_acc: [hgnc_id, ...]}."""
    by_id, by_acc = {}, defaultdict(list)
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    for row in reader:
        if row.get("status", "Approved") != "Approved":
            continue
        hid = row["hgnc_id"].strip()
        by_id[hid] = (row["symbol"].strip(), row["name"].strip())
        for acc in split_ids(row.get("uniprot_ids", "")):
            by_acc[acc].append(hid)
    return by_id, by_acc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="human_swissprot_protein_gene_map.tsv")
    ap.add_argument("--uniprot-file", help="local UniProt TSV with the fields listed in UNIPROT_FIELDS, in that order")
    ap.add_argument("--hgnc-file", help="local hgnc_complete_set.txt")
    args = ap.parse_args()

    uni_text = read_text(args.uniprot_file) if args.uniprot_file else fetch_text(
        UNIPROT_STREAM, {"query": UNIPROT_QUERY, "format": "tsv",
                         "fields": ",".join(UNIPROT_FIELDS), "compressed": "true"})
    try:
        hgnc_text = read_text(args.hgnc_file) if args.hgnc_file else fetch_text(HGNC_URL)
        hgnc_by_id, hgnc_by_acc = load_hgnc(hgnc_text)
    except Exception as exc:  # keep going with UniProt data alone
        print(f"WARNING: HGNC data unavailable ({exc}); Gene_Name will be empty.", file=sys.stderr)
        hgnc_by_id, hgnc_by_acc = {}, {}

    rows = list(csv.reader(io.StringIO(uni_text), delimiter="\t"))
    header, body = rows[0], [r for r in rows[1:] if r]
    if len(header) != len(UNIPROT_FIELDS):
        sys.exit(f"Unexpected UniProt columns: {header}")

    out_cols = ["UniProt_Accession", "Entry_Name", "Recommended_Protein_Name", "All_Protein_Names",
                "Gene_Symbol", "Gene_Name", "Gene_Synonyms", "HGNC_ID"]
    with_name = 0
    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(out_cols)
        for acc, entry, pnames, gsym, gnames, hgnc in sorted(body, key=lambda r: r[0]):
            hids = split_ids(hgnc) or hgnc_by_acc.get(acc, [])
            hids = [h for h in dict.fromkeys(hids) if h in hgnc_by_id] or list(dict.fromkeys(hids))
            gene_names = "; ".join(hgnc_by_id[h][1] for h in hids if h in hgnc_by_id)
            symbol = gsym.strip() or "; ".join(hgnc_by_id[h][0] for h in hids if h in hgnc_by_id)
            primaries = {s.strip() for s in symbol.split(";")}
            synonyms = "; ".join(g for g in gnames.split() if g not in primaries)
            with_name += bool(gene_names)
            w.writerow([acc, entry, recommended_name(pnames), pnames,
                        symbol, gene_names, synonyms, "; ".join(hids)])

    print(f"Wrote {len(body):,} entries to {args.output} "
          f"({with_name:,} with an HGNC gene name).", file=sys.stderr)


if __name__ == "__main__":
    main()