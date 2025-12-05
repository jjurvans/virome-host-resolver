#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 09:59:34 2025

@author: jaanajurvansuu

End-to-end host inference and taxonomy resolution pipeline.

Steps:
1) Read detection tables, pool accessions/labels, infer most probable host using NCBI (Entrez).
2) Extract unique host names.
3) Resolve host names to canonical scientific names + taxonomy using GNverifier + GBIF.

Requires:
  - pandas
  - biopython (Bio.Entrez)
  - requests
  - tqdm
"""

import re
import time
import collections
from typing import List, Tuple, Optional, Dict

import pandas as pd
import requests
from tqdm import tqdm
from Bio import Entrez

# =========================
# USER SETTINGS
# =========================

# Input virus detection tables (TSV)
DETECTIONS_00 = "detected_viruses.tsv"

# Entrez credentials (NCBI)
ENTREZ_EMAIL = "X"  # <-- set this
ENTREZ_API_KEY = "Y"
NCBI_DELAY = 0.34                        # seconds between NCBI calls
RETMAX_SEARCH = 5                        # max nuccore records in label search

# Output files
POOLED_HOSTS_TSV = "most_probable_hosts_pooled1.tsv"
POOL_MEMBERS_TSV = "pool_members1.tsv"
UNIQUE_HOSTS_TXT = "unique_most_probable_hosts1.txt"
RESOLVED_HOSTS_CSV = "resolved_hosts_full_taxonomy_final1.csv"

# Column name hints (script auto-detects first present)
ACC_COLS  = ["accession", "Accession", "nuccore", "genbank_accession"]
NAME_COLS = ["Name", "species", "label", "description", "Virus"]
ASM_COLS  = ["Assembly", "assembly", "set", "Set"]

# External APIs
GNVERIFIER_URL = "https://verifier.globalnames.org/api/v1/verifications"
GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"

# =========================
# SMALL HELPERS
# =========================

def pause():
    time.sleep(NCBI_DELAY)


def to_text(x):
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return ""
    return str(x).strip()


def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


ACCN_RE = re.compile(r'^[A-Z]{1,4}_?\d{3,9}(\.\d+)?$', re.I)

# =========================
# NCBI HELPERS
# =========================

Entrez.email = ENTREZ_EMAIL
if ENTREZ_API_KEY and ENTREZ_API_KEY.strip():
    Entrez.api_key = ENTREZ_API_KEY


def efetch_nuccore_xml(uid_or_acc: str):
    try:
        h = Entrez.efetch(db="nuccore", id=uid_or_acc, rettype="gb", retmode="xml")
        gb = Entrez.read(h)
        pause()
        return gb
    except Exception:
        return None


def esearch_nuccore_ids_by_label(label: str, retmax=5) -> List[str]:
    q = f'"{label}"[All Fields] AND txid10239[Organism:exp]'
    try:
        r = Entrez.read(Entrez.esearch(db="nuccore", term=q, retmax=retmax))
        pause()
        return r.get("IdList", [])
    except Exception:
        return []


def taxonomy_host_taxid(name: str) -> Optional[str]:
    if not name:
        return None
    try:
        r = Entrez.read(Entrez.esearch(db="taxonomy", term=f'"{name}"[All Names]', retmax=1))
        pause()
        ids = r.get("IdList", [])
        return ids[0] if ids else None
    except Exception:
        return None


def extract_hosts_from_gbset(gbset) -> List[str]:
    """Extract host/isolation_source from GenBank XML record."""
    hosts = []
    if not gbset:
        return hosts
    for rec in gbset:
        for feat in rec.get("GBSeq_feature-table", []):
            if feat.get("GBFeature_key") != "source":
                continue
            quals = {
                q.get("GBQualifier_name"): to_text(q.get("GBQualifier_value"))
                for q in feat.get("GBFeature_quals", []) if q
            }
            host = quals.get("host", "")
            isol = quals.get("isolation_source", "")
            if host:
                hosts.append(host)
            elif isol:
                hosts.append(isol)
    return hosts


def majority_vote(names: List[str]) -> Tuple[str, float]:
    names = [n for n in (to_text(x) for x in names) if n]
    if not names:
        return "", 0.0
    c = collections.Counter(n.lower() for n in names)
    winner, n = c.most_common(1)[0]
    conf = n / max(1, sum(c.values()))
    nice = next((x for x in names if x.lower() == winner), winner.title())
    return nice, conf

# =========================
# GROUPING / POOLING
# =========================

def make_group_id(row, acc_col, asm_col, name_col) -> str:
    asm = to_text(row.get(asm_col, "")) if asm_col else ""
    if asm.lower().startswith("set:"):
        return f"SET:{asm.split(':', 1)[1]}"
    acc = to_text(row.get(acc_col, "")) if acc_col else ""
    if ACCN_RE.match(acc):
        return f"ACC:{acc}"
    label = to_text(row.get(name_col, "")) if name_col else ""
    return f"LAB:{label}" if label else f"ROW:{row.name}"


def collect_group_records(df, acc_col, asm_col, name_col):
    groups = {}
    members = []  # for mapping file
    for i, r in df.iterrows():
        gid = make_group_id(r, acc_col, asm_col, name_col)
        d = groups.setdefault(gid, {"labels": set(), "accessions": set(), "rows": []})
        label = to_text(r.get(name_col, "")) if name_col else ""
        acc = to_text(r.get(acc_col, "")) if acc_col else ""
        if label:
            d["labels"].add(label)
        if ACCN_RE.match(acc):
            d["accessions"].add(acc)
        d["rows"].append(i)
        members.append({
            "group_id": gid,
            "row_index": i,
            "label": label,
            "accession": acc,
            "assembly": to_text(r.get(asm_col, "")) if asm_col else ""
        })
    return groups, pd.DataFrame(members)


def infer_group_host(gid: str, ginfo: Dict) -> Tuple[str, str, float, str, int, List[str]]:
    """
    Returns:
      host_name, host_taxid, confidence, evidence_tier, n_accessions_used, all_hosts_list
    T1: pooled accession evidence across all accessions
    T2: label-based search (majority across a few nuccore records)
    T0: unknown
    """
    accs = sorted(ginfo["accessions"])
    host_votes: List[str] = []
    used = 0

    # Tier 1: accessions
    for acc in accs:
        gb = efetch_nuccore_xml(acc)
        if gb:
            used += 1
            host_votes.extend(extract_hosts_from_gbset(gb))

    if host_votes:
        host, conf = majority_vote(host_votes)
        return host, (taxonomy_host_taxid(host) or ""), conf, "T1_accession_pool", used, host_votes

    # Tier 2: label-based
    label = next(iter(ginfo["labels"]), "")
    if label:
        ids = esearch_nuccore_ids_by_label(label, retmax=RETMAX_SEARCH)
        for uid in ids:
            gb = efetch_nuccore_xml(uid)
            host_votes.extend(extract_hosts_from_gbset(gb))
        if host_votes:
            host, conf = majority_vote(host_votes)
            return host, (taxonomy_host_taxid(host) or ""), conf, "T2_label_search", used, host_votes

    # Tier 0: unknown
    return "", "", 0.0, "T0_unknown", used, []

# =========================
# HOST NAME RESOLUTION
# =========================

HUMAN_HINTS = [
    "homo sapiens", "human", "child", "boy", "girl", "man", "woman",
    "infant", "patient", "nasopharyngeal", "cervical", "vulvar",
    "nasal lavage", "csf", "cerebrospinal", "polio", "diabetes",
    "skin", "wart", "keratosis", "carcinoma", "sarcoma", "neurofibroma",
]


def map_to_human_if_applicable(name_lower: str) -> Optional[str]:
    for hint in HUMAN_HINTS:
        if hint in name_lower:
            return "Homo sapiens"
    return None


COMMON_NAME_MAP = {
    # cattle
    "cow": "Bos taurus",
    "cow ": "Bos taurus",
    "cattle": "Bos taurus",
    "cattle blood": "Bos taurus",
    "dairy cattle": "Bos taurus",
    "calf": "Bos taurus",
    "korean native calf": "Bos taurus",
    # pigs
    "pig": "Sus scrofa domesticus",
    "pig ": "Sus scrofa domesticus",
    "pig faeces": "Sus scrofa domesticus",
    "pig feces": "Sus scrofa domesticus",
    "piglet": "Sus scrofa domesticus",
    "swine": "Sus scrofa domesticus",
    "porcine": "Sus scrofa domesticus",
    "fecal pig": "Sus scrofa domesticus",
    "domestic pig": "Sus scrofa domesticus",
    "wild boar": "Sus scrofa",
    # dogs
    "dog": "Canis lupus familiaris",
    "dog ": "Canis lupus familiaris",
    "mongrel dog": "Canis lupus familiaris",
    "domestic dog": "Canis lupus familiaris",
    "dog (rottweiler mix adult)": "Canis lupus familiaris",
    "canine": "Canis lupus familiaris",
    "canine lupus familaris": "Canis lupus familiaris",
    # cats
    "cat": "Felis catus",
    "cat ": "Felis catus",
    "domestic cat": "Felis catus",
    "domestic cat with diarrhea": "Felis catus",
    "domestic cat; sex: male": "Felis catus",
    "feline": "Felis catus",
    # goats / sheep / equids
    "goat": "Capra hircus",
    "lamb": "Ovis aries",
    "donkey": "Equus africanus asinus",
    "foal (yili horse)": "Equus caballus",
    # other mammals
    "amur tiger": "Panthera tigris altaica",
    "raccoon dog": "Nyctereutes procyonoides",
    "rat": "Rattus norvegicus",
    # birds
    "chicken": "Gallus gallus domesticus",
    "broiler chicken": "Gallus gallus domesticus",
    "local chicken": "Gallus gallus domesticus",
    "chicken breed luqin": "Gallus gallus domesticus",
    "chicken; breed: yellow feather broiler": "Gallus gallus domesticus",
    "fowl": "Gallus gallus domesticus",
    "duck": "Anas platyrhynchos",
    "turkey": "Meleagris gallopavo",
    "pigeon": "Columba livia",
    "pigeon ": "Columba livia",
    "carrier pigeon": "Columba livia",
    "fancy pigeon": "Columba livia",
    # fish
    "rainbow trout tissue": "Oncorhynchus mykiss",
    # invertebrates
    "house cricket": "Acheta domesticus",
    # plants (subset)
    "barley": "Hordeum vulgare",
    "wheat": "Triticum aestivum",
    "maize": "Zea mays",
    "cucumber": "Cucumis sativus",
    "muskmelon": "Cucumis melo",
    "pea": "Pisum sativum",
    "pepper": "Capsicum annuum",
    "bell pepper": "Capsicum annuum",
    "grapevine": "Vitis vinifera",
    "grapevine cultivar pinot gris": "Vitis vinifera",
    "groundnut": "Arachis hypogaea",
    "strawberry": "Fragaria × ananassa",
    "sugar beet": "Beta vulgaris",
    "common reed": "Phragmites australis",
    "red currant cv. holandsky cerveny": "Ribes rubrum",
    "tomato; cultivar=candela": "Solanum lycopersicum",
    "pumpkin sp.": "Cucurbita pepo",
    "squash": "Cucurbita pepo",
    "zucchini": "Cucurbita pepo",
}

SAMPLE_ENV_KEYWORDS = [
    "sewage", "wastewater", "influent", "oxidation pond",
    "sewage contaminated water", "sewage oxidation pond",
    "nest material", "mud",
    "intestinal content",
    "pooled weed plants",
]


def is_sample_environment_only(name_lower: str) -> bool:
    for kw in SAMPLE_ENV_KEYWORDS:
        if kw in name_lower:
            return True
    if name_lower.strip() in {"stool", "feces", "faeces", "blood"}:
        return True
    return False


def gnverifier_canonical(name: str) -> Optional[str]:
    payload = {"nameStrings": [name]}
    try:
        resp = requests.post(GNVERIFIER_URL, json=payload, timeout=20)
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        return None

    names_list = data.get("names") or []
    if not names_list:
        return None

    name_record = names_list[0]
    best = name_record.get("bestResult")
    if not best:
        return None

    canonical = (
        best.get("currentCanonicalSimple")
        or best.get("matchedCanonicalSimple")
        or best.get("currentCanonicalFull")
        or best.get("matchedCanonicalFull")
        or name_record.get("name")
        or name_record.get("nameString")
    )
    return canonical


def gbif_taxonomy(scientific_name: str) -> Optional[dict]:
    try:
        resp = requests.get(GBIF_MATCH_URL, params={"name": scientific_name}, timeout=20)
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        return None

    if not data.get("usageKey"):
        return None

    confidence = data.get("confidence", 0)
    match_type = data.get("matchType", "NONE")
    if match_type == "NONE" or confidence < 60:
        return None

    return {
        "scientificName": data.get("canonicalName") or data.get("scientificName"),
        "rank": data.get("rank", ""),
        "kingdom": data.get("kingdom", ""),
        "phylum": data.get("phylum", ""),
        "class": data.get("class", ""),
        "order": data.get("order", ""),
        "family": data.get("family", ""),
        "genus": data.get("genus", ""),
        "species": data.get("species", ""),
        "confidence": confidence,
        "matchType": match_type,
        "status": data.get("status", ""),
        "synonym": data.get("synonym", False),
    }

# =========================
# MAIN PIPELINE
# =========================

def main():
    # --- 1) Load and pool detections ---
    detections = pd.read_csv(DETECTIONS_00, sep="\t")

    acc_col = pick_col(detections, ACC_COLS)
    name_col = pick_col(detections, NAME_COLS)
    asm_col = pick_col(detections, ASM_COLS)

    print(f"Detected columns → accession:{acc_col or '-'}  name:{name_col or '-'}  assembly:{asm_col or '-'}")

    groups, members = collect_group_records(detections, acc_col, asm_col, name_col)
    print(f"Pooled groups: {len(groups)}  (from {len(detections)} rows)")

    rows = []
    for k, (gid, ginfo) in enumerate(groups.items(), 1):
        host, taxid, conf, tier, nacc, votes = infer_group_host(gid, ginfo)
        all_hosts_pretty = "; ".join(sorted({to_text(v) for v in votes if to_text(v)}))
        unique_host_count = len({to_text(v).lower() for v in votes if to_text(v)})
        rows.append({
            "group_id": gid,
            "representative_label": next(iter(ginfo["labels"]), ""),
            "n_rows": len(ginfo["rows"]),
            "n_accessions": len(ginfo["accessions"]),
            "n_accessions_fetched": nacc,
            "accessions": ";".join(sorted(ginfo["accessions"])) if ginfo["accessions"] else "",
            "most_probable_host": host,
            "host_taxid": taxid,
            "confidence": round(conf, 3),
            "evidence_tier": tier,
            "unique_host_count": unique_host_count,
            "all_hosts_seen": all_hosts_pretty,
        })
        if k % 25 == 0 or k == len(groups):
            print(f"\rProcessed groups: {k}/{len(groups)}", end="")
    print()

    pooled = pd.DataFrame(rows).sort_values(
        ["evidence_tier", "confidence", "unique_host_count", "group_id"],
        ascending=[True, False, True, True]
    )
    pooled.to_csv(POOLED_HOSTS_TSV, sep="\t", index=False)
    members.to_csv(POOL_MEMBERS_TSV, sep="\t", index=False)

    print(f"Wrote: {POOLED_HOSTS_TSV}")
    print(f"Wrote: {POOL_MEMBERS_TSV}")
    print(pooled["evidence_tier"].value_counts(dropna=False).to_string())

    # --- 2) Extract unique host names ---
    hosts = (
        pooled["most_probable_host"]
        .astype(str)
        .str.strip()
        .replace("", pd.NA)
        .dropna()
        .unique()
    )
    hosts_sorted = sorted(hosts, key=lambda x: x.lower())

    with open(UNIQUE_HOSTS_TXT, "w", encoding="utf-8") as f:
        for h in hosts_sorted:
            f.write(h + "\n")

    print(f"Unique hosts: {len(hosts_sorted)}")
    print(f"Saved: {UNIQUE_HOSTS_TXT}")
    print("Preview:")
    for h in hosts_sorted[:20]:
        print(" -", h)

    # --- 3) Resolve host taxonomy via GNverifier + GBIF ---
    rows_tax = []
    for raw_name in tqdm(hosts_sorted, desc="Resolving host names", unit="name"):
        norm = raw_name.strip()
        lower = norm.lower()

        if is_sample_environment_only(lower):
            rows_tax.append({
                "input_name": raw_name,
                "resolved_name": "",
                "resolved_rank": "",
                "flag": "sample_descriptor",
                "kingdom": "",
                "phylum": "",
                "class": "",
                "order": "",
                "family": "",
                "genus": "",
                "species": "",
            })
            continue

        human_name = map_to_human_if_applicable(lower)
        if human_name:
            sci_name = human_name
        else:
            if lower in COMMON_NAME_MAP:
                sci_name = COMMON_NAME_MAP[lower]
            else:
                sci_name = gnverifier_canonical(norm)
                if not sci_name:
                    rows_tax.append({
                        "input_name": raw_name,
                        "resolved_name": "",
                        "resolved_rank": "",
                        "flag": "cannot_interpret",
                        "kingdom": "",
                        "phylum": "",
                        "class": "",
                        "order": "",
                        "family": "",
                        "genus": "",
                        "species": "",
                    })
                    continue

        gbif_info = gbif_taxonomy(sci_name)
        if gbif_info is None:
            rows_tax.append({
                "input_name": raw_name,
                "resolved_name": "",
                "resolved_rank": "",
                "flag": "cannot_interpret",
                "kingdom": "",
                "phylum": "",
                "class": "",
                "order": "",
                "family": "",
                "genus": "",
                "species": "",
            })
            continue

        family = gbif_info.get("family", "") or ""
        flag = "ok" if family else "no_family"

        rows_tax.append({
            "input_name": raw_name,
            "resolved_name": gbif_info.get("scientificName", "") or "",
            "resolved_rank": gbif_info.get("rank", "") or "",
            "flag": flag,
            "kingdom": gbif_info.get("kingdom", "") or "",
            "phylum": gbif_info.get("phylum", "") or "",
            "class": gbif_info.get("class", "") or "",
            "order": gbif_info.get("order", "") or "",
            "family": family,
            "genus": gbif_info.get("genus", "") or "",
            "species": gbif_info.get("species", "") or "",
        })

        time.sleep(0.05)

    df_tax = pd.DataFrame(rows_tax)
    df_tax.to_csv(RESOLVED_HOSTS_CSV, index=False)
    print(f"Saved taxonomy table to {RESOLVED_HOSTS_CSV}")
    print(df_tax["flag"].value_counts())


if __name__ == "__main__":
    main()
