#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 11:21:11 2026

@author: jaanajurvansuu

End-to-end host inference and taxonomy resolution pipeline.

Steps:
1) Read detection tables, pool accessions/labels, infer most probable host using NCBI (Entrez).
2) Extract unique host names.
3) Resolve host names, including manually curated higher taxa, through NCBI Taxonomy.

Requires:
  - pandas
  - biopython (Bio.Entrez)
  - requests
  - tqdm
"""
import re
import time
import collections
from pathlib import Path
from functools import lru_cache
from typing import List, Tuple, Optional, Dict

import pandas as pd
import requests
from tqdm import tqdm
from Bio import Entrez

# =========================
# USER SETTINGS
# =========================

# Input virus detection tables (TSV)
DETECTIONS_00 = Path.home() / "Desktop" / "detected_viruses_combined_segments.tsv"

# Entrez credentials (NCBI)
ENTREZ_EMAIL = X
ENTREZ_API_KEY = Z
NCBI_DELAY = 0.34                        # seconds between NCBI calls
RETMAX_SEARCH = 5                        # max nuccore records in label search
# Label searches are slow and may retrieve unrelated virus records.  Leave off
# unless you specifically need to recover hosts absent from accession records.
ENABLE_LABEL_SEARCH = False

# Output files
POOLED_HOSTS_TSV = Path.home() / "Desktop" / "most_probable_hosts_pooled1.tsv"
POOL_MEMBERS_TSV = Path.home() / "Desktop" / "pool_members1.tsv"
UNIQUE_HOSTS_TXT = Path.home() / "Desktop" / "unique_most_probable_hosts1.txt"
RESOLVED_HOSTS_CSV = Path.home() / "Desktop" / "resolved_hosts_full_taxonomy_final1.csv"
FINAL_OUTPUT_TSV = Path.home() / "Desktop" / "detected_viruses_combined_segments_with_hosts.tsv"

# Column name hints (script auto-detects first present)
ACC_COLS  = ["accession", "Accession", "nuccore", "genbank_accession"]
NAME_COLS = ["Name", "species", "label", "description", "Virus"]
ASM_COLS  = ["Assembly", "assembly", "set", "Set"]

# External APIs
GNVERIFIER_URL = "https://verifier.globalnames.org/api/v1/verifications"
GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"

# Accession-specific decisions from manual review.  The taxon can be a species
# or any higher taxon recognised by NCBI Taxonomy.
MANUAL_HOSTS_BY_TAXON = {
    "Chiroptera": "JF938078.1 JF938080.1 MT734814.1 MW597741.1 NC_032426.1 NC_034547.1 NC_038393.1 NC_038394.1 NC_038395.1 NC_077079.1 NC_077082.1 OK032372.1 JQ412755.1 JQ412756.1 JQ412757.1 JQ412758.1 JQ412759.1 JQ412760.1 JQ412761.1 JQ412762.1 JQ412763.1 JQ412764.1",
    "Oomycetes": "NC_022089.1",
    "Gallus gallus": "JQ308210.1" "JQ308213.1" "JQ308214.1" "NC_017091.1", 
    "Cynodon": "NC_076284.1", 
    "Anas": "MH453841.1 MH453848.1", 
    "Erinaceus": "PP421220.1", 
    "Sida": "MH939366.1 MH939411.1", 
    "Aves": "KT894018.1", 
    "Galleria mellonella": "NC_004286.1", 
    "Casphalia extranea": "NC_004288.1", 
    "Blattella germanica": "NC_005041.2", 
    "Mythimna loreyi": "NC_005341.1", 
    "Culex": "NC_009942.1 OQ190161.1 OQ190162.1 OQ190163.1 OQ190164.1 OQ190165.1 OQ190166.1 OQ190167.1 OQ190168.1 OQ190169.1 OQ190170.1 OQ190171.1 OQ190172.1", 
    "Culex pipiens": "NC_024503.1", 
    "Magallana gigas": "NC_032115.1", 
    "Odonata": "MG779479.1", 
    "Chloebia gouldiae": "NC_008522.1", 
    "Larus argentatus": "NC_008521.1", 
    "Paguroidea": "NC_027789.1 NC_027790.1", 
    "Macaca": "MH539648.1", 
    "Macaca fascicularis": "KY120185.1 KY120186.1 KY120190.1 KY120191.1 KY751022.1 KY751024.1", 
    "Culicidae": "NC_015874.1 NC_031687.1 NC_031688.1 NC_038481.1 MT096521.1 MT096523.1", 
    "Nycticorax nycticorax": "NC_016994.1", 
    "Psittacus erithacus": "MH188863.1", 
    "Phasianus colchicus": "KT954542.1 KT954543.1 KT954544.1 KT954545.1 KT954546.1 KT954547.1 KT954548.1 KT954549.1 KT954550.1 KT954551.1", 
    "Primates": "KT581587.1", 
    "Rodentia": "KX156156.1 MK530648.1 MF139053.1 MH892425.1 MH976711.1 NC_055465.1", 
    "Rattus norvegicus":"KT724051.1 KT724051.1",
    "Myodes rufocanus": "KY432931.1 NC_075977.1", 
    "Passer": "MG737712.1 MG812376.1 MG812377.1 MG812378.1 NC_016992.1", 
    "Ixodida": "MG972981.1 ON456881.1 ON456882.1 ON456885.1 ON456888.1 ON456889.1 ON456890.1 ON456894.1", 
    "Bos taurus": "NC_006259.1 NC_001442.1 MW655626.1", 
    "Zosteropidae": "NC_016991.1", 
    "Charadriiformes": "EU030978.1",
    "Bovidae": "OR682617.1",
    "Capra hircus":"NC_005336.1",
    "Macaca mulatta":"NC_012532.1",
    "Canis lupus familiaris":"KU871312.1 KU871313.1 NC_006564.1",
    "Oryctolagus cuniculus":"MT150090.1",
    "Ailuropoda melanoleuca":"NC_075330.1 NC_075331.1 NC_075332.1 NC_075333.1 NC_075334.1 NC_075335.1 NC_075336.1",
    "Homo sapiens": "KJ427799.1 NC_040640.1 DQ358805.1 DQ358806.1 DQ358807.1 DQ358811.1 EF554093.1 EF554094.1 EF554097.1 EF554101.1 EF554104.1 EF554105.1 EF554106.1 EF554107.1 EF554108.1 EF554109.1 EF554110.1 EF554111.1 EF554112.1 EF554113.1 EF554114.1 EF554126.1 EF554127.1 EF554128.1 EF554129.1 EF554130.1 EF554131.1 EF554132.1 EF554133.1 EF554134.1 EF554135.1 EF554136.1 AY167123.1 DQ854714.1 KP036481.1 AB369688.1 KF779353.1 KT447526.1 KT447527.1 NC_076138.1 AB433629.1 DQ028632.1 DQ028633.1 FM178558.1 KF541281.1 KF541282.1 KF541283.1 KF541284.1 KF541285.1 KF541286.1 KF541287.1 KF541288.1 KF541289.1 KF541290.1 KF541291.1 KT879915.1 KT879917.1 KT879919.1 KT879920.1 KT879921.1 KT879923.1 KT879925.1 KT879927.1 KT879928.1 KT879929.1 KT879930.1 KT879931.1 NC_011400.1 NC_039476.1 OQ744533.1 MG148341.1",
}

ACCESSION_HOST_MAP = {acc: taxon for taxon, accs in MANUAL_HOSTS_BY_TAXON.items() for acc in accs.split()}

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


def split_accessions(value) -> List[str]:
    """Return individual valid accessions from one cell, including segment lists."""
    parts = re.split(r"[;,|]", to_text(value))
    return [part.strip() for part in parts if ACCN_RE.fullmatch(part.strip())]

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


@lru_cache(maxsize=None)
def taxonomy_host_taxid(name: str) -> Optional[str]:
    if not name:
        return None
    try:
        r = Entrez.read(Entrez.esearch(
            db="taxonomy", term=f'"{name}"[All Names]', retmax=1
        ))
        pause()
        ids = r.get("IdList", [])
        return ids[0] if ids else None
    except Exception as exc:
        print(f"NCBI taxonomy search failed for {name!r}: {exc}")
        return None


HOST_RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies")


@lru_cache(maxsize=None)
def ncbi_host_taxonomy(name: str) -> dict:
    """Return NCBI taxonomy, retaining the supplied taxon's actual rank."""
    empty = {"host_taxid": "", "host_scientific_name": "", "host_rank": "",
             **{f"host_{rank}": "" for rank in HOST_RANKS}}
    taxid = taxonomy_host_taxid(name)
    if not taxid:
        return empty
    try:
        record = Entrez.read(Entrez.efetch(db="taxonomy", id=taxid, retmode="xml"))[0]
        pause()
    except Exception:
        return empty
    result = empty | {"host_taxid": str(record.get("TaxId", taxid)),
                      "host_scientific_name": str(record.get("ScientificName", "")),
                      "host_rank": str(record.get("Rank", "")).lower()}
    nodes = list(record.get("LineageEx", [])) + [record]
    for node in nodes:
        rank = str(node.get("Rank", "")).lower()
        if rank in HOST_RANKS:
            result[f"host_{rank}"] = str(node.get("ScientificName", ""))
    return result


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
        d["accessions"].update(split_accessions(acc))
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
    manual_votes = [ACCESSION_HOST_MAP[acc] for acc in accs if acc in ACCESSION_HOST_MAP]
    if manual_votes:
        host, conf = majority_vote(manual_votes)
        return host, "", conf, "T0_manual_accession", 0, manual_votes
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
        return host, "", conf, "T1_accession_pool", used, host_votes

    # Tier 2: label-based
    label = next(iter(ginfo["labels"]), "")
    if ENABLE_LABEL_SEARCH and label:
        ids = esearch_nuccore_ids_by_label(label, retmax=RETMAX_SEARCH)
        for uid in ids:
            gb = efetch_nuccore_xml(uid)
            host_votes.extend(extract_hosts_from_gbset(gb))
        if host_votes:
            host, conf = majority_vote(host_votes)
            return host, "", conf, "T2_label_search", used, host_votes

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
    "hbsag carrier", "elite suppressor", "commercial sex worker",
    "neonate", "immunodeficient individual", "54 year old female","vaccine isolate",
    "tracheal aspirate", "Peripheral blood mononuclear cell", "IRF Isolate",
    "peripheral blood mononuclear cells", "cervicovaginal cells",
]


def map_to_human_if_applicable(name_lower: str) -> Optional[str]:
    for hint in HUMAN_HINTS:
        if hint.lower() in name_lower:
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
    "healthy cattle blood": "Bos taurus",
    "holstein dairy cattle": "Bos taurus",
    "pi calf": "Bos taurus",
    "steer": "Bos taurus",
    "supermarket ground beef": "Bos taurus",
    "three week old calf with oral lesions, then passaged in primary lamb kidney cells": "Bos taurus",
    "beef": "Bos taurus",
    "buffalo": "Bubalus bubalis",
    "bovine adenovirus 2 culture atcc vr-314": "Bos taurus",

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
    "42-day-old pig": "Sus scrofa domesticus",
    "asymptomatic piglets": "Sus scrofa domesticus",
    "domestic pig (sus scrofa domestica)": "Sus scrofa domesticus",
    "domestic pig serum": "Sus scrofa domesticus",
    "fecal sample of experimentally infected swine": "Sus scrofa domesticus",
    "sucking piglets": "Sus scrofa domesticus",
    "swine feces": "Sus scrofa domesticus",
    "tibetan pig": "Sus scrofa domesticus",
    "tibetan pigs": "Sus scrofa domesticus",
    "wild boar": "Sus scrofa",
    "wild boar liver": "Sus scrofa",

    # dogs and other canids
    "dog": "Canis lupus familiaris",
    "dog ": "Canis lupus familiaris",
    "mongrel dog": "Canis lupus familiaris",
    "domestic dog": "Canis lupus familiaris",
    "dog (rottweiler mix adult)": "Canis lupus familiaris",
    "dog (breed island mix)": "Canis lupus familiaris",
    "canine": "Canis lupus familiaris",
    "canine lupus familaris": "Canis lupus familiaris",
    "african wild dog": "Lycaon pictus",
    "raccoon dog": "Nyctereutes procyonoides",
    "side-striped jackal": "Lupulella adusta",
    "wolf": "Canis lupus",
    "fox": "Vulpes vulpes",

    # cats and other felids
    "cat": "Felis catus",
    "cat ": "Felis catus",
    "domestic cat": "Felis catus",
    "domestic cat with diarrhea": "Felis catus",
    "domestic cat; sex: male": "Felis catus",
    "feline": "Felis catus",
    "jejunum of a kitten with feline infectious peritonitis": "Felis catus",
    "amur tiger": "Panthera tigris altaica",
    "lion": "Panthera leo",
    "leopard": "Panthera pardus",

    # goats, sheep, and equids
    "goat": "Capra hircus",
    "black goat": "Capra hircus",
    "goat kid": "Capra hircus",
    "goat; asymptomatic animal": "Capra hircus",
    "yunling goat": "Capra hircus",
    "lamb": "Ovis aries",
    "domestic sheep": "Ovis aries",
    "ovine": "Ovis aries",
    "donkey": "Equus africanus asinus",
    "foal (yili horse)": "Equus caballus",
    "equine": "Equus caballus",
    "horse": "Equus caballus",
    "horse cutaneous papilloma": "Equus caballus",

    # other mammals
    "african green monkey": "Chlorocebus aethiops",
    "alpaca": "Vicugna pacos",
    "american white-tailed deer fibroma": "Odocoileus virginianus",
    "baikal seal": "Pusa sibirica",
    "chimpanzee": "Pan troglodytes",
    "c57bl/6 mouse": "Mus musculus",
    "ferret": "Mustela putorius furo",
    "japanese macaque": "Macaca fuscata",
    "lapine": "Oryctolagus cuniculus",
    "masked palm civet": "Paguma larvata",
    "mouse": "Mus musculus",
    "mouse strain km": "Mus musculus",
    "norway rat (rattus norvegicus)": "Rattus norvegicus",
    "northern red-backed vole": "Myodes rutilus",
    "rabbit": "Oryctolagus cuniculus",
    "rabbit kidney cells rk-13; atcc:ccl-37": "Oryctolagus cuniculus",
    "rat": "Rattus norvegicus",
    "rhesus macaque": "Macaca mulatta",
    "rhesus monkey": "Macaca mulatta",
    "striped skunk": "Mephitis mephitis",
    "tissues obtained from a naturally infected mouse": "Mus musculus",
    "water buffalo": "Bubalus bubalis",
    "white-tail deer": "Odocoileus virginianus",
    "yak": "Bos grunniens",
    "mink": "Neogale vison",
    "bat colony of molossus molossus and tadarida brasiliensis": "Chiroptera",
    "feces of an interferon-gamma receptor knock-out mouse": "Mus musculus",
    "rag1 knock-out mouse": "Mus musculus",
    "spleen of nod-scid": "Mus musculus",
    "roe deer": "Capreolus",

    # birds
    "chicken": "Gallus gallus domesticus",
    "Gallus gallus":"Gallus gallus domesticus",
    "poultry":"Gallus gallus domesticus",
    "broiler chicken": "Gallus gallus domesticus",
    "broiler": "Gallus gallus domesticus",
    "broiler chicken voucher m05833_0083": "Gallus gallus domesticus",
    "local chicken": "Gallus gallus domesticus",
    "layer chicken": "Gallus gallus domesticus",
    "chicken meat": "Gallus gallus domesticus",
    "silkie": "Gallus gallus domesticus",
    "chicken breed luqin": "Gallus gallus domesticus",
    "chicken; breed: yellow feather broiler": "Gallus gallus domesticus",
    "chicken (broiler breeder)": "Gallus gallus domesticus",
    "chicken; brabanter": "Gallus gallus domesticus",
    "fowl": "Gallus gallus domesticus",
    "duck": "Anas platyrhynchos",
    "cherry valley breeder duck": "Anas platyrhynchos",
    "mallard": "Anas platyrhynchos",
    "pekin duck": "Anas platyrhynchos",
    "muscovy duck": "Cairina moschata",
    "green-winged teal": "Anas carolinensis",
    "northern shoveler": "Spatula clypeata",
    "turkey": "Meleagris gallopavo",
    "turkey small intestine": "Meleagris gallopavo",
    "pigeon": "Columba livia",
    "pigeon ": "Columba livia",
    "carrier pigeon": "Columba livia",
    "fancy pigeon": "Columba livia",
    "feral pigeon": "Columba livia",
    "racing pigeon": "Columba livia",
    "rock dove": "Columba livia",
    "agapomis roseicollis": "Agapornis roseicollis",
    "alexandrine parakeet": "Psittacula eupatria",
    "black-tailed gull": "Larus crassirostris",
    "budgerigar": "Melopsittacus undulatus",
    "common magpie": "Pica pica",
    "eurasian bullfinch (pyrrhula pyrrhula griseiventris)": "Pyrrhula pyrrhula",
    "eurasian jackdaw (corvus monedula)": "Corvus monedula",
    "european starling (sturnus vulgaris)": "Sturnus vulgaris",
    "grey-backed thrush": "Turdus hortulorum",
    "house crow": "Corvus splendens",
    "red-necked avocet": "Recurvirostra novaehollandiae",
    "red-whiskered bulbul": "Pycnonotus jocosus",
    "serinus canaria": "Serinus canaria",
    "spotted dove": "Spilopelia chinensis",
    "white-rumped munia": "Lonchura striata",
    "goose": "Anser anser",
    "indian parakeet": "Psittacula krameri",
    "mule duck": "Anatidae",
    "magpie-robin":"Copsychus",
    
    #reptiles
    "corn snake (elaphe guttata)": "Pantherophis guttatus",

    # fish
    "rainbow trout tissue": "Oncorhynchus mykiss",
    "atlantic herring": "Clupea harengus",
    "atlantic salmon": "Salmo salar",
    "farmed atlantic salmon": "Salmo salar",
    "haddock tissue": "Melanogrammus aeglefinus",
    "largemouth bass": "Micropterus salmoides",

    # invertebrates
    "house cricket": "Acheta domesticus",
    "black fly (simulium ubiquitum)": "Simulium ubiquitum",
    "drosophila suzukii; adult": "Drosophila suzukii",
    "eastern lubber grasshopper (romulea microptera)": "Romalea microptera",
    "honey bee": "Apis mellifera",
    "solenopsis invicta; worker ants": "Solenopsis invicta",
    "spider (cybaeus signifer)": "Cybaeus signifer",
    "spider (segestria pacifica)": "Segestria pacifica",
    "tropiladaps clareae": "Tropilaelaps clareae",
    "bee larvae": "Apis mellifera",
    "mud crab": "Scylla serrata",
    "palmetto bug":"Eurycotis floridana",
    "cricket": "Orthoptera",
    "feces from pool of neophema elegans and neophema splendida": "Neophema",
    "large pool of wild-caught drosophila (mixed d. ananassae, d. melanogaster, d. malerkotliana, and scaptodrosophila latifasciaeformis)": "Drosophilidae",
    
    # plants
    "barley": "Hordeum vulgare",
    "wheat": "Triticum aestivum",
    "winter wheat": "Triticum aestivum",
    "spring wheat": "Triticum aestivum",
    "maize": "Zea mays",
    "cucumber": "Cucumis sativus",
    "muskmelon": "Cucumis melo",
    "melon": "Cucumis melo",
    "snake melon": "Cucumis melo",
    "pea": "Pisum sativum",
    "pepper": "Capsicum annuum",
    "bell pepper": "Capsicum annuum",
    "grapevine": "Vitis vinifera",
    "grapevine cultivar pinot gris": "Vitis vinifera",
    "grapevine 'pinot gris'": "Vitis vinifera",
    "groundnut": "Arachis hypogaea",
    "strawberry": "Fragaria × ananassa",
    "sugar beet": "Beta vulgaris",
    "common reed": "Phragmites australis",
    "red currant cv. holandsky cerveny": "Ribes rubrum",
    "tomato; cultivar=candela": "Solanum lycopersicum",
    "tomato": "Solanum lycopersicum",
    "eggplant": "Solanum melongena",
    "black nightshade": "Solanum nigrum",
    "tobacco": "Nicotiana tabacum",
    "soybean": "Glycine max",
    "chickpea": "Cicer arietinum",
    "broad bean": "Vicia faba",
    "faba bean": "Vicia faba",
    "radish": "Raphanus sativus",
    "turnip": "Brassica rapa",
    "broccoli (brassica oleracea var. italica)": "Brassica oleracea",
    "cauliflower": "Brassica oleracea",
    "oilseed rape": "Brassica napus",
    "watermelon": "Citrullus lanatus",
    "bitter apple": "Citrullus colocynthis",
    "bottle gourd": "Lagenaria siceraria",
    "sweet cherry": "Prunus avium",
    "white clover": "Trifolium repens",
    "saguaro": "Carnegiea gigantea",
    "cynodon dactylon l. (bermuda grass)": "Cynodon dactylon",
    "sorghum sudangrass hybrid (sorghum bicolor (l.) moench x sorghum bicolor var. sudanense)": "Sorghum bicolor",
    "garlic": "Allium sativum",
    "pumpkin sp.": "Cucurbita pepo",
    "squash": "Cucurbita pepo",
    "zucchini": "Cucurbita pepo",
    "cotton": "Gossypium hirsutum",
    "(citrus unshiu x citrus sinensis) x citrus reticulata cv. setoka": "Citrus",
    "cherry": "Prunus",
    "grape": "Vitis",
    "pumpkin": "Cucurbita",
    "raspberry": "Rubus",
    "vitis labrusca x vitis vinifera 'shine muscat'": "Vitis",
    "vitis labruscana x vitis vinifera cv shine muscat": "Vitis",
    "vitis vinifera x vitis berlandieri cv. 41b": "Vitis",
    "walnut": "Juglans",
    "marigold":"Tagetes",
}

SAMPLE_ENV_KEYWORDS = [
    "sewage", "wastewater", "influent", "oxidation pond",
    "sewage contaminated water", "sewage oxidation pond",
    "secondary treated and disinfected wastewater effluent",
    "benthic sediment", "freshwater lake", "ocean water",
    "irrigation water", "environment - desk top", "nest material",
    "pooled weed plants","leaves","pooled weed plants gathered at the immediate surroundings of tomato farms",
    "macrophytes collected in freshwater springs",
]

def is_sample_environment_only(name_lower: str) -> bool:
    for kw in SAMPLE_ENV_KEYWORDS:
        if kw in name_lower:
            return True
    if name_lower.strip() in {"stool", "feces", "faeces", "blood"}:
        return True
    return False


def resolve_host_name(raw_name: str) -> Tuple[str, str]:
    """Choose a taxonomic lookup name while preserving all resolution paths.

    Returns ``(lookup_name, status)``. An empty lookup name means that the
    label is an environmental/sample descriptor and must not be sent to NCBI.
    """
    norm = to_text(raw_name)
    lower = norm.lower()
    human_name = map_to_human_if_applicable(lower)
    if human_name:
        return human_name, "ok"
    if lower in COMMON_NAME_MAP:
        return COMMON_NAME_MAP[lower], "ok"
    if is_sample_environment_only(lower):
        return "", "sample_descriptor"
    # This was present in the original pipeline and is essential for labels
    # containing common names, strain details, parentheses, or cultivars.
    return gnverifier_canonical(norm) or norm, "ok"


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


def apply_gbif_taxonomy_fallback(ncbi_taxonomy: dict, lookup_name: str) -> Tuple[dict, str]:
    """Fill only NCBI-missing taxonomy levels from GBIF.

    NCBI remains the authoritative source for ``host_taxid`` and ``host_rank``
    whenever it resolves a taxon. GBIF is queried only to complete blank rank
    fields, which is important for lineages such as Oomycota that lack an NCBI
    node explicitly ranked ``kingdom``.
    """
    tax = dict(ncbi_taxonomy)
    gbif_lookup_name = tax.get("host_scientific_name") or lookup_name
    gbif_info = gbif_taxonomy(gbif_lookup_name)

    if not gbif_info:
        return tax, "NCBI" if tax.get("host_taxid") else "none"

    filled_from_gbif = False
    for rank in HOST_RANKS:
        gbif_value = gbif_info.get(rank, "") or ""
        host_column = f"host_{rank}"
        if not tax.get(host_column) and gbif_value:
            tax[host_column] = str(gbif_value)
            filled_from_gbif = True

    # Use GBIF's name/rank only if NCBI did not resolve the taxon at all.
    if not tax.get("host_scientific_name"):
        tax["host_scientific_name"] = str(gbif_info.get("scientificName", "") or "")
    if not tax.get("host_rank"):
        tax["host_rank"] = str(gbif_info.get("rank", "") or "").lower()

    if tax.get("host_taxid"):
        source = "NCBI+GBIF" if filled_from_gbif else "NCBI"
    else:
        source = "GBIF"
    return tax, source

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

    # --- 2) Extract unique host names and resolve through NCBI Taxonomy ---
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

    # NCBI supports species and higher ranks (e.g. Chiroptera, Culex).
    rows_tax = []
    taxonomy_by_lookup_name = {}
    for raw_name in tqdm(hosts_sorted, desc="Resolving host names", unit="name"):
        lookup_name, status = resolve_host_name(raw_name)
        if lookup_name:
            if lookup_name not in taxonomy_by_lookup_name:
                ncbi_tax = ncbi_host_taxonomy(lookup_name)
                taxonomy_by_lookup_name[lookup_name] = apply_gbif_taxonomy_fallback(
                    ncbi_tax, lookup_name
                )
            tax, taxonomy_source = taxonomy_by_lookup_name[lookup_name]
        else:
            tax = {"host_taxid": "", "host_scientific_name": "", "host_rank": "",
                   **{f"host_{rank}": "" for rank in HOST_RANKS}}
            taxonomy_source = "environmental_source"
        flag = "ok" if taxonomy_source in {"NCBI", "NCBI+GBIF", "GBIF"} else \
               status if status != "ok" else "cannot_interpret"
        rows_tax.append({"input_name": raw_name, "name_used_for_lookup": lookup_name,
                         "flag": flag, "taxonomy_source": taxonomy_source, **tax})

    df_tax = pd.DataFrame(rows_tax)
    
    # Label environmental/sample descriptors consistently in host fields.
    environmental_mask = df_tax["flag"].eq("sample_descriptor")

    environmental_columns = [
        "host_scientific_name",
        "host_kingdom",
        "host_phylum",
        "host_class",
        "host_order",
        "host_family",
        "host_genus",
        "host_species",
        "host_subspecies",
        ]

    df_tax.loc[environmental_mask, environmental_columns] = "Environmental source"
    df_tax.loc[environmental_mask, "host_rank"] = "environmental_source"
    df_tax.loc[environmental_mask, "host_taxid"] = ""
    df_tax.loc[environmental_mask, "taxonomy_source"] = "environmental_source"
  
    df_tax.to_csv(RESOLVED_HOSTS_CSV, index=False)
    print(f"Saved taxonomy table to {RESOLVED_HOSTS_CSV}")
    print(df_tax["flag"].value_counts())

    # --- 3) Add host assignment and taxonomy to every combined-segment row ---
    taxonomy_for_merge = df_tax.rename(columns={"input_name": "most_probable_host"})
    pooled_for_merge = pooled.drop(columns=["host_taxid"], errors="ignore").merge(
        taxonomy_for_merge, on="most_probable_host", how="left"
    )
    host_columns = ["group_id", "most_probable_host", "host_taxid", "confidence",
                    "evidence_tier", "all_hosts_seen", "host_scientific_name", "host_rank",
                    "taxonomy_source"] + \
                   [f"host_{rank}" for rank in HOST_RANKS]
    detections["group_id"] = detections.apply(
        lambda row: make_group_id(row, acc_col, asm_col, name_col), axis=1
    )
    final = detections.merge(pooled_for_merge[host_columns], on="group_id", how="left")
    final.to_csv(FINAL_OUTPUT_TSV, sep="\t", index=False)
    print(f"Saved combined detections with host taxonomy to {FINAL_OUTPUT_TSV}")


if __name__ == "__main__":
    main()
