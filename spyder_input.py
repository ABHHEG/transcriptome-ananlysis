"""File/schema configuration and organism detection for the Spyder pipeline."""
from pathlib import Path
import re

import numpy as np
import pandas as pd
import requests

KEGG_ORGANISM_LIST_URL = "https://rest.kegg.jp/list/organism"


def _read_table(path):
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, dtype=str)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype=str)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, dtype=str)
    raise ValueError("Supported input files: CSV, TSV/TXT, XLSX, or XLS.")


def _normalized_name(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def _choose(items, prompt, st, default=0):
    if not items:
        raise ValueError("No columns are available for that selection.")
    return st.selectbox(prompt, options=items, index=default)


def _likely_gene_column(columns):
    aliases = {
        "geneid", "geneidentifier", "gene", "locus", "locustag",
        "featureid", "systematicname", "orf", "accession", "identifier",
        "transcriptid", "ensemblgeneid", "ensemblid", "id"
    }
    matches = [col for col in columns if _normalized_name(col) in aliases]
    return matches


def _numeric_expression_candidates(frame, gene_column):
    excluded = (
        "pvalue", "padj", "adjustedp", "fdr", "qvalue", "log2fc",
        "foldchange", "chromosome", "position", "start", "end", "length",
        "replicate", "batch", "taxid", "samplecount", "numreads", "timepoint", "timehr"
    )
    found = []
    minimum_valid = max(2, int(np.ceil(len(frame) * 0.5)))
    for column in frame.columns:
        if column == gene_column:
            continue
        name = _normalized_name(column)
        if any(token in name for token in excluded):
            continue
        values = pd.to_numeric(frame[column], errors="coerce")
        if int(values.notna().sum()) >= minimum_valid:
            found.append(column)
    return found


def _suggested_expression_columns(columns):
    tokens = ("tpm", "fpkm", "rpkm", "cpm", "expression", "abundance", "count")
    suggested = [col for col in columns if any(token in _normalized_name(col) for token in tokens)]
    return suggested if len(suggested) >= 2 else columns


def _pick_columns(columns, prompt, default_columns):
    print(f"\n{prompt}")
    for index, column in enumerate(columns, start=1):
        marker = " [suggested]" if column in default_columns else ""
        print(f"  {index}. {column}{marker}")
    default_text = ",".join(str(columns.index(col) + 1) for col in default_columns)
    raw = input(f"Enter comma-separated column numbers [default {default_text}]: ").strip()
    if not raw:
        selected = list(default_columns)
    else:
        try:
            indices = [int(piece.strip()) for piece in raw.split(",") if piece.strip()]
            if not indices or any(index < 1 or index > len(columns) for index in indices):
                raise ValueError
            selected = [columns[index - 1] for index in indices]
        except ValueError:
            raise ValueError("Enter valid numeric column numbers separated by commas.")
    if len(selected) < 2:
        raise ValueError("Select at least two expression columns/time points.")
    return selected


def _infer_time_label(column):
    text = str(column)
    match = re.search(r"(?i)(\d+(?:\.\d+)?)\s*(hours?|hrs?|h|days?|d|minutes?|mins?|min|weeks?|w)\b", text)
    if match:
        value = match.group(1)
        unit = match.group(2).lower()
        if unit.startswith("d"):
            suffix = "day"
        elif unit.startswith("m"):
            suffix = "min"
        elif unit.startswith("w"):
            suffix = "week"
        else:
            suffix = "hr"
        return f"{value} {suffix}"
    match = re.search(r"(?i)(?:^|[_\-\s])(?:timepoint|time|tp|t)[_\-\s]*(\d+(?:\.\d+)?)(?:$|[_\-\s])", text)
    if match:
        return match.group(1)
    return ""


def _time_sort_key(label):
    match = re.match(r"\s*(-?\d+(?:\.\d+)?)", str(label))
    return (0, float(match.group(1))) if match else (1, str(label).lower())


def _ordered_unique(values):
    unique = list(dict.fromkeys(str(value).strip() for value in values if str(value).strip()))
    if unique and all(_time_sort_key(value)[0] == 0 for value in unique):
        unique.sort(key=_time_sort_key)
    return unique


def _expression_unit(columns, st):
    text = " ".join(str(col).lower() for col in columns)
    if "tpm" in text:
        suggestion = "TPM"
    elif "fpkm" in text:
        suggestion = "FPKM"
    elif "rpkm" in text:
        suggestion = "RPKM"
    elif "cpm" in text:
        suggestion = "CPM"
    elif "count" in text:
        suggestion = "Counts"
    else:
        suggestion = "Normalized expression / other"
    options = ["TPM", "FPKM", "RPKM", "CPM", "Counts", "Normalized expression / other"]
    return _choose(options, f"Expression value type (suggested: {suggestion})", st, options.index(suggestion))


def configure_input(st, default_path):
    """Read and normalize wide or long expression input into a gene-by-time matrix."""
    default_path = Path(default_path).resolve()
    supplied = input(f"Input expression file path [{default_path}]: ").strip()
    path = Path(supplied).expanduser() if supplied else default_path
    if not path.is_absolute():
        path = Path.cwd() / path
    path = path.resolve()
    if not path.exists():
        raise FileNotFoundError(f"Input file was not found: {path}")
    frame = _read_table(path)
    if frame.empty:
        raise ValueError("The selected input file has no data rows.")

    gene_matches = _likely_gene_column(frame.columns)
    gene_default = gene_matches[0] if gene_matches else frame.columns[0]
    gene_column = _choose(list(frame.columns), "Select the gene identifier column", st, list(frame.columns).index(gene_default))
    layout_options = [
        "Wide: one row per gene, with a separate expression column per time point/sample",
        "Long: one row per gene/time point measurement"
    ]
    layout = _choose(layout_options, "Select the input layout", st, 0)
    expression_candidates = _numeric_expression_candidates(frame, gene_column)
    if not expression_candidates:
        raise ValueError("No numeric expression columns were detected. Check the selected gene ID and input file.")

    if layout.startswith("Wide:"):
        suggestions = _suggested_expression_columns(expression_candidates)
        chosen_columns = _pick_columns(expression_candidates, "Select the expression columns", suggestions)
        labels_by_column = {}
        for column in chosen_columns:
            guessed = _infer_time_label(column)
            suffix = f" [{guessed}]" if guessed else ""
            entered = input(f"Time point label for '{column}'{suffix}: ").strip()
            label = entered or guessed
            if not label:
                raise ValueError(f"Please provide a time point label for column '{column}'.")
            labels_by_column[column] = label

        work = pd.DataFrame({"Gene_ID": frame[gene_column].astype(str).str.strip()})
        for index, column in enumerate(chosen_columns):
            work[f"source::{index}"] = pd.to_numeric(frame[column], errors="coerce")
        grouped_columns = {}
        for index, column in enumerate(chosen_columns):
            grouped_columns.setdefault(labels_by_column[column], []).append(f"source::{index}")
        # Use the actual chosen source columns to correctly aggregate replicate columns sharing a time label.
        values = pd.DataFrame({"Gene_ID": work["Gene_ID"]})
        for label, source_columns in grouped_columns.items():
            values[label] = work[source_columns].mean(axis=1, skipna=True)
        values = values.groupby("Gene_ID", as_index=False)[list(grouped_columns)].mean(numeric_only=True)
        labels = _ordered_unique(grouped_columns)
    else:
        time_candidates = [col for col in frame.columns if col != gene_column]
        time_column = _choose(time_candidates, "Select the time point / condition column", st, 0)
        expression_column = _choose(expression_candidates, "Select the expression value column", st, 0)
        work = pd.DataFrame({
            "Gene_ID": frame[gene_column].astype(str).str.strip(),
            "_time": frame[time_column].astype(str).str.strip(),
            "_value": pd.to_numeric(frame[expression_column], errors="coerce"),
        })
        work = work.dropna(subset=["_value"])
        if work.empty:
            raise ValueError("The selected long-format expression column has no numeric values.")
        values = work.pivot_table(index="Gene_ID", columns="_time", values="_value", aggfunc="mean").reset_index()
        labels = _ordered_unique([column for column in values.columns if column != "Gene_ID"])
        values = values[["Gene_ID"] + labels]

    if len(labels) < 2:
        raise ValueError("At least two distinct time points are required.")
    unit = _expression_unit(chosen_columns if layout.startswith("Wide:") else [expression_column], st)

    # Standardize names for downstream analysis without requiring fixed input headers.
    normalized = pd.DataFrame({"Gene_ID": values["Gene_ID"].astype(str).str.strip()})
    timepoint_columns = {}
    for index, label in enumerate(labels):
        standardized_column = f"Expression_{index + 1}"
        normalized[standardized_column] = pd.to_numeric(values[label], errors="coerce")
        timepoint_columns[label] = standardized_column
    normalized = normalized.dropna(subset=["Gene_ID"])
    normalized = normalized[normalized["Gene_ID"] != ""]
    normalized = normalized.groupby("Gene_ID", as_index=False)[list(timepoint_columns.values())].mean(numeric_only=True)
    normalized = normalized.dropna(subset=list(timepoint_columns.values()), how="all")
    if normalized.empty:
        raise ValueError("No usable gene expression rows remained after cleaning.")
    if (normalized[list(timepoint_columns.values())].fillna(0) < 0).any().any():
        raise ValueError("Expression values must be non-negative for log₂ fold-change analysis.")
    return normalized, timepoint_columns, unit, path, gene_column, frame


def _read_organism_records():
    try:
        response = requests.get(KEGG_ORGANISM_LIST_URL, timeout=20)
        response.raise_for_status()
    except requests.RequestException as exc:
        print(f"Could not retrieve KEGG organism list: {exc}")
        return []
    records = []
    for line in response.text.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        entry_id, code, name = fields[:3]
        if not re.fullmatch(r"[a-z0-9]{3,4}", code.strip().lower()):
            continue
        taxonomy_match = re.search(r"(?:TAX:|taxid[:= ]*)(\d+)", line, re.IGNORECASE)
        records.append({
            "entry_id": entry_id.strip(),
            "code": code.strip().lower(),
            "name": name.strip(),
            "taxid": taxonomy_match.group(1) if taxonomy_match else None,
        })
    return records


def _organism_hint(gene_ids, source_frame=None, gene_column=None):
    if source_frame is not None:
        for column in source_frame.columns:
            if _normalized_name(column) in {"organism", "species", "organismname", "speciesname", "scientificname"}:
                values = source_frame[column].dropna().astype(str).str.strip().unique()
                if len(values) == 1 and values[0]:
                    return values[0], "organism metadata column"
    values = [str(value).strip() for value in list(gene_ids)[:1000]]
    prefixes = []
    for value in values:
        match = re.match(r"^([A-Za-z0-9]{3,4}):", str(value).strip())
        if match:
            prefixes.append(match.group(1).lower())
    unique_prefixes = list(dict.fromkeys(prefixes))
    if (len(unique_prefixes) == 1 and values
            and len(prefixes) / len(values) >= 0.8):
        return unique_prefixes[0], "KEGG-style gene ID prefix"
    if values and sum(value.upper().startswith("PAS_") for value in values) / len(values) >= 0.8:
        return "ppa", "PAS_ locus ID pattern"
    ensembl_prefixes = {
        "ENSMUSG": "Mus musculus", "ENSDARG": "Danio rerio", "ENSGALG": "Gallus gallus",
        "ENSBTAG": "Bos taurus", "ENSRNOG": "Rattus norvegicus",
        "ENSPTRG": "Pan troglodytes", "ENSCAFG": "Canis lupus familiaris", "ENSG": "Homo sapiens",
    }
    for prefix, species in ensembl_prefixes.items():
        if values and sum(value.upper().startswith(prefix) for value in values) / len(values) >= 0.8:
            return species, f"{prefix} identifier pattern"
    return "", "no reliable organism signature in the gene IDs"


def configure_organism(st, gene_ids, source_frame=None, gene_column=None):
    """Infer a likely organism, let the user confirm/correct it, and resolve KEGG metadata."""
    records = _read_organism_records()
    hint, hint_source = _organism_hint(gene_ids, source_frame, gene_column)
    if hint:
        print(f"Likely organism from {hint_source}: {hint}")
    else:
        print(f"Organism could not be inferred safely: {hint_source}.")

    if source_frame is not None:
        for column in source_frame.columns:
            if _normalized_name(column) in {"taxid", "taxonomyid", "ncbitaxid"}:
                taxids = source_frame[column].dropna().astype(str).str.extract(r"(\d+)", expand=False).unique()
                if len(taxids) == 1:
                    metadata_taxid = str(taxids[0])
                    break
        else:
            metadata_taxid = None
    else:
        metadata_taxid = None

    prompt_default = hint or ""
    entered = input(f"Confirm organism or enter scientific name / KEGG code [{prompt_default or 'blank to skip KEGG'}]: ").strip()
    query = entered or hint
    if not query:
        print("No organism selected. Organism-specific KEGG mapping will be skipped.")
        return {"name": "Unspecified", "kegg_code": None, "taxid": metadata_taxid, "hint": hint_source}

    record = None
    query_lower = query.lower()
    for item in records:
        if query_lower == item["code"] or query_lower == item["entry_id"].lower():
            record = item
            break
    if record is None:
        matches = [item for item in records if query_lower == item["name"].lower()]
        if not matches:
            matches = [item for item in records if query_lower in item["name"].lower()]
        if len(matches) == 1:
            record = matches[0]
        elif len(matches) > 1:
            print("Multiple KEGG organisms match that name:")
            for item in matches[:20]:
                print(f"  {item['code']} — {item['name']}")
            followup = input("Enter the exact KEGG code to use (blank to skip KEGG): ").strip().lower()
            record = next((item for item in matches if item["code"] == followup), None)

    # An explicit code can be used even when the KEGG organism-list endpoint is unavailable.
    explicit_code = re.fullmatch(r"[a-zA-Z0-9]{3,4}", query)
    if record is None and explicit_code:
        record = {"code": query.lower(), "name": query, "entry_id": "", "taxid": None}
    if record is None:
        print(f"Could not resolve '{query}' to one KEGG organism. KEGG mapping will be skipped.")
        return {"name": query, "kegg_code": None, "taxid": metadata_taxid, "hint": hint_source}

    taxid = metadata_taxid or record.get("taxid")
    if not taxid and record.get("entry_id"):
        try:
            detail = requests.get(f"https://rest.kegg.jp/get/{record['entry_id']}", timeout=20)
            if detail.ok:
                match = re.search(r"Taxonomy\s+TAX:(\d+)", detail.text)
                if match:
                    taxid = match.group(1)
        except requests.RequestException:
            pass
    print(f"Using organism: {record['name']} | KEGG code: {record['code']} | Taxonomy ID: {taxid or 'not resolved'}")
    return {"name": record["name"], "kegg_code": record["code"], "taxid": taxid, "hint": hint_source}





