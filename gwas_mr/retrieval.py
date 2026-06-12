"""GWAS and eQTL data retrieval.

Refactored from gwas_eqtl_data_retrieval_pipeline.py into a callable
function ``retrieve_gwas_data()`` that accepts all parameters explicitly
(no interactive ``input()`` calls, no module-level side-effects).
"""

import gzip
import json
import os
import re
import shutil
import tempfile
import time
from difflib import get_close_matches
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import pandas as pd
import requests
from requests.exceptions import ChunkedEncodingError, ConnectionError, Timeout

# ── Constants ──────────────────────────────────────────────────────────────

CHUNK_SIZE = 1024 * 1024
HTTP_TIMEOUT = 30
HTTP_RETRIES = 5
RETRY_SLEEP = 3
MIN_HARMONISED_FILE_TYPES = 2

EURO_POS = ["european", "eur", "white british", "british", "irish", "ceu", "nfe"]
NON_EURO_NEG = ["asian", "african", "hispanic", "latino", "admixed"]

# Tiered dataset selection — evaluate multiple candidates and pick the
# one with the best eQTL instrument overlap rather than largest N alone.
TIERED_EVAL_TOP_K = 5
MIN_GWAS_N = 50_000
EQTL_P_THRESH_INSTRUMENT = 5e-8
EQTL_FSTAT_THRESH_INSTRUMENT = 10
OVERLAP_EVAL_LIMIT = 3  # max candidates to download for overlap evaluation

# ── Logging ────────────────────────────────────────────────────────────────


def log(msg: str) -> None:
    print(msg, flush=True)


# ── String helpers ─────────────────────────────────────────────────────────


def safe_name(s: str) -> str:
    """Normalise a string into a safe filesystem-friendly key."""
    return re.sub(r"[^\w]+", "_", str(s).lower())[:150].strip("_")


# ── HTTP helpers ───────────────────────────────────────────────────────────


def safe_request_get(
    url: str,
    stream: bool = False,
    timeout: int = HTTP_TIMEOUT,
    retries: int = HTTP_RETRIES,
) -> Optional[requests.Response]:
    """GET with retries and simple back-off.

    Returns ``None`` for permanent failures (404) instead of retrying.
    """
    last_err: Optional[Exception] = None
    for i in range(1, retries + 1):
        try:
            r = requests.get(url, stream=stream, timeout=timeout)
            if r.status_code == 404:
                log(f"[INFO] 404 Not Found — skipping: {url}")
                return None
            r.raise_for_status()
            return r
        except Exception as e:
            last_err = e
            log(f"[WARN] GET failed ({i}/{retries}) for {url} -> {e}")
            time.sleep(RETRY_SLEEP * i)
    raise RuntimeError(
        f"Failed GET after {retries} retries: {url} | Last error: {last_err}"
    )


def get_remote_file_size(url: str) -> int:
    """Return file size in bytes via HTTP HEAD. Returns 0 on failure."""
    try:
        r = requests.head(url, timeout=HTTP_TIMEOUT, allow_redirects=True)
        return int(r.headers.get("Content-Length", 0))
    except Exception:
        return 0


# ── Sample-size extraction ─────────────────────────────────────────────────


def extract_n_better(text: Any) -> Optional[int]:
    """
    Conservative N extraction from free-text sample descriptions.
    Prefers ``total N=xxxxx``; falls back to max number found.
    """
    if not isinstance(text, str):
        return None
    t = text.lower()

    m = re.search(r"(total\s*)?n\s*=\s*([\d,]+)", t)
    if m:
        return int(m.group(2).replace(",", ""))

    nums = re.findall(r"\d[\d,]*", text)
    vals = [int(n.replace(",", "")) for n in nums] if nums else []
    return max(vals) if vals else None


# ── Ancestry helpers ───────────────────────────────────────────────────────


def clean_discovery_ancestry(val: Any) -> str:
    if not isinstance(val, str):
        return ""
    return val.split(",")[0].strip()


def is_eur(text: Any) -> bool:
    if not isinstance(text, str):
        return False
    t = text.lower()
    return any(x in t for x in EURO_POS) and not any(x in t for x in NON_EURO_NEG)


# ── GTEx biosample + sample counts ────────────────────────────────────────


def load_gtex_table(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df.columns = [c.strip() for c in df.columns]
    return df


def _parse_tissue_from_eqtl_filename(filename: str) -> str:
    """
    Extract a human-readable tissue name from a GTEx eQTL filename.

    Examples:
      'Kidney_Cortex.v10.eGenes.txt'                          → 'kidney cortex'
      'Brain_Anterior_cingulate_cortex_BA24.v10.eGenes.txt'   → 'brain anterior cingulate cortex ba24'
      'Cells_EBV-transformed_lymphocytes.v10.eGenes.txt'      → 'cells ebv-transformed lymphocytes'
    """
    name = Path(filename).name
    # Strip everything from '.v<digits>' onward (version + extension)
    name = re.sub(r'\.v\d+.*$', '', name)
    return name.replace('_', ' ').lower().strip()


def find_manual_eqtl_file(
    biosample_input: str,
    eqtl_root: str,
) -> Optional[tuple]:
    """
    Check whether the user passed a GTEx eQTL filename (or path) directly.

    Accepts:
      - Bare filename:  'Kidney_Cortex.v10.eGenes.txt'
      - Full path:      '/path/to/GTEx_eQTL_TISSUE_EXPRESSION/Kidney_Cortex.v10.eGenes.txt'
      - Relative path:  'GTEx_eQTL_TISSUE_EXPRESSION/Kidney_Cortex.v10.eGenes.txt'

    Returns ``(Path, tissue_name_str)`` on success, or ``None`` if the input
    does not look like a file reference.
    """
    eqtl_root_path = Path(eqtl_root)

    # 1. Direct absolute / relative path that exists on disk
    p = Path(biosample_input)
    if p.exists() and p.is_file():
        return p, _parse_tissue_from_eqtl_filename(p.name)

    # 2. Bare filename inside eqtl_root
    candidate = eqtl_root_path / biosample_input
    if candidate.exists() and candidate.is_file():
        return candidate, _parse_tissue_from_eqtl_filename(candidate.name)

    # 3. Case-insensitive filename scan inside eqtl_root
    biosample_lower = biosample_input.lower()
    for p in eqtl_root_path.rglob("*.txt"):
        if p.name.lower() == biosample_lower:
            return p, _parse_tissue_from_eqtl_filename(p.name)

    return None


def _resolve_biosample_string(user_norm: str, gtex_biosamples_set: Set[str]) -> Optional[str]:
    """
    Pure string-based biosample resolver against the GTEx biosample CSV set.

    Stages:
      1. Exact match
      2. User term is a substring of a GTEx tissue name  (e.g. 'blood' → 'whole blood')
      3. Fuzzy character match with cutoff=0.6
    """
    if user_norm in gtex_biosamples_set:
        return user_norm

    for t in sorted(gtex_biosamples_set):
        if user_norm in t:
            return t

    matches = get_close_matches(user_norm, sorted(gtex_biosamples_set), n=1, cutoff=0.6)
    return matches[0] if matches else None


def _resolve_biosample_llm(
    user_biosample: str,
    gtex_biosamples_set: Set[str],
) -> Optional[str]:
    """
    Use GPT-4o-mini at temperature=0 to map any biological sample concept to the
    most appropriate GTEx tissue.

    The model receives the exact GTEx tissue list and no hardcoded examples —
    it uses its own biomedical knowledge to reason about the relationship between
    the sample type and available tissues.  It returns a ranked list of up to 3
    candidates; we accept the first one whose name appears verbatim in the GTEx set.
    """
    from decouple import config, UndefinedValueError
    try:
        api_key = config("OPENAI_API_KEY")
    except UndefinedValueError:
        raise RuntimeError("OPENAI_API_KEY not set — cannot use LLM for biosample resolution.")

    from openai import OpenAI
    client = OpenAI(api_key=api_key)

    tissue_list = "\n".join(f"- {t}" for t in sorted(gtex_biosamples_set))

    prompt = (
        "You are a biomedical expert. Your task is to map a biological sample type "
        "to the most appropriate GTEx tissue(s) from the provided list.\n\n"
        "STRICT RULES:\n"
        "1. Every tissue name in your response MUST be copied VERBATIM from the list "
        "below — do not alter spacing, capitalisation, or spelling.\n"
        "2. Apply your own biomedical knowledge to reason freshly for each input. "
        "Do NOT use any pre-programmed examples or rules.\n"
        "3. Return your top 3 choices ranked from most to least relevant.\n"
        "4. Temperature is 0 — be fully deterministic.\n\n"
        "BIOLOGICAL RELEVANCE RULE — this is critical:\n"
        "Set 'is_biologically_relevant': true ONLY when there is a DIRECT biological "
        "relationship between the sample type and the tissue:\n"
        "  - The tissue physically produces or secretes the sample "
        "(e.g. kidney produces urine, liver produces bile, brain produces CSF)\n"
        "  - The sample is directly derived from or contained within the tissue "
        "(e.g. PBMC from blood, synovial fluid from synovial membrane)\n"
        "  - The tissue is the primary biological origin of the cells or molecules "
        "in the sample\n"
        "Set 'is_biologically_relevant': false for:\n"
        "  - Anatomical proximity alone (organ is nearby but does not produce the sample)\n"
        "  - Indirect relationships: provides oxygen, receives waste, shares blood supply\n"
        "  - Functional co-regulation without direct sample derivation\n"
        "If NO tissue has a direct biological relationship, set all to false — "
        "do NOT pick a distant tissue just to give an answer.\n\n"
        f"Available GTEx tissues (copy names exactly):\n{tissue_list}\n\n"
        f"Biological sample type to map: \"{user_biosample}\"\n\n"
        "Return JSON with this exact structure:\n"
        '{\n'
        '  "ranked_tissues": [\n'
        '    {"tissue": "<exact name from list>", "reasoning": "<one sentence>", "is_biologically_relevant": true},\n'
        '    {"tissue": "<exact name from list>", "reasoning": "<one sentence>", "is_biologically_relevant": false},\n'
        '    {"tissue": "<exact name from list>", "reasoning": "<one sentence>", "is_biologically_relevant": false}\n'
        '  ]\n'
        '}'
    )

    resp = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.0,
        response_format={"type": "json_object"},
    )
    data = json.loads(resp.choices[0].message.content)
    ranked = data.get("ranked_tissues", [])

    for item in ranked:
        if not item.get("is_biologically_relevant", False):
            log(
                f"[LLM-BIOSAMPLE] Skipping '{item.get('tissue')}' — not directly relevant: "
                f"{item.get('reasoning', '')}"
            )
            continue
        tissue = str(item.get("tissue", "")).strip().lower()
        reasoning = item.get("reasoning", "")
        if tissue in gtex_biosamples_set:
            log(f"[LLM-BIOSAMPLE] '{user_biosample}' → '{tissue}' | {reasoning}")
            return tissue
        # Model may return mixed-case — normalise and retry
        if tissue and tissue in {t.lower() for t in gtex_biosamples_set}:
            matched = next(t for t in gtex_biosamples_set if t.lower() == tissue)
            log(f"[LLM-BIOSAMPLE] '{user_biosample}' → '{matched}' | {reasoning}")
            return matched

    log(
        f"[WARN] LLM found no biologically relevant tissue for '{user_biosample}'. "
        f"No indirect/proximity matches accepted. "
        f"LLM returned: {[i.get('tissue') for i in ranked]}"
    )
    return None


def _select_gwas_trait_llm(
    available_traits: List[str],
    disease_name: str,
    biosample_type: str,
) -> Optional[str]:
    """
    Use GPT-4o-mini at temperature=0 to pick the single most biologically
    appropriate GWAS trait from ``available_traits``, given both the disease
    name and the biosample/tissue context.

    Passing biosample context matters: e.g. disease='diabetes' + biosample='urine'
    should prefer 'diabetic nephropathy' over generic 'type 2 diabetes mellitus'
    when both are available.
    """
    from decouple import config, UndefinedValueError
    try:
        api_key = config("OPENAI_API_KEY")
    except UndefinedValueError:
        raise RuntimeError("OPENAI_API_KEY not set — cannot use LLM for GWAS trait selection.")

    from openai import OpenAI
    client = OpenAI(api_key=api_key)

    trait_list = "\n".join(f"- {t}" for t in sorted(set(available_traits)))

    prompt = (
        "You are a biomedical expert selecting the most appropriate GWAS study trait.\n\n"
        "CONTEXT:\n"
        f"  Target disease: {disease_name}\n"
        f"  Biosample / tissue type: {biosample_type}\n\n"
        "TASK:\n"
        "From the list of available GWAS traits below, choose the single trait that "
        "is most biologically relevant to the target disease AND consistent with the "
        "biosample/tissue type provided.\n\n"
        "STRICT RULES:\n"
        "1. Your answer MUST be copied VERBATIM from the list — do not alter spelling or case.\n"
        "2. Use your own biomedical knowledge to reason about the match. "
        "Consider disease subtypes, tissue specificity, and clinical context. "
        "Do NOT rely on any pre-programmed rules — reason freshly.\n"
        "3. If several traits are equally close, prefer the one most specific to the "
        "given biosample/tissue context.\n"
        "4. Temperature is 0 — give a single deterministic answer.\n\n"
        f"Available GWAS traits:\n{trait_list}\n\n"
        "Return JSON:\n"
        '{"selected_trait": "<exact trait name from list>", "reasoning": "<one sentence>"}'
    )

    resp = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.0,
        response_format={"type": "json_object"},
    )
    data = json.loads(resp.choices[0].message.content)
    selected = str(data.get("selected_trait", "")).strip()
    reasoning = data.get("reasoning", "")

    # Build a case-insensitive lookup
    trait_lower_map = {t.lower(): t for t in available_traits}
    matched = trait_lower_map.get(selected.lower())
    if matched:
        log(
            f"[LLM-TRAIT] disease='{disease_name}' + biosample='{biosample_type}' "
            f"→ '{matched}' | {reasoning}"
        )
        return matched

    log(
        f"[WARN] LLM selected trait '{selected}' not found verbatim in the "
        f"available traits list — discarding."
    )
    return None


def resolve_biosample(
    user_biosample: str,
    gtex_biosamples_set: Set[str],
    use_llm: bool = True,
) -> Optional[str]:
    """
    Resolve any biological sample description to a valid GTEx tissue name.

    Strategy (in order):
      1. String matching  — exact → substring → fuzzy ≥0.6
      2. LLM ranked mapping at temperature=0  (skipped if use_llm=False)
      3. Return None  (caller handles fallback / warning)
    """
    if not user_biosample:
        return None
    user_norm = user_biosample.lower().strip()

    result = _resolve_biosample_string(user_norm, gtex_biosamples_set)
    if result:
        return result

    if use_llm:
        try:
            result = _resolve_biosample_llm(user_biosample, gtex_biosamples_set)
            if result:
                return result
        except Exception as exc:
            log(f"[WARN] LLM biosample resolver failed: {exc}")

    return None


def gtex_sample_counts(
    gtex_df: pd.DataFrame, biosample_lower: str
) -> tuple:
    """
    Pull RNA-seq and genotype sample counts from the GTEx biosample CSV.
    Returns ``(rnaseq_n, genotype_n)`` or ``(None, None)``.
    """
    if "Tissue" not in gtex_df.columns:
        return (None, None)

    sub = gtex_df[gtex_df["Tissue"].astype(str).str.lower().str.strip() == biosample_lower]
    if sub.empty:
        sub = gtex_df[
            gtex_df["Tissue"]
            .astype(str)
            .str.lower()
            .str.contains(re.escape(biosample_lower), na=False)
        ]
        if sub.empty:
            return (None, None)
        sub = sub.head(1)

    row = sub.iloc[0].to_dict()

    rnaseq_col: Optional[str] = None
    geno_col: Optional[str] = None
    for c in gtex_df.columns:
        cl = c.lower()
        if rnaseq_col is None and (
            ("rna" in cl and ("seq" in cl or "sequencing" in cl)) or ("rnaseq" in cl)
        ):
            rnaseq_col = c
        if geno_col is None and ("geno" in cl or "genotype" in cl):
            geno_col = c

    def _to_int(v: Any) -> Optional[int]:
        try:
            if pd.isna(v):
                return None
            if isinstance(v, str):
                v = v.replace(",", "").strip()
            return int(float(v))
        except Exception:
            return None

    rnaseq_n = _to_int(row.get(rnaseq_col)) if rnaseq_col else None
    genotype_n = _to_int(row.get(geno_col)) if geno_col else None
    return rnaseq_n, genotype_n


# ── eQTL file finding ─────────────────────────────────────────────────────


def find_eqtl_file(eqtl_root: str, biosample: str) -> Optional[Path]:
    key = safe_name(biosample)
    for p in Path(eqtl_root).rglob("*"):
        if p.is_file() and key in safe_name(p.name):
            return p
    return None


# ── LLM disease-name normalisation ────────────────────────────────────────


def _normalize_disease_llm(user_query: str) -> tuple:
    from decouple import config, UndefinedValueError

    try:
        api_key = config("OPENAI_API_KEY")
    except UndefinedValueError:
        raise RuntimeError(
            "OPENAI_API_KEY not found. Set it in a .env file in your "
            "working directory, as a system environment variable, or in "
            "a settings.ini file."
        )

    from openai import OpenAI  # deferred import – package is optional

    client = OpenAI(api_key=api_key)
    prompt = (
        "You are a biomedical disease-name normalizer.\n\n"
        "Rules:\n"
        "- Correct spelling errors\n"
        "- If already valid, return unchanged\n"
        "- Output standard clinical disease name\n"
        "- If unknown, return UNKNOWN\n\n"
        "Return JSON:\n"
        '{"normalized_name": "...", "synonyms": ["...", "..."]}\n\n'
        f"Input:\n{user_query}"
    )
    resp = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.0,
        response_format={"type": "json_object"},
    )
    data = json.loads(resp.choices[0].message.content)
    return data.get("normalized_name", user_query), data.get("synonyms", [])


def normalize_disease(user_query: str, use_llm: bool = True) -> tuple:
    if not use_llm:
        return user_query, []
    try:
        return _normalize_disease_llm(user_query)
    except Exception as e:
        log(f"[WARN] LLM normalization failed, fallback to raw input. Reason: {e}")
        return user_query, []


# ── Option A+ exact phenotype matching ─────────────────────────────────────


def filter_exact_with_fallback(
    hits_df: pd.DataFrame,
    primary_name: str,
    synonyms: List[str],
    wait_seconds: int = 5,
) -> pd.DataFrame:
    """
    Exact-match filter with synonym fallback.
    Prevents accidental inclusion of unrelated phenotype subtypes.
    """
    if hits_df.empty:
        return hits_df

    df = hits_df.copy()
    df["efo_clean"] = df["efoTraits"].astype(str).str.lower().str.strip()

    target = primary_name.lower().strip()
    exact = df[df["efo_clean"] == target]
    if not exact.empty:
        log(f"[OK] Exact GWAS trait match found for '{primary_name}'")
        return exact

    if synonyms:
        log(
            f"[INFO] No exact GWAS for '{primary_name}'. "
            f"Trying synonyms after {wait_seconds}s..."
        )
        time.sleep(wait_seconds)
        for syn in synonyms:
            syn_clean = str(syn).lower().strip()
            if not syn_clean:
                continue
            syn_hits = df[df["efo_clean"] == syn_clean]
            if not syn_hits.empty:
                log(f"[OK] Exact GWAS trait match found using synonym '{syn}'")
                return syn_hits

    log("[STOP] No exact GWAS trait match found for disease or its synonyms.")
    try:
        sample_traits = sorted(set(df["efoTraits"].astype(str)))[:10]
        log("Nearby available traits (sample):")
        for t in sample_traits:
            log(f"  - {t}")
    except Exception:
        pass
    return df.iloc[0:0]


# ── Download + decompress ──────────────────────────────────────────────────


def download_resume(url: str, out: str, max_retries: int = 5) -> None:
    """Resume-safe downloader for large files with progress indicator."""
    out = Path(out)
    retries = 0

    while retries < max_retries:
        try:
            headers: Dict[str, str] = {}
            if out.exists():
                existing_size = out.stat().st_size
                headers["Range"] = f"bytes={existing_size}-"
                log(f"Resuming download at byte {existing_size:,}")
            else:
                existing_size = 0
                log("Starting new download")

            with requests.get(url, stream=True, headers=headers, timeout=HTTP_TIMEOUT) as r:
                if r.status_code == 416:
                    log("Server rejected range request. Restarting download.")
                    out.unlink(missing_ok=True)
                    continue

                r.raise_for_status()
                # Guard: if we sent a Range header but the server returned 200
                # (full content) instead of 206, appending would corrupt the file.
                # Delete the partial file and restart from byte 0.
                if existing_size > 0 and r.status_code == 200:
                    log(
                        "[WARN] Server ignored Range request (returned 200, not 206). "
                        "Deleting partial file and restarting from byte 0 to avoid corruption."
                    )
                    out.unlink(missing_ok=True)
                    existing_size = 0
                remaining = int(r.headers.get("content-length", 0))
                full_size = existing_size + remaining
                downloaded = 0
                mode = "ab" if existing_size > 0 else "wb"
                with open(out, mode) as f:
                    for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                        if chunk:
                            f.write(chunk)
                            downloaded += len(chunk)
                            if full_size > 0:
                                pct = (existing_size + downloaded) / full_size * 100
                                done_mb = (existing_size + downloaded) / 1e6
                                total_mb = full_size / 1e6
                                print(
                                    f"\r  Downloading... {pct:.1f}%  "
                                    f"({done_mb:.1f} / {total_mb:.1f} MB)",
                                    end="", flush=True,
                                )
                if full_size > 0:
                    print()

            log(f"Download completed: {out}")
            return

        except (ChunkedEncodingError, ConnectionError, Timeout) as e:
            retries += 1
            log(f"[WARN] Download interrupted ({retries}/{max_retries}). Retrying... {e}")
            time.sleep(RETRY_SLEEP)

        except Exception as e:
            retries += 1
            log(f"[WARN] Unexpected download error ({retries}/{max_retries}): {e}")
            time.sleep(RETRY_SLEEP)

    raise RuntimeError(f"Failed to download after {max_retries} retries: {url}")


def gunzip_file(gz_path: Path, out_path: Path) -> bool:
    """Unzip .gz to *out_path* (overwrite-safe).

    Removes the partial output file on failure so downstream code
    never picks up a corrupt/incomplete TSV.
    """
    try:
        with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        return True
    except Exception as e:
        log(f"[WARN] Failed to unzip {gz_path} -> {out_path}: {e}")
        Path(out_path).unlink(missing_ok=True)
        return False


# ── Harmonised folder mirroring ────────────────────────────────────────────


def list_harmonised_files(harmonised_url: str) -> List[str]:
    """Parse directory listing HTML for href links."""
    try:
        resp = safe_request_get(harmonised_url, stream=False)
        if resp is None:
            return []
        files = re.findall(r'href="([^"]+)"', resp.text)
        files = [f for f in files if f and not f.startswith("?") and not f.startswith("/")]
        return sorted(set(files))
    except Exception as e:
        log(f"[WARN] Could not list harmonised files: {harmonised_url} | {e}")
        return []


def mirror_harmonised_folder(
    harmonised_url: str,
    out_dir: Path,
    preloaded_gz: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    Download all files from harmonised folder into ``out_dir/harmonised/``.
    Unzips the first ``.h.tsv.gz`` encountered.

    If ``preloaded_gz`` is provided and points to an existing file, the
    ``.h.tsv.gz`` is moved from that path instead of being re-downloaded,
    saving several GB of bandwidth.
    """
    harm_dir = out_dir / "harmonised"
    harm_dir.mkdir(parents=True, exist_ok=True)

    files = list_harmonised_files(harmonised_url)
    if not files:
        return {"harmonised_dir": str(harm_dir), "files": [], "gwas_gz": None, "gwas_tsv": None}

    downloaded: List[str] = []
    gwas_gz: Optional[Path] = None
    gwas_tsv: Optional[Path] = None

    for f in files:
        if f.endswith("/"):
            continue
        out_path = harm_dir / f

        if f.endswith(".h.tsv.gz") and gwas_gz is None:
            # Reuse the file already downloaded during overlap evaluation
            if preloaded_gz is not None and preloaded_gz.exists():
                if not out_path.exists():
                    log(f"[REUSE] Moving pre-evaluated file to final destination (skipping re-download): {f}")
                    shutil.move(str(preloaded_gz), str(out_path))
                else:
                    log(f"[SKIP] {f} already present at destination")
                downloaded.append(str(out_path))
                gwas_gz = out_path
                tsv_out = harm_dir / f.replace(".gz", "")
                if gunzip_file(gwas_gz, tsv_out):
                    gwas_tsv = tsv_out
                    log(f"Unzipped -> {tsv_out}")
                continue

        url = harmonised_url.rstrip("/") + "/" + f
        try:
            log(f"Downloading harmonised file: {f}")
            download_resume(url, str(out_path), max_retries=HTTP_RETRIES)
            downloaded.append(str(out_path))

            if f.endswith(".h.tsv.gz") and gwas_gz is None:
                gwas_gz = out_path
                tsv_out = harm_dir / f.replace(".gz", "")
                if gunzip_file(gwas_gz, tsv_out):
                    gwas_tsv = tsv_out
                    log(f"Unzipped -> {tsv_out}")
        except Exception as e:
            log(f"[WARN] Failed download for {url}: {e}")
            continue

    return {
        "harmonised_dir": str(harm_dir),
        "files": downloaded,
        "gwas_gz": str(gwas_gz) if gwas_gz else None,
        "gwas_tsv": str(gwas_tsv) if gwas_tsv else None,
    }


# ── Reporting ──────────────────────────────────────────────────────────────


def write_dataset_excel(out_dir: Path, dataset_name: str, row_dict: dict) -> Path:
    """Write an Excel report inside each dataset folder."""
    xlsx_path = out_dir / f"{dataset_name}.xlsx"
    df = pd.DataFrame([row_dict])
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="dataset")
    return xlsx_path


# ── Cache detection ────────────────────────────────────────────────────────


def find_existing_local_dataset(out_root: Path, normalized_trait: str) -> Optional[Path]:
    """Detect cached GWAS data (layout-agnostic)."""
    trait_key = safe_name(normalized_trait)
    out_root = Path(out_root)

    if not out_root.exists():
        return None

    for d in out_root.iterdir():
        if not d.is_dir():
            continue
        if safe_name(d.name) != trait_key:
            continue

        harm = d / "harmonised"
        if harm.exists() and any(harm.iterdir()):
            return d

        files = list(d.glob("*"))
        if any(
            f.is_file()
            and (
                f.name.lower().endswith(".meta.yaml")
                or f.name.lower().endswith(".yaml")
                or f.name.lower().endswith(".tsv")
                or f.name.lower().endswith(".tsv.gz")
                or f.name.lower().endswith(".gz")
            )
            for f in files
        ):
            return d

    return None


def detect_cached_gwas_files(existing_dataset_dir: Path) -> Dict[str, Any]:
    """
    Search for GWAS + meta files in both the dataset directory
    and its ``harmonised/`` subfolder.
    """
    existing_dataset_dir = Path(existing_dataset_dir)
    search_dirs = [existing_dataset_dir, existing_dataset_dir / "harmonised"]

    gwas_tsv: Optional[Path] = None
    gwas_gz: Optional[Path] = None
    meta: Optional[Path] = None

    for sd in search_dirs:
        if not sd.exists():
            continue
        if gwas_tsv is None:
            gwas_tsv = next(sd.glob("*.h.tsv"), None)
        if gwas_gz is None:
            gwas_gz = next(sd.glob("*.h.tsv.gz"), None)
        if meta is None:
            meta = next(sd.glob("*.meta.yaml"), None)

    if (existing_dataset_dir / "harmonised").exists() and any(
        (existing_dataset_dir / "harmonised").iterdir()
    ):
        harmonised_dir_used = existing_dataset_dir / "harmonised"
        layout = "newstyle_harmonised_subfolder"
    else:
        harmonised_dir_used = existing_dataset_dir
        layout = "legacy_or_root_files"

    if gwas_tsv:
        status = f"exists_unzipped_{layout}"
    elif gwas_gz:
        status = f"exists_gz_{layout}"
    elif meta:
        status = f"exists_meta_only_{layout}"
    else:
        status = f"exists_partial_only_{layout}"

    return {
        "gwas_tsv": str(gwas_tsv) if gwas_tsv else None,
        "gwas_gz": str(gwas_gz) if gwas_gz else None,
        "meta": str(meta) if meta else None,
        "harmonised_dir_used": str(harmonised_dir_used),
        "status": status,
    }


# ── Tiered dataset selection helpers ────────────────────────────────────────


def load_eqtl_instrument_snps(
    eqtl_path: str,
    p_thresh: float = EQTL_P_THRESH_INSTRUMENT,
    fstat_thresh: float = EQTL_FSTAT_THRESH_INSTRUMENT,
) -> Set[str]:
    """Return rsIDs from an eQTL file that pass instrument-selection filters.

    Mirrors the R pipeline thresholds:
    ``pval_nominal < EQTL_P_THRESH`` AND ``F_stat > FSTAT_THRESH``.
    """
    df = pd.read_csv(eqtl_path, sep="\t", dtype=str)
    required = ("pval_nominal", "slope", "slope_se", "rs_id_dbSNP155_GRCh38p13")
    for col in required:
        if col not in df.columns:
            log(f"[WARN] eQTL file missing column '{col}'; "
                "cannot extract instruments")
            return set()

    df["pval_nominal"] = pd.to_numeric(df["pval_nominal"], errors="coerce")
    df["slope"] = pd.to_numeric(df["slope"], errors="coerce")
    df["slope_se"] = pd.to_numeric(df["slope_se"], errors="coerce")

    df = df[df["pval_nominal"].notna() & (df["pval_nominal"] < p_thresh)]
    df["F_stat"] = (df["slope"] / df["slope_se"]) ** 2
    df = df[df["F_stat"] > fstat_thresh]

    rsids = df["rs_id_dbSNP155_GRCh38p13"].dropna().astype(str)
    return set(rsids[rsids.str.match(r"^rs\d+$")])


def extract_gwas_rsids(path: Path) -> Set[str]:
    """Extract the set of valid rsIDs from a harmonised GWAS file.

    Handles ``.tsv`` and ``.tsv.gz`` (pandas infers compression).
    """
    try:
        df = pd.read_csv(
            path, sep="\t", usecols=["rsid"],
            dtype={"rsid": str}, on_bad_lines="skip",
        )
        rsids = df["rsid"].dropna()
        return set(rsids[rsids.str.startswith("rs")])
    except ValueError:
        hdr = pd.read_csv(path, sep="\t", nrows=0)
        rs_col = next(
            (c for c in hdr.columns if c.lower() in ("rsid", "rs_id", "snp")),
            None,
        )
        if rs_col is None:
            return set()
        df = pd.read_csv(
            path, sep="\t", usecols=[rs_col],
            dtype={rs_col: str}, on_bad_lines="skip",
        )
        return set(df[rs_col].dropna()[df[rs_col].str.startswith("rs")])
    except Exception:
        return set()


def evaluate_candidate_overlap(
    harmonised_url: str,
    gz_filename: str,
    eqtl_snps: Set[str],
) -> Dict[str, Any]:
    """Download a candidate ``.h.tsv.gz`` and compute eQTL overlap.

    On success the downloaded file is kept so the caller can reuse it for the
    final download step (avoiding a second download of the same file).
    The caller is responsible for cleaning up ``tmp_dir`` when done.
    On failure the temp dir is cleaned up internally and tmp_path/tmp_dir are None.
    """
    gz_url = harmonised_url.rstrip("/") + "/" + gz_filename
    tmp_dir = Path(tempfile.mkdtemp(prefix="gwas_eval_"))
    tmp_path = tmp_dir / gz_filename

    try:
        download_resume(gz_url, str(tmp_path))
        gwas_rsids = extract_gwas_rsids(tmp_path)
        overlap = gwas_rsids & eqtl_snps
        return {
            "variant_count": len(gwas_rsids),
            "overlap_count": len(overlap),
            "tmp_path": tmp_path,
            "tmp_dir": tmp_dir,
        }
    except Exception as exc:
        log(f"[WARN] Overlap evaluation failed for {gz_filename}: {exc}")
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return {"variant_count": 0, "overlap_count": 0, "tmp_path": None, "tmp_dir": None}


# ═══════════════════════════════════════════════════════════════════════════
#  Main retrieval function
# ═══════════════════════════════════════════════════════════════════════════


def retrieve_gwas_data(
    disease_name: str,
    biosample_type: str,
    out_root: str,
    *,
    tsv_path: Optional[str] = None,
    eqtl_root: Optional[str] = None,
    gtex_biosample_csv: Optional[str] = None,
    use_llm: bool = True,
    top_k: int = 1,
) -> List[Dict[str, Any]]:
    """
    Retrieve GWAS + eQTL data for a disease / biosample type combination.

    This is the programmatic equivalent of running
    ``gwas_eqtl_data_retrieval_pipeline.py`` interactively.

    Parameters
    ----------
    disease_name : str
        Disease or phenotype name (e.g. ``"Breast Cancer"``).
    biosample_type : str
        GTEx biosample type (e.g. ``"Whole Blood"``, ``"Pancreas"``, ``"PBMC"``).
    out_root : str
        Output directory for downloaded GWAS data.
    tsv_path : str, optional
        Path to the GWAS catalog TSV.  Defaults to
        ``gwas_mr_reference/GWAS-DATABASE-FTP.tsv``.
    eqtl_root : str, optional
        Path to the GTEx eQTL biosample expression directory.  Defaults to
        ``gwas_mr_reference/GTEx_eQTL_TISSUE_EXPRESSION/``.
    gtex_biosample_csv : str, optional
        Path to the GTEx biosample sample-count CSV.  Defaults to
        ``gwas_mr_reference/GTEx_EQTL_SAMPLE_COUNTS/EQTL-SAMPLE-COUNT-GTEx Portal.csv``.
    use_llm : bool
        Use OpenAI LLM for disease-name normalisation (default ``True``).
        The API key is read via ``python-decouple`` (checks ``.env``,
        environment variables, and ``settings.ini`` automatically).
    top_k : int
        Number of top GWAS rows per trait by sample size (default ``1``).

    Returns
    -------
    list[dict]
        One dict per processed trait containing paths, sample sizes, and
        metadata needed to feed the downstream MR pipeline.
    """
    from .defaults import (
        DEFAULT_TSV_PATH,
        DEFAULT_EQTL_ROOT,
        DEFAULT_GTEX_BIOSAMPLE_CSV,
    )

    tsv_path = tsv_path or DEFAULT_TSV_PATH
    eqtl_root = eqtl_root or DEFAULT_EQTL_ROOT
    gtex_biosample_csv = gtex_biosample_csv or DEFAULT_GTEX_BIOSAMPLE_CSV
    # ── Disease normalisation ──────────────────────────────────────────────
    normalized_name, synonyms = normalize_disease(disease_name, use_llm=use_llm)
    log(f"Normalized disease: {normalized_name}")
    if synonyms:
        log(f"Synonyms: {', '.join(synonyms[:10])}")

    # ── Load GWAS table ────────────────────────────────────────────────────
    log("\nLoading GWAS FTP table (this can take time)...")
    df = pd.read_csv(tsv_path, sep="\t", dtype=str, on_bad_lines="skip")

    df["discoverySampleAncestry"] = df["discoverySampleAncestry"].apply(
        clean_discovery_ancestry
    )
    df["N"] = df["initialSampleDescription"].apply(extract_n_better)
    df["is_EUR"] = df["initialSampleDescription"].apply(is_eur) | df[
        "discoverySampleAncestry"
    ].apply(is_eur)

    def _search_terms(terms: List[str]) -> pd.DataFrame:
        mask = (
            df["efoTraits"]
            .fillna("")
            .str.lower()
            .apply(
                lambda x: any(
                    t.lower() in x for t in terms if isinstance(t, str) and t.strip()
                )
            )
        )
        return df[mask & df["summaryStatistics"].notna() & df["is_EUR"]]

    hits = _search_terms([normalized_name])
    if hits.empty and synonyms:
        hits = _search_terms(synonyms)
    broad_hits = hits.copy()  # keep for LLM fallback before exact filtering narrows it

    if hits.empty:
        log("No GWAS dataset found via broad match — will try LLM trait selection from all EUR studies.")

    log(f"Found {len(hits)} GWAS rows (broad match).")

    # ── Option A+ narrowing ────────────────────────────────────────────────
    hits = filter_exact_with_fallback(
        hits_df=hits,
        primary_name=normalized_name,
        synonyms=synonyms,
        wait_seconds=5,
    )

    # ── LLM trait fallback (two levels) ───────────────────────────────────
    if hits.empty and use_llm:
        # Level 1: pick from the broad substring-match results
        if not broad_hits.empty:
            broad_traits = sorted(set(broad_hits["efoTraits"].dropna().astype(str).tolist()))
            log(
                f"[LLM-TRAIT] Exact match failed. Asking LLM to pick closest trait "
                f"from {len(broad_traits)} broad-match candidates "
                f"(disease='{normalized_name}', biosample='{biosample_type}')..."
            )
            try:
                selected = _select_gwas_trait_llm(broad_traits, normalized_name, biosample_type)
                if selected:
                    hits = broad_hits[
                        broad_hits["efoTraits"].astype(str).str.lower().str.strip()
                        == selected.lower().strip()
                    ]
            except Exception as exc:
                log(f"[WARN] LLM trait selection (level 1) failed: {exc}")

        # Level 2: broad search was also empty — try from ALL EUR studies
        if hits.empty:
            log(
                f"[LLM-TRAIT] Broad search empty or LLM level-1 failed. "
                f"Expanding to all EUR GWAS studies..."
            )
            try:
                all_eur = df[df["summaryStatistics"].notna() & df["is_EUR"]]
                trait_series = all_eur["efoTraits"].dropna().astype(str)
                # Sort by study count descending so the 300-cap keeps the most-studied
                # traits (covers all common diseases regardless of starting letter).
                trait_counts = trait_series.value_counts()
                all_traits = trait_counts.index.tolist()[:300]
                log(
                    f"[LLM-TRAIT] Asking LLM to pick from {len(all_traits)} EUR traits "
                    f"(disease='{normalized_name}', biosample='{biosample_type}')..."
                )
                selected = _select_gwas_trait_llm(all_traits, normalized_name, biosample_type)
                if selected:
                    hits = all_eur[
                        all_eur["efoTraits"].astype(str).str.lower().str.strip()
                        == selected.lower().strip()
                    ]
            except Exception as exc:
                log(f"[WARN] LLM trait selection (level 2) failed: {exc}")

    if hits.empty:
        log("\nPipeline stopped: could not find a matching GWAS trait (exact, synonym, or LLM).")
        return []

    log(f"Proceeding with phenotype: '{hits.iloc[0]['efoTraits']}'")

    # ── Biosample resolution ─────────────────────────────────────────────
    gtex_df = load_gtex_table(gtex_biosample_csv)
    gtex_biosamples = set(gtex_df["Tissue"].astype(str).str.lower().str.strip())

    # Path A — user passed an eQTL filename or file path directly
    #   e.g.  biosample_type = "Kidney_Cortex.v10.eGenes.txt"
    #   or    biosample_type = "GTEx_eQTL_TISSUE_EXPRESSION/Kidney_Cortex.v10.eGenes.txt"
    manual_result = find_manual_eqtl_file(biosample_type, eqtl_root)
    if manual_result:
        eqtl_file, parsed_tissue = manual_result
        log(f"[MANUAL] eQTL file provided directly: {eqtl_file.name}")
        log(f"[MANUAL] Parsed tissue name: '{parsed_tissue}'")
        # Resolve parsed tissue name to the GTEx biosample CSV for sample counts
        user_biosample = resolve_biosample(parsed_tissue, gtex_biosamples, use_llm=use_llm) or parsed_tissue
        log(f"[MANUAL] Mapped to GTEx biosample entry: '{user_biosample}'")
    else:
        # Path B — user passed a tissue name / biosample description; resolve it
        biosample_lower = biosample_type.strip().lower()
        resolved_biosample = resolve_biosample(biosample_lower, gtex_biosamples, use_llm=use_llm)
        if not resolved_biosample:
            sorted_tissues = sorted(gtex_biosamples)
            log(
                f"\n[ERROR] Could not map biosample '{biosample_type}' to any GTEx tissue.\n"
                f"  String matching and LLM both returned no biologically relevant match.\n"
                f"  To avoid using an irrelevant tissue, the pipeline cannot proceed.\n"
                f"  Use MANUAL mode with an exact GTEx filename, e.g.:\n"
                f"    -Biosample \"Kidney_Cortex.v10.eGenes.txt\"\n"
                f"  Available GTEx tissues:\n"
                + "\n".join(f"    - {t}" for t in sorted_tissues)
            )
            return []
        elif resolved_biosample != biosample_lower:
            log(f"[OK] Biosample '{biosample_type}' resolved to GTEx tissue: '{resolved_biosample}'")

        user_biosample = resolved_biosample
        eqtl_file = find_eqtl_file(eqtl_root, user_biosample)

    rnaseq_n, genotype_n = gtex_sample_counts(gtex_df, user_biosample)
    log(f"Selected biosample type: {user_biosample}")
    log(f"Selected eQTL file: {eqtl_file}")
    log(f"GTEx RNA-seq N: {rnaseq_n} | Genotype N: {genotype_n}")

    # ── Load eQTL instrument SNPs for tiered selection ─────────────────────
    eqtl_instrument_snps: Set[str] = set()
    if eqtl_file:
        try:
            eqtl_instrument_snps = load_eqtl_instrument_snps(str(eqtl_file))
            log(
                f"[TIERED] Loaded {len(eqtl_instrument_snps):,} eQTL "
                f"instrument SNPs (p<{EQTL_P_THRESH_INSTRUMENT}, "
                f"F>{EQTL_FSTAT_THRESH_INSTRUMENT})"
            )
        except Exception as exc:
            log(f"[WARN] Could not load eQTL instruments: {exc}")
    else:
        log("[INFO] No eQTL file; tiered overlap selection will be skipped.")

    # ── Create output root ─────────────────────────────────────────────────
    out_root_path = Path(out_root)
    out_root_path.mkdir(parents=True, exist_ok=True)
    master_rows: List[Dict[str, Any]] = []

    log("\nStarting GWAS selection & download...")

    # ── Trait loop (offline-first + single pass) ───────────────────────────
    for trait, sub in hits.groupby("efoTraits"):
        log("\n==================================================")
        log(f"Processing trait: {trait}")

        harmonised_url: Optional[str] = None
        existing_dataset_dir = find_existing_local_dataset(out_root_path, trait)

        if existing_dataset_dir:
            out_dir = Path(existing_dataset_dir)
            data_status = "cached"
            log(f"[CACHE] Found existing local dataset folder: {out_dir}")

            cached_info = detect_cached_gwas_files(out_dir)
            log(f"[CACHE] Layout status: {cached_info['status']}")

            gwas_tsv = cached_info.get("gwas_tsv")
            gwas_gz = cached_info.get("gwas_gz")
            harm_dir_used = Path(cached_info["harmonised_dir_used"])

            if gwas_tsv is None and gwas_gz:
                gwas_tsv_path = Path(gwas_gz).with_suffix("")
                if not gwas_tsv_path.exists():
                    if gunzip_file(Path(gwas_gz), gwas_tsv_path):
                        log(f"[CACHE] Unzipped cached GWAS -> {gwas_tsv_path}")
                    else:
                        log(
                            f"[WARN] Corrupt .gz detected. Deleting "
                            f"{gwas_gz} so next run will re-download."
                        )
                        Path(gwas_gz).unlink(missing_ok=True)
                        gwas_gz = None
                gwas_tsv = str(gwas_tsv_path) if gwas_tsv_path.exists() else None

            if not gwas_tsv and not gwas_gz:
                log(
                    "[WARN] Cached dataset found but no GWAS TSV/GZ present "
                    "(metadata-only / partial). Delete the cache folder and "
                    "re-run to trigger a fresh download."
                )

            mirror_info: Dict[str, Any] = {
                "harmonised_dir": str(harm_dir_used),
                "files": [],
                "gwas_gz": gwas_gz,
                "gwas_tsv": gwas_tsv,
            }

            best = sub.iloc[0]
            try:
                harmonised_url = (
                    str(best["summaryStatistics"]).rstrip("/") + "/harmonised/"
                )
            except Exception:
                pass

        else:
            out_dir = out_root_path / safe_name(trait)
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "harmonised").mkdir(parents=True, exist_ok=True)

            data_status = "downloaded"
            sub = sub.copy()
            sub["N_num"] = pd.to_numeric(sub["N"], errors="coerce")
            sub = sub.sort_values("N_num", ascending=False)

            # ── Tiered Phase 1: minimum sample-size gate ───────────
            qualified = sub[sub["N_num"] >= MIN_GWAS_N]
            if qualified.empty:
                log(
                    f"[TIERED] No studies meet min-N threshold "
                    f"({MIN_GWAS_N:,}); using all {len(sub)} candidates."
                )
                eval_pool = sub.head(TIERED_EVAL_TOP_K)
            else:
                log(
                    f"[TIERED] {len(qualified)} studies meet min-N "
                    f"threshold ({MIN_GWAS_N:,})"
                )
                eval_pool = qualified.head(TIERED_EVAL_TOP_K)

            # ── Tiered Phase 2: harmonised-folder quality gate ─────
            candidates: List[Dict[str, Any]] = []
            for _, r in eval_pool.iterrows():
                harm_url = (
                    str(r["summaryStatistics"]).rstrip("/") + "/harmonised/"
                )
                log(f"Checking harmonised folder: {harm_url}")
                try:
                    resp = safe_request_get(harm_url, stream=False)
                    if resp is None:
                        log(f"[SKIP] 404 for harmonised URL")
                        continue
                    html = resp.text
                except Exception as e:
                    log(f"[WARN] Cannot open harmonised URL: {e}")
                    continue

                harm_score = sum([
                    ".h.tsv.gz" in html,
                    "meta.yaml" in html,
                    "README" in html,
                ])
                if harm_score < MIN_HARMONISED_FILE_TYPES:
                    log(
                        f"[SKIP] Harmonised score {harm_score} < "
                        f"{MIN_HARMONISED_FILE_TYPES}"
                    )
                    continue

                gz_files = re.findall(
                    r'href="([^"]*\.h\.tsv\.gz)"', html
                )
                gz_fname = gz_files[0] if gz_files else None

                # Fast file-size probe via HEAD request
                file_size = 0
                if gz_fname:
                    gz_url = harm_url.rstrip("/") + "/" + gz_fname
                    file_size = get_remote_file_size(gz_url)
                    log(
                        f"  {r.get('accessionId', '?')}: "
                        f".h.tsv.gz = {file_size / 1e6:.1f} MB"
                    )

                candidates.append({
                    "row": r,
                    "harm_url": harm_url,
                    "harm_score": harm_score,
                    "gz_filename": gz_fname,
                    "file_size": file_size,
                    "n_num": (
                        float(r["N_num"]) if pd.notna(r["N_num"]) else 0
                    ),
                    "overlap_count": -1,
                    "variant_count": 0,
                    "tmp_path": None,
                    "tmp_dir": None,
                })

            if not candidates:
                log(
                    "[WARN] No harmonised GWAS found online "
                    "(score threshold not met)."
                )
                continue

            # Sort by file size descending so the most promising
            # candidates are evaluated first for overlap.
            candidates.sort(key=lambda c: c["file_size"], reverse=True)

            # ── Tiered Phase 3: eQTL overlap evaluation ────────────
            #    Download .h.tsv.gz for the top candidates (by file
            #    size) and compute actual rsID overlap with eQTL
            #    instruments.  Limited to OVERLAP_EVAL_LIMIT to avoid
            #    excessive bandwidth.
            if eqtl_instrument_snps and len(candidates) > 1:
                n_inst = len(eqtl_instrument_snps)
                eval_subset = candidates[:OVERLAP_EVAL_LIMIT]
                log(
                    f"\n[TIERED] Evaluating top {len(eval_subset)} of "
                    f"{len(candidates)} candidates for eQTL instrument "
                    f"overlap ({n_inst:,} instrument SNPs)..."
                )

                for cand in eval_subset:
                    if not cand["gz_filename"]:
                        cand["overlap_count"] = 0
                        continue
                    accession = cand["row"].get("accessionId", "?")
                    log(
                        f"  Evaluating {accession} "
                        f"(N={cand['n_num']:,.0f}, "
                        f"size={cand['file_size'] / 1e6:.1f} MB)..."
                    )
                    info = evaluate_candidate_overlap(
                        cand["harm_url"],
                        cand["gz_filename"],
                        eqtl_instrument_snps,
                    )
                    cand["variant_count"] = info["variant_count"]
                    cand["overlap_count"] = info["overlap_count"]
                    cand["tmp_path"] = info.get("tmp_path")
                    cand["tmp_dir"] = info.get("tmp_dir")
                    pct = cand["overlap_count"] / n_inst * 100
                    log(
                        f"    variants={cand['variant_count']:,}, "
                        f"eQTL overlap={cand['overlap_count']:,} "
                        f"({pct:.1f}%)"
                    )

                # ── Tiered Phase 4: rank & select ─────────────────
                # Primary key: overlap count (more instruments = more
                # genes testable with robust MR methods).
                # Tie-break: sample size (better per-SNP precision).
                candidates.sort(
                    key=lambda c: (c["overlap_count"], c["n_num"]),
                    reverse=True,
                )
                selected = candidates[0]

                if selected["overlap_count"] > 0:
                    runner = (
                        candidates[1] if len(candidates) > 1 else None
                    )
                    msg = (
                        f"[TIERED] Selected "
                        f"{selected['row'].get('accessionId', '?')}: "
                        f"overlap={selected['overlap_count']:,}, "
                        f"N={selected['n_num']:,.0f}"
                    )
                    if runner and runner["overlap_count"] >= 0:
                        msg += (
                            f" (over "
                            f"{runner['row'].get('accessionId', '?')}: "
                            f"overlap={runner.get('overlap_count', 0):,}"
                            f", N={runner.get('n_num', 0):,.0f})"
                        )
                    log(msg)
                else:
                    log(
                        "[TIERED] No candidate has eQTL overlap. "
                        "Falling back to largest sample size."
                    )
                    candidates.sort(
                        key=lambda c: c["n_num"], reverse=True
                    )
                    selected = candidates[0]

                best = selected["row"]

            else:
                # Single candidate or no eQTL instruments — fall back
                # to best harmonised score / largest N.
                if len(candidates) == 1:
                    log(
                        f"[INFO] Single candidate: "
                        f"{candidates[0]['row'].get('accessionId', '?')}"
                    )
                else:
                    log(
                        "[INFO] eQTL instruments unavailable; "
                        "selecting by sample size."
                    )
                candidates.sort(
                    key=lambda c: (c["harm_score"], c["n_num"]),
                    reverse=True,
                )
                selected = candidates[0]
                best = selected["row"]

            # Reuse the file already downloaded during overlap evaluation for
            # the winning candidate, and clean up all loser temp dirs now.
            winner_tmp_gz: Optional[Path] = selected.get("tmp_path")
            winner_tmp_dir: Optional[Path] = selected.get("tmp_dir")
            for cand in candidates:
                td = cand.get("tmp_dir")
                if td is not None and td != winner_tmp_dir:
                    shutil.rmtree(str(td), ignore_errors=True)

            harmonised_url = (
                str(best["summaryStatistics"]).rstrip("/") + "/harmonised/"
            )
            log(f"[DOWNLOAD] Downloading GWAS from {harmonised_url}")
            mirror_info = mirror_harmonised_folder(
                harmonised_url, out_dir, preloaded_gz=winner_tmp_gz
            )
            # Winner's temp dir is now safe to remove (file was moved or unused)
            if winner_tmp_dir is not None:
                shutil.rmtree(str(winner_tmp_dir), ignore_errors=True)

        # ── Build result row ───────────────────────────────────────────────
        n_eqtl_for_mr = genotype_n if genotype_n is not None else rnaseq_n
        gwas_for_cmd = mirror_info.get("gwas_tsv") or mirror_info.get("gwas_gz") or ""

        n_gwas_val = best.get("N") if best is not None else None
        try:
            n_gwas_val = int(n_gwas_val) if n_gwas_val is not None else None
        except Exception:
            pass

        dataset_row: Dict[str, Any] = {
            "trait": trait,
            "accession": best.get("accessionId") if best is not None else None,
            "N_gwas": n_gwas_val,
            "ancestry": (
                best.get("discoverySampleAncestry") if best is not None else None
            ),
            "eqtl_biosample_type": user_biosample,
            "eqtl_path": str(eqtl_file) if eqtl_file else None,
            "n_eqtl_rnaseq": rnaseq_n,
            "n_eqtl_genotype": genotype_n,
            "n_eqtl_for_mr": n_eqtl_for_mr,
            "harmonised_url": harmonised_url,
            "harmonised_dir": mirror_info.get("harmonised_dir"),
            "gwas_path_gz": mirror_info.get("gwas_gz"),
            "gwas_path_tsv": mirror_info.get("gwas_tsv"),
            "gwas_path": str(gwas_for_cmd) if gwas_for_cmd else None,
            "dataset_folder": str(out_dir),
            "dataset_report_xlsx": None,
            "data_status": data_status,
        }

        log(f"[MASTER] Dataset recorded ({data_status}) -> {out_dir.name}")

        try:
            dataset_name = out_dir.name
            xlsx_path = write_dataset_excel(out_dir, dataset_name, dataset_row)
            dataset_row["dataset_report_xlsx"] = str(xlsx_path)
            log(f"[OK] Dataset report saved: {xlsx_path}")
        except Exception as e:
            log(f"[WARN] Failed to write dataset Excel: {e}")

        master_rows.append(dataset_row)

    # ── Write master report ────────────────────────────────────────────────
    if master_rows:
        out_master_tsv = out_root_path / "MASTER_REPORT.tsv"
        out_master_xlsx = out_root_path / "MASTER_REPORT.xlsx"
        try:
            mdf = pd.DataFrame(master_rows)
            mdf.to_csv(out_master_tsv, sep="\t", index=False)
            with pd.ExcelWriter(out_master_xlsx, engine="openpyxl") as writer:
                mdf.to_excel(writer, index=False, sheet_name="master")
            log(f"\nMASTER_REPORT saved: {out_master_tsv}")
        except Exception as e:
            log(f"[ERROR] Failed to write master reports: {e}")
    else:
        log("\n[INFO] No datasets were recorded.")

    return master_rows
