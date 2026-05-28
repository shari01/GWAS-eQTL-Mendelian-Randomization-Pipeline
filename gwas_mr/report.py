"""HTML report generator for MR pipeline results.

Produces a self-contained HTML file with embedded logos, interactive
tables, and matplotlib charts (base64-encoded).
"""

import base64
import io
import json
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LOGOS_DIR = Path(__file__).resolve().parent / "logos"


def _generate_narrative(
    disease_name: str,
    biosample_type: str,
    n_genes: int,
    n_sig: int,
    n_risk: int,
    n_protect: int,
    tier_counts: dict,
    top_hits: pd.DataFrame,
    methods_summary: dict,
) -> Optional[str]:
    """Use OpenAI to draft a short narrative summary of the MR results."""
    try:
        from decouple import config as decouple_config
        api_key = decouple_config("OPENAI_API_KEY")
    except Exception:
        return None

    top_genes_text = ""
    if top_hits is not None and not top_hits.empty:
        rows = []
        for _, r in top_hits.head(8).iterrows():
            gene = r.get("gene", "?")
            pv = r.get("pval")
            or_val = r.get("OR")
            dirn = r.get("direction", "")
            tier = r.get("confidence_tier", "")
            pv_str = f"p={pv:.2e}" if pd.notna(pv) else ""
            or_str = f"OR={or_val:.2f}" if pd.notna(or_val) else ""
            rows.append(f"  {gene}: {or_str}, {pv_str}, {dirn}, tier={tier}")
        top_genes_text = "\n".join(rows)

    tier_text = ", ".join(f"{k}: {v}" for k, v in tier_counts.items()) if tier_counts else "none"
    methods_text = ", ".join(f"{k} ({v})" for k, v in methods_summary.items()) if methods_summary else "unknown"

    prompt = f"""You are a genomics research assistant. Write a concise plain-text summary (4-6 short bullet points) of these Mendelian Randomization results. No headers, no markdown, just bullet points starting with "•". Keep each bullet to 1-2 sentences max.

Disease: {disease_name}
Biosample Type: {biosample_type}
Genes tested: {n_genes}
Significant (p<0.05): {n_sig}
Risk-increasing: {n_risk}
Protective: {n_protect}
Confidence tiers: {tier_text}
MR methods: {methods_text}
Top gene hits:
{top_genes_text}

Focus on: what was tested, key significant findings (gene names, direction, OR), overall evidence strength, and a brief clinical interpretation note. Be factual and precise."""

    try:
        from openai import OpenAI
        client = OpenAI(api_key=api_key)
        resp = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[{"role": "user", "content": prompt}],
            max_tokens=400,
            temperature=0.3,
        )
        return resp.choices[0].message.content.strip()
    except Exception:
        return None


def _img_to_base64(path: Path) -> str:
    if not path.is_file():
        return ""
    data = path.read_bytes()
    ext = path.suffix.lower()
    mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg"}.get(
        ext.lstrip("."), "image/png"
    )
    return f"data:{mime};base64,{base64.b64encode(data).decode()}"


def _fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    return f"data:image/png;base64,{base64.b64encode(buf.read()).decode()}"


def _safe_read_csv(path: Path, **kwargs) -> Optional[pd.DataFrame]:
    if not path.is_file():
        return None
    try:
        df = pd.read_csv(path, **kwargs)
        if df.empty or (len(df.columns) == 1 and df.columns[0] in ("note",)):
            return None
        return df
    except Exception:
        return None


def _find_tsv_by_suffix(directory: Path, suffix: str) -> Optional[Path]:
    for f in directory.iterdir():
        if f.name.endswith(suffix):
            return f
    return None


# ── Chart generators ───────────────────────────────────────────────────────

def _chart_volcano(mr_main: pd.DataFrame) -> Optional[str]:
    if "pval" not in mr_main.columns or "b" not in mr_main.columns:
        return None

    df = mr_main.dropna(subset=["pval", "b"]).copy()
    df = df[df["pval"] > 0]
    df["neg_log_p"] = -np.log10(df["pval"])

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#e74c3c" if p < 0.05 else "#bdc3c7" for p in df["pval"]]
    ax.scatter(df["b"], df["neg_log_p"], c=colors, alpha=0.7, s=40, edgecolors="none")
    ax.axhline(-np.log10(0.05), color="#2c3e50", ls="--", lw=0.8, alpha=0.6)
    ax.axvline(0, color="#7f8c8d", ls="-", lw=0.5, alpha=0.5)
    ax.set_xlabel("Effect size (beta)", fontsize=11)
    ax.set_ylabel("-log10(p-value)", fontsize=11)
    ax.set_title("Volcano Plot — MR Effect Sizes", fontsize=13, fontweight="bold")

    sig = df[df["pval"] < 0.05]
    if "gene" in sig.columns and len(sig) <= 15:
        for _, row in sig.iterrows():
            ax.annotate(
                row["gene"], (row["b"], row["neg_log_p"]),
                fontsize=7, ha="center", va="bottom",
                xytext=(0, 4), textcoords="offset points",
            )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    result = _fig_to_base64(fig)
    plt.close(fig)
    return result


def _chart_forest(mr_interp: pd.DataFrame) -> Optional[str]:
    needed = {"gene", "OR", "CI_low", "CI_high", "pval"}
    if not needed.issubset(mr_interp.columns):
        return None

    df = mr_interp.dropna(subset=["OR", "CI_low", "CI_high"]).copy()
    df = df.sort_values("pval")
    df = df.head(20)
    if df.empty:
        return None

    df = df.sort_values("OR")

    fig, ax = plt.subplots(figsize=(8, max(3, len(df) * 0.35 + 1)))
    y_pos = range(len(df))

    colors = []
    for _, row in df.iterrows():
        if row["pval"] < 0.05 and row["OR"] > 1:
            colors.append("#e74c3c")
        elif row["pval"] < 0.05 and row["OR"] <= 1:
            colors.append("#27ae60")
        else:
            colors.append("#7f8c8d")

    for i, (_, row) in enumerate(df.iterrows()):
        ax.plot([row["CI_low"], row["CI_high"]], [i, i], color=colors[i], lw=1.5, alpha=0.7)
        ax.scatter(row["OR"], i, color=colors[i], s=50, zorder=3, edgecolors="white", lw=0.5)

    ax.axvline(1, color="#2c3e50", ls="--", lw=0.8, alpha=0.6)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(df["gene"].tolist(), fontsize=9)
    ax.set_xlabel("Odds Ratio (95% CI)", fontsize=11)
    ax.set_title("Forest Plot — Top Genes by Significance", fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    result = _fig_to_base64(fig)
    plt.close(fig)
    return result


def _chart_tier_donut(tier_counts: dict) -> Optional[str]:
    if not tier_counts:
        return None

    order = ["Strong", "Moderate", "Weak", "Not_significant"]
    palette = {"Strong": "#27ae60", "Moderate": "#f39c12", "Weak": "#e67e22", "Not_significant": "#bdc3c7"}
    labels, sizes, colors = [], [], []
    for t in order:
        if t in tier_counts:
            labels.append(t.replace("_", " "))
            sizes.append(tier_counts[t])
            colors.append(palette.get(t, "#95a5a6"))

    if not sizes:
        return None

    fig, ax = plt.subplots(figsize=(5, 5))
    wedges, texts, autotexts = ax.pie(
        sizes, labels=labels, colors=colors, autopct="%1.0f%%",
        startangle=90, pctdistance=0.78, wedgeprops=dict(width=0.4, edgecolor="white"),
        textprops={"fontsize": 10},
    )
    for t in autotexts:
        t.set_fontsize(9)
    ax.set_title("Confidence Tier Distribution", fontsize=13, fontweight="bold", pad=16)
    plt.tight_layout()
    result = _fig_to_base64(fig)
    plt.close(fig)
    return result


def _chart_gene_evidence(gene_ev: pd.DataFrame) -> Optional[str]:
    if gene_ev is None or gene_ev.empty:
        return None
    if "Gene_Genetic_Confidence_Score" not in gene_ev.columns:
        return None

    df = gene_ev.sort_values("Gene_Genetic_Confidence_Score", ascending=True).tail(20).copy()
    strength_colors = {
        "Very_Strong": "#1a5276", "Strong": "#27ae60",
        "Moderate": "#f39c12", "Weak": "#e67e22", "Low": "#bdc3c7",
    }
    colors = [strength_colors.get(s, "#bdc3c7") for s in df.get("Evidence_Strength", [])]

    fig, ax = plt.subplots(figsize=(8, max(3, len(df) * 0.35 + 1)))
    ax.barh(range(len(df)), df["Gene_Genetic_Confidence_Score"], color=colors, edgecolor="white", height=0.7)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Gene_Symbol"].tolist(), fontsize=9)
    ax.set_xlabel("Genetic Confidence Score", fontsize=11)
    ax.set_title("Gene-Level Genetic Evidence", fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    result = _fig_to_base64(fig)
    plt.close(fig)
    return result


# ── HTML template ──────────────────────────────────────────────────────────

_CSS = """
:root {
    --primary: #1a5276;
    --accent: #e74c3c;
    --bg: #f8f9fa;
    --card-bg: #ffffff;
    --text: #2c3e50;
    --muted: #7f8c8d;
    --border: #dee2e6;
    --success: #27ae60;
    --warning: #f39c12;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: 'Segoe UI', system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); line-height: 1.6; }
.container { max-width: 1100px; margin: 0 auto; padding: 20px; }

.header {
    background: linear-gradient(135deg, var(--primary) 0%, #2c3e50 100%);
    color: white; padding: 24px 40px; display: flex;
    align-items: center; justify-content: space-between;
    border-radius: 12px; margin-bottom: 28px;
}
.header-logo { height: 70px; max-width: 180px; object-fit: contain; }
.header-center { text-align: center; flex: 1; padding: 0 20px; }
.header-center h1 { font-size: 22px; font-weight: 700; letter-spacing: 0.5px; }
.header-center .subtitle { font-size: 13px; opacity: 0.85; margin-top: 4px; }

.meta-bar {
    display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 12px; margin-bottom: 24px;
}
.meta-item {
    background: var(--card-bg); border-radius: 8px; padding: 14px 18px;
    border-left: 4px solid var(--primary); box-shadow: 0 1px 3px rgba(0,0,0,0.06);
}
.meta-item .label { font-size: 11px; text-transform: uppercase; color: var(--muted); font-weight: 600; letter-spacing: 0.5px; }
.meta-item .value { font-size: 16px; font-weight: 600; margin-top: 2px; }

.card {
    background: var(--card-bg); border-radius: 10px; padding: 24px;
    margin-bottom: 20px; box-shadow: 0 1px 4px rgba(0,0,0,0.06);
}
.card h2 { font-size: 17px; font-weight: 700; color: var(--primary); margin-bottom: 16px; padding-bottom: 8px; border-bottom: 2px solid var(--border); }

.stats-grid {
    display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
    gap: 14px; margin-bottom: 16px;
}
.stat-box {
    text-align: center; padding: 16px 10px; border-radius: 8px; background: var(--bg);
}
.stat-box .num { font-size: 28px; font-weight: 800; }
.stat-box .desc { font-size: 11px; color: var(--muted); text-transform: uppercase; letter-spacing: 0.3px; margin-top: 2px; }
.stat-box.risk .num { color: var(--accent); }
.stat-box.protect .num { color: var(--success); }
.stat-box.highlight .num { color: var(--primary); }

.chart-row { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 20px; }
.chart-row.single { grid-template-columns: 1fr; }
.chart-box { text-align: center; }
.chart-box img { max-width: 100%; border-radius: 8px; border: 1px solid var(--border); }

table { width: 100%; border-collapse: collapse; font-size: 13px; }
th { background: var(--primary); color: white; padding: 10px 12px; text-align: left; font-weight: 600; font-size: 11px; text-transform: uppercase; letter-spacing: 0.3px; }
td { padding: 8px 12px; border-bottom: 1px solid var(--border); }
tr:hover td { background: #eef2f7; }
.badge { display: inline-block; padding: 2px 10px; border-radius: 12px; font-size: 11px; font-weight: 600; }
.badge-strong { background: #d4efdf; color: #1e8449; }
.badge-moderate { background: #fdebd0; color: #b9770e; }
.badge-weak { background: #fce4d6; color: #c0392b; }
.badge-ns { background: #eaecee; color: #7f8c8d; }
.badge-risk { background: #fadbd8; color: #c0392b; }
.badge-protect { background: #d5f5e3; color: #1e8449; }

.footer { text-align: center; padding: 20px; color: var(--muted); font-size: 12px; }

@media print {
    body { background: white; }
    .container { max-width: 100%; }
    .card { break-inside: avoid; box-shadow: none; border: 1px solid var(--border); }
}
@media (max-width: 768px) {
    .chart-row { grid-template-columns: 1fr; }
    .header { flex-direction: column; text-align: center; gap: 12px; }
}
"""


def _tier_badge(tier: str) -> str:
    cls = {
        "Strong": "badge-strong", "Moderate": "badge-moderate",
        "Weak": "badge-weak",
    }.get(tier, "badge-ns")
    return f'<span class="badge {cls}">{tier.replace("_", " ")}</span>'


def _direction_badge(d: str) -> str:
    if not isinstance(d, str):
        return ""
    if "higher risk" in d.lower():
        return '<span class="badge badge-risk">Risk</span>'
    if "protective" in d.lower():
        return '<span class="badge badge-protect">Protective</span>'
    return ""


def _df_to_html_table(df: pd.DataFrame, columns: list, formatters: dict = None) -> str:
    formatters = formatters or {}
    rows = []
    rows.append("<table><thead><tr>" +
                "".join(f"<th>{c}</th>" for c in columns) +
                "</tr></thead><tbody>")
    for _, row in df.iterrows():
        cells = []
        for c in columns:
            val = row.get(c, "")
            if c in formatters:
                cells.append(f"<td>{formatters[c](val, row)}</td>")
            elif isinstance(val, float):
                if abs(val) < 0.001 or abs(val) > 10000:
                    cells.append(f"<td>{val:.2e}</td>")
                else:
                    cells.append(f"<td>{val:.4f}</td>")
            else:
                cells.append(f"<td>{val}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    rows.append("</tbody></table>")
    return "\n".join(rows)


# ── Main generator ─────────────────────────────────────────────────────────

def generate_report(
    mr_output_dir: str,
    disease_name: str,
    biosample_type: str,
    gwas_accession: Optional[str] = None,
    n_gwas: Optional[int] = None,
    n_eqtl: Optional[int] = None,
) -> Optional[Path]:
    """Generate a self-contained HTML summary report.

    Returns the path to the written ``.html`` file, or ``None`` if
    critical output files are missing.
    """
    base = Path(mr_output_dir)
    d12 = base / "summary_tables"
    d00 = base / "00_config"

    if not d12.is_dir():
        return None

    # ── Load data ──────────────────────────────────────────────────────
    mr_main = _safe_read_csv(d12 / "MR_MAIN_RESULTS_ALL_GENES.csv")
    mr_interp = _safe_read_csv(d12 / "MR_INTERPRETATION_SUMMARY.csv")
    gene_ev_path = _find_tsv_by_suffix(d12, "_GeneLevel_GeneticEvidence.tsv")
    gene_ev = _safe_read_csv(gene_ev_path, sep="\t") if gene_ev_path else None
    coloc = _safe_read_csv(d12 / "COLOC_SUMMARY.csv")

    config = None
    if (d00 / "run_config.json").is_file():
        try:
            config = json.loads((d00 / "run_config.json").read_text())
        except Exception:
            pass

    if mr_main is None:
        return None

    # ── Compute stats ──────────────────────────────────────────────────
    n_genes = mr_main["gene"].nunique() if "gene" in mr_main.columns else len(mr_main)
    sig = mr_main[mr_main["pval"] < 0.05] if "pval" in mr_main.columns else pd.DataFrame()
    n_sig = sig["gene"].nunique() if "gene" in sig.columns and not sig.empty else 0

    n_risk = 0
    n_protect = 0
    if not sig.empty and "direction" in sig.columns:
        n_risk = sig[sig["direction"].str.contains("higher risk", case=False, na=False)]["gene"].nunique()
        n_protect = sig[sig["direction"].str.contains("protective", case=False, na=False)]["gene"].nunique()

    tier_counts = {}
    if mr_interp is not None and "confidence_tier" in mr_interp.columns:
        tier_counts = mr_interp["confidence_tier"].value_counts().to_dict()

    n_methods = mr_main["method"].nunique() if "method" in mr_main.columns else 0
    methods_summary = mr_main["method"].value_counts().to_dict() if "method" in mr_main.columns else {}

    # ── LLM narrative ──────────────────────────────────────────────────
    top_for_narrative = None
    if mr_interp is not None and not mr_interp.empty:
        tier_order = {"Strong": 0, "Moderate": 1, "Weak": 2, "Not_significant": 3}
        _ranked = mr_interp.copy()
        if "confidence_tier" in _ranked.columns:
            _ranked["_r"] = _ranked["confidence_tier"].map(tier_order).fillna(4)
            _ranked = _ranked.sort_values(["_r", "pval"] if "pval" in _ranked.columns else ["_r"])
        top_for_narrative = _ranked.head(8)

    narrative = _generate_narrative(
        disease_name=disease_name,
        biosample_type=biosample_type,
        n_genes=n_genes,
        n_sig=n_sig,
        n_risk=n_risk,
        n_protect=n_protect,
        tier_counts=tier_counts,
        top_hits=top_for_narrative,
        methods_summary=methods_summary,
    )

    # ── Logos ──────────────────────────────────────────────────────────
    logo_left = _img_to_base64(LOGOS_DIR / "ARI-Ayass-Research-Institute.png")
    logo_right = _img_to_base64(LOGOS_DIR / "Laboratory-Testing-RF.jpg")

    # ── Charts ─────────────────────────────────────────────────────────
    volcano_img = _chart_volcano(mr_main)
    forest_img = _chart_forest(mr_interp) if mr_interp is not None else None
    tier_img = _chart_tier_donut(tier_counts)
    gene_ev_img = _chart_gene_evidence(gene_ev)

    # ── Build HTML ─────────────────────────────────────────────────────
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    parts = []
    parts.append(f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>MR Report — {disease_name}</title>
<style>{_CSS}</style>
</head>
<body>
<div class="container">

<div class="header">
  <img src="{logo_left}" class="header-logo" alt="ARI">
  <div class="header-center">
    <h1>Mendelian Randomization Analysis Report</h1>
    <div class="subtitle">{disease_name} &mdash; {biosample_type}</div>
  </div>
  <img src="{logo_right}" class="header-logo" alt="Ayass Bioscience">
</div>
""")

    # ── Meta bar ───────────────────────────────────────────────────────
    meta_items = [
        ("Disease / Outcome", disease_name),
        ("Biosample Type", biosample_type),
    ]
    if gwas_accession:
        meta_items.append(("GWAS Accession", gwas_accession))
    if n_gwas:
        meta_items.append(("GWAS N", f"{int(n_gwas):,}"))
    if n_eqtl:
        meta_items.append(("eQTL N", f"{int(n_eqtl):,}"))

    parts.append('<div class="meta-bar">')
    for label, value in meta_items:
        parts.append(f'<div class="meta-item"><div class="label">{label}</div><div class="value">{value}</div></div>')
    parts.append("</div>")

    # ── Overview stats ─────────────────────────────────────────────────
    parts.append("""<div class="card"><h2>Overview</h2><div class="stats-grid">""")
    parts.append(f'<div class="stat-box highlight"><div class="num">{n_genes}</div><div class="desc">Genes Tested</div></div>')
    parts.append(f'<div class="stat-box highlight"><div class="num">{n_sig}</div><div class="desc">Significant (p&lt;0.05)</div></div>')
    parts.append(f'<div class="stat-box risk"><div class="num">{n_risk}</div><div class="desc">Risk Genes</div></div>')
    parts.append(f'<div class="stat-box protect"><div class="num">{n_protect}</div><div class="desc">Protective Genes</div></div>')

    for tier in ["Strong", "Moderate", "Weak"]:
        if tier in tier_counts:
            parts.append(f'<div class="stat-box"><div class="num">{tier_counts[tier]}</div><div class="desc">{tier} Confidence</div></div>')

    parts.append("</div></div>")

    # ── Narrative summary ──────────────────────────────────────────────
    if narrative:
        bullet_html = ""
        for line in narrative.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith("•"):
                line = line[1:].strip()
            bullet_html += f"<li>{line}</li>\n"
        if bullet_html:
            parts.append(f'<div class="card"><h2>Summary of Findings</h2><ul style="margin-left:16px;line-height:1.8;">{bullet_html}</ul></div>')

    # ── Charts: Volcano + Tier ─────────────────────────────────────────
    if volcano_img or tier_img:
        if volcano_img and tier_img:
            parts.append('<div class="chart-row">')
            parts.append(f'<div class="chart-box"><img src="{volcano_img}" alt="Volcano Plot"></div>')
            parts.append(f'<div class="chart-box"><img src="{tier_img}" alt="Tier Distribution"></div>')
            parts.append("</div>")
        else:
            img = volcano_img or tier_img
            parts.append(f'<div class="chart-row single"><div class="chart-box"><img src="{img}" alt="Chart"></div></div>')

    # ── Charts: Forest + Gene evidence ─────────────────────────────────
    if forest_img or gene_ev_img:
        if forest_img and gene_ev_img:
            parts.append('<div class="chart-row">')
            parts.append(f'<div class="chart-box"><img src="{forest_img}" alt="Forest Plot"></div>')
            parts.append(f'<div class="chart-box"><img src="{gene_ev_img}" alt="Gene Evidence"></div>')
            parts.append("</div>")
        else:
            img = forest_img or gene_ev_img
            parts.append(f'<div class="chart-row single"><div class="chart-box"><img src="{img}" alt="Chart"></div></div>')

    # ── Top hits table ─────────────────────────────────────────────────
    if mr_interp is not None and not mr_interp.empty:
        parts.append('<div class="card"><h2>Top Gene Results</h2>')

        display_df = mr_interp.copy()
        tier_order = {"Strong": 0, "Moderate": 1, "Weak": 2, "Not_significant": 3}
        if "confidence_tier" in display_df.columns:
            display_df["_rank"] = display_df["confidence_tier"].map(tier_order).fillna(4)
            if "pval" in display_df.columns:
                display_df = display_df.sort_values(["_rank", "pval"])
            else:
                display_df = display_df.sort_values("_rank")

        cols = ["gene", "method", "nsnp", "b", "OR", "pval", "direction", "confidence_tier"]
        cols = [c for c in cols if c in display_df.columns]
        formatters = {
            "confidence_tier": lambda v, r: _tier_badge(v) if isinstance(v, str) else "",
            "direction": lambda v, r: _direction_badge(v),
        }
        parts.append(_df_to_html_table(display_df, cols, formatters))
        parts.append("</div>")

    # ── Gene evidence table ────────────────────────────────────────────
    if gene_ev is not None and not gene_ev.empty:
        parts.append('<div class="card"><h2>Gene-Level Genetic Evidence</h2>')
        ev_cols = ["Gene_Symbol", "n_eQTL_variants", "Mean_eQTL_Z", "Mean_GWAS_Z",
                   "Min_GWAS_p", "Gene_Genetic_Confidence_Score", "Evidence_Strength"]
        ev_cols = [c for c in ev_cols if c in gene_ev.columns]
        ev_sorted = gene_ev.sort_values("Gene_Genetic_Confidence_Score", ascending=False)

        def _ev_badge(v, r):
            cls_map = {"Very_Strong": "badge-strong", "Strong": "badge-strong",
                       "Moderate": "badge-moderate", "Weak": "badge-weak", "Low": "badge-ns"}
            cls = cls_map.get(v, "badge-ns")
            return f'<span class="badge {cls}">{str(v).replace("_"," ")}</span>'

        parts.append(_df_to_html_table(ev_sorted, ev_cols, {"Evidence_Strength": _ev_badge}))
        parts.append("</div>")

    # ── Colocalization ─────────────────────────────────────────────────
    parts.append('<div class="card"><h2>Colocalization Summary</h2>')
    if coloc is not None and "PP4_shared" in coloc.columns:
        strong_coloc = coloc[coloc["PP4_shared"] >= 0.75]
        parts.append(f"<p>Genes tested: <strong>{len(coloc)}</strong> &nbsp;|&nbsp; Strong coloc (PP4 &ge; 0.75): <strong>{len(strong_coloc)}</strong></p>")
        if not strong_coloc.empty:
            coloc_cols = [c for c in ["gene", "PP4_shared", "PP3_distinct", "nsnp_region_common"] if c in coloc.columns]
            parts.append(_df_to_html_table(strong_coloc, coloc_cols))
    else:
        parts.append("<p>No colocalization results available (requires &ge; 2 shared SNPs per gene).</p>")
    parts.append("</div>")

    # ── Methods & params ───────────────────────────────────────────────
    parts.append('<div class="card"><h2>Methods &amp; Parameters</h2>')
    if methods_summary:
        parts.append("<p><strong>MR methods applied:</strong></p><ul>")
        for m, c in sorted(methods_summary.items(), key=lambda x: -x[1]):
            parts.append(f"<li>{m}: {c} tests</li>")
        parts.append("</ul>")
    if config:
        parts.append("<p><strong>Key parameters:</strong></p><ul>")
        param_keys = [
            ("EQTL_P_THRESH", "eQTL p-value threshold"),
            ("FSTAT_THRESH", "F-statistic threshold"),
            ("CLUMP_R2_STRICT", "Clumping r\u00b2 (strict)"),
            ("CLUMP_R2_SMALLN", "Clumping r\u00b2 (small N)"),
            ("CLUMP_KB", "Clumping window (kb)"),
            ("HARMONISE_ACTION", "Harmonise action"),
            ("N_EQTL_DEFAULT", "N eQTL"),
            ("N_GWAS_DEFAULT", "N GWAS"),
        ]
        for key, label in param_keys:
            if key in config:
                parts.append(f"<li>{label}: <code>{config[key]}</code></li>")
        parts.append("</ul>")
    parts.append("</div>")

    # ── Footer ─────────────────────────────────────────────────────────
    parts.append(f"""
<div class="footer">
  Generated by <strong>gwas_mr</strong> pipeline &mdash; {now}<br>
  Ayass Research Institute &amp; Ayass Bioscience Laboratory Testing
</div>

</div>
</body>
</html>""")

    html = "\n".join(parts)

    report_path = base / "MR_PIPELINE_REPORT.html"
    report_path.write_text(html, encoding="utf-8")
    return report_path
