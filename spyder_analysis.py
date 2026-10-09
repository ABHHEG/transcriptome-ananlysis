"""Descriptive time-course analysis and plots for configurable expression data."""
from pathlib import Path
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def analyze_time_course(data, timepoint_columns, baseline_label, threshold, pseudocount=0.01):
    """Compare one selected time point with every other point and summarize trends."""
    labels = list(timepoint_columns)
    if baseline_label not in timepoint_columns:
        raise ValueError(f"Unknown baseline time point: {baseline_label}")

    baseline_col = timepoint_columns[baseline_label]
    contrast_frames = []
    for comparison_label in labels:
        if comparison_label == baseline_label:
            continue
        comparison_col = timepoint_columns[comparison_label]
        baseline_values = pd.to_numeric(data[baseline_col], errors="coerce")
        comparison_values = pd.to_numeric(data[comparison_col], errors="coerce")
        fold_change = (comparison_values + pseudocount) / (baseline_values + pseudocount)
        log2fc = np.log2(fold_change)
        regulation = np.select(
            [log2fc >= threshold, log2fc <= -threshold],
            ["Upregulated", "Downregulated"],
            default="Below threshold",
        )
        contrast_frames.append(pd.DataFrame({
            "Gene_ID": data["Gene_ID"].astype(str),
            "Baseline": baseline_label,
            "Comparison": comparison_label,
            "Baseline_Expression": baseline_values,
            "Comparison_Expression": comparison_values,
            "Fold_Change": fold_change,
            "Log2_Fold_Change": log2fc,
            "Regulation": regulation,
            "P_Value": np.nan,
            "Adjusted_P_Value": np.nan,
        }))

    pairwise = pd.concat(contrast_frames, ignore_index=True)
    expression_columns = [timepoint_columns[label] for label in labels]
    log_expression = np.log2(data[expression_columns].apply(pd.to_numeric, errors="coerce") + pseudocount)
    net_change = log_expression.iloc[:, -1] - log_expression.iloc[:, 0]
    max_change = log_expression.sub(log_expression.iloc[:, 0], axis=0).max(axis=1)
    min_change = log_expression.sub(log_expression.iloc[:, 0], axis=0).min(axis=1)

    def classify_trend(index):
        profile = log_expression.iloc[index].to_numpy(dtype=float)
        if not np.isfinite(profile).all():
            return "Insufficient data"
        differences = np.diff(profile)
        required_steps = math.ceil(0.75 * len(differences))
        mostly_up = int(np.sum(differences >= -0.05)) >= required_steps
        mostly_down = int(np.sum(differences <= 0.05)) >= required_steps
        change = float(net_change.iloc[index])
        if change >= threshold and mostly_up:
            return "Upregulated"
        if change <= -threshold and mostly_down:
            return "Downregulated"
        if max_change.iloc[index] >= threshold and min_change.iloc[index] <= -threshold:
            return "Mixed trend"
        if max_change.iloc[index] >= threshold:
            return "Transient increase"
        if min_change.iloc[index] <= -threshold:
            return "Transient decrease"
        return "No strong trend"

    trends = [classify_trend(i) for i in range(len(data))]
    strongest_contrast = pairwise.groupby("Gene_ID")["Log2_Fold_Change"].apply(lambda x: x.abs().max())
    has_threshold_change = strongest_contrast.reindex(data["Gene_ID"].astype(str)).to_numpy() >= threshold
    trend_summary = pd.DataFrame({
        "Gene_ID": data["Gene_ID"].astype(str).to_numpy(),
        "Log2_Fold_Change": net_change.to_numpy(),
        "Regulation": trends,
        "Trend": trends,
        "Any_Contrast_Above_Threshold": has_threshold_change,
        "Mean_Expression": data[expression_columns].mean(axis=1).to_numpy(),
        "P_Value": np.nan,
        "Adjusted_P_Value": np.nan,
    })
    candidates = trend_summary.loc[trend_summary["Any_Contrast_Above_Threshold"]].copy()
    candidates = candidates.merge(data[["Gene_ID"] + expression_columns], on="Gene_ID", how="left")
    return pairwise, trend_summary, candidates, log_expression


def save_time_course_outputs(pairwise, trend_summary, candidates, log_expression,
                             gene_ids, timepoint_columns, output_dir, expression_unit="expression", max_heatmap_genes=75):
    """Save comparison tables, an expression heatmap, and volcano availability note."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pairwise.to_csv(output_dir / "timepoint_comparisons.csv", index=False)
    trend_summary.to_csv(output_dir / "trend_summary_all_genes.csv", index=False)
    candidates.to_csv(output_dir / "candidate_genes.csv", index=False)

    if len(candidates):
        trend_order = {"Upregulated": 0, "Downregulated": 1, "Transient increase": 2, "Transient decrease": 3, "Mixed trend": 4, "No strong trend": 5, "Insufficient data": 6}
        heatmap_candidates = candidates.assign(_trend_order=candidates["Trend"].map(trend_order).fillna(99)).sort_values(["_trend_order", "Gene_ID"])
        candidate_ids = heatmap_candidates["Gene_ID"].astype(str).tolist()
        positions = {str(gene): i for i, gene in enumerate(gene_ids)}
        rows = [positions[gene] for gene in candidate_ids if gene in positions]
        matrix = log_expression.iloc[rows].to_numpy(dtype=float)
        row_means = np.nanmean(matrix, axis=1, keepdims=True)
        row_stds = np.nanstd(matrix, axis=1, keepdims=True)
        row_stds[row_stds == 0] = 1.0
        z_matrix = (matrix - row_means) / row_stds
        labels = [str(gene_ids[i]) for i in rows]
        if len(labels) > max_heatmap_genes:
            priority = np.nanmax(matrix, axis=1) - np.nanmin(matrix, axis=1)
            keep = sorted(np.argsort(priority)[-max_heatmap_genes:])
            z_matrix = z_matrix[keep]
            labels = [labels[i] for i in keep]

        height = max(5, min(24, 0.24 * len(labels) + 2))
        fig, ax = plt.subplots(figsize=(10, height))
        image = ax.imshow(z_matrix, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-2, vmax=2)
        ax.set_xticks(range(len(timepoint_columns)))
        ax.set_xticklabels(list(timepoint_columns), rotation=0)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_title(f"Candidate gene expression trends (top {len(labels)} genes)")
        ax.set_xlabel("Time point")
        ax.set_ylabel("Gene_ID")
        fig.colorbar(image, ax=ax, label=f"Per-gene z-score of log₂({expression_unit} + 0.01)")
        fig.tight_layout()
        fig.savefig(output_dir / "expression_heatmap.png", dpi=180, bbox_inches="tight")
        plt.show()
    else:
        print("No genes meet the fold-change threshold; heatmap not generated.")

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.axis("off")
    ax.text(0.5, 0.60, "A statistical volcano plot cannot be calculated from this input.",
            ha="center", va="center", fontsize=13, weight="bold", wrap=True)
    ax.text(0.5, 0.36,
            "This workflow calculates descriptive fold-change screens, not replicate-aware p-values/FDR. "
            "A statistical volcano plot needs replicate-level data analyzed with an appropriate model. "
            "If replicate sample columns were supplied, they were averaged for the time-course screen.",
            ha="center", va="center", fontsize=10, wrap=True)
    fig.tight_layout()
    fig.savefig(output_dir / "volcano_plot_unavailable.png", dpi=180, bbox_inches="tight")
    plt.show()


def render_colored_kegg_map(url, output_dir):
    """Fetch the KEGG PNG/KGML, overlay mapped gene colors, and show it in Spyder."""
    from io import BytesIO
    from urllib.parse import parse_qs, urlparse
    import xml.etree.ElementTree as ET
    import requests
    from matplotlib.patches import Rectangle

    query = parse_qs(urlparse(url).query)
    pathway_id = query.get("map", [""])[0]
    if not pathway_id:
        raise ValueError("Could not determine the KEGG pathway identifier from the URL.")

    colors = {}
    for line in query.get("multi_query", [""])[0].splitlines():
        fields = line.split()
        if len(fields) >= 2:
            color_values = fields[1].split(",")
            colors[fields[0]] = (color_values[0], color_values[-1])

    base_url = "https://rest.kegg.jp/get/" + pathway_id
    image_response = requests.get(base_url + "/image", timeout=60)
    image_response.raise_for_status()
    image = plt.imread(BytesIO(image_response.content), format="png")

    kgml_response = requests.get(base_url + "/kgml", timeout=60)
    kgml_response.raise_for_status()
    root = ET.fromstring(kgml_response.content)

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.imshow(image, origin="upper")
    image_height, image_width = image.shape[:2]
    colored_boxes = 0
    for entry in root.findall("entry"):
        if entry.get("type") != "gene":
            continue
        matched_color = None
        for gene_id in entry.get("name", "").split():
            if gene_id in colors:
                matched_color = colors[gene_id]
                break
        if not matched_color:
            continue
        background, foreground = matched_color
        for graphic in entry.findall("graphics"):
            try:
                center_x = float(graphic.get("x"))
                center_y = float(graphic.get("y"))
                width = float(graphic.get("width", 46))
                height = float(graphic.get("height", 17))
            except (TypeError, ValueError):
                continue
            ax.add_patch(Rectangle(
                (center_x - width / 2, center_y - height / 2),
                width,
                height,
                facecolor=background,
                edgecolor=foreground,
                linewidth=1.5,
                alpha=0.62,
                zorder=3,
            ))
            colored_boxes += 1

    ax.set_xlim(0, image_width)
    ax.set_ylim(image_height, 0)
    ax.axis("off")
    ax.set_title(f"{pathway_id}: mapped expression trends ({colored_boxes} gene boxes colored)")
    fig.tight_layout(pad=0.2)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{pathway_id}_colored_pathway.png", dpi=160, bbox_inches="tight")
    plt.show()
    print(f"Colored KEGG pathway image saved to {output_dir / f'{pathway_id}_colored_pathway.png'}")
    if colors and not colored_boxes:
        print("No dataset genes matched KGML gene entries; the uncolored reference map is shown.")
    return fig




