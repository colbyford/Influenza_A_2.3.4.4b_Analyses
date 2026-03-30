"""
Paired slopegraph — one column per WT/mutant pair, shared y-axis.
Binding metric: van der Waals (vdW) energy.
Data: single CSV file; pairs defined explicitly in CONFIG.

──────────────────────────────────────────────────────────────────
EXPECTED CSV FORMAT
──────────────────────────────────────────────────────────────────
  data.csv
  ┌──────────────┬─────────┬────────┐
  │ accession    │ type    │ vdw    │
  ├──────────────┼─────────┼────────┤
  │ 2FK0_avian   │ WT      │ -28.19 │
  │ 2FK0_avian   │ mutant  │ -26.68 │
  │ 2FK0_mammal  │ WT      │ -37.52 │
  │ 2FK0_mammal  │ mutant  │ -21.95 │
  │ 2HU4_NA      │ WT      │ -31.10 │
  │ 3CL0_NA      │ mutant  │ -29.92 │
  └──────────────┴─────────┴────────┘
──────────────────────────────────────────────────────────────────
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

# ── CONFIG ────────────────────────────────────────────────────────────────────
CSV_FILE    = "final_output_w_9NR.csv"          # path to your CSV

COL_ACC     = "accession"         # column: accession/sample name
COL_TYPE    = "type"              # column: condition label
COL_VDW     = "vdw"              # column: vdW energy value

LABEL_WT    = "WT"                # value in COL_TYPE that means wild-type
LABEL_MUT   = "mutant"           # value in COL_TYPE that means mutant

# Define pairs explicitly as (pair_label, wt_accession, mutant_accession).
# This handles cases where WT and mutant have different accession names.
## 2FK0
# PAIRS = [
#     ("HA_Avian",   "2FK0_avian",   "2FK0_avian"),
#     ("HA_Mammal",  "2FK0_mammal",  "2FK0_mammal"),
#     ("NA", "2HU4_NA",      "3CL0_NA"),
# ]

## 9NR*
PAIRS = [
    ("HA_Avian",   "9NR2_avian",   "9NR5_avian"),
    ("HA_Mammal",  "9NR2_mammal",  "9NR5_mammal"),
    ("NA", "2HU4_NA",      "3CL0_NA"),
]

YLABEL      = "van der Waals Energy (kcal/mol)"

# Colors
COLOR_WT     = "#2C3E50"   # dark slate — WT dots
COLOR_WEAK   = "#C0392B"   # red        — mutant less negative (weaker interaction)
COLOR_STRONG = "#2980B9"   # blue       — mutant more negative (stronger interaction)
COLOR_LINE   = "#AAAAAA"   # grey       — connecting line
# ─────────────────────────────────────────────────────────────────────────────

# ── LOAD DATA ─────────────────────────────────────────────────────────────────
df = pd.read_csv(CSV_FILE)

def get_vdw(accession, condition):
    row = df[(df[COL_ACC] == accession) & (df[COL_TYPE] == condition)]
    if row.empty:
        raise ValueError(f"No row found for accession='{accession}', type='{condition}'")
    return row[COL_VDW].values[0]

# Resolve vdW values for each pair
pair_data = []
for label, wt_acc, mut_acc in PAIRS:
    wt_val  = get_vdw(wt_acc,  LABEL_WT)
    mut_val = get_vdw(mut_acc, LABEL_MUT)
    pair_data.append((label, wt_val, mut_val))

n_pairs = len(pair_data)

# ── FIGURE LAYOUT ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(
    1, n_pairs,
    figsize=(2.2 * n_pairs, 5),
    sharey=True,
)
if n_pairs == 1:
    axes = [axes]

x_wt  = 0
x_mut = 1

for ax, (name, wt, mut) in zip(axes, pair_data):
    mut_color = COLOR_WEAK if mut > wt else COLOR_STRONG

    # Connecting line
    ax.plot(
        [x_wt, x_mut], [wt, mut],
        color=COLOR_LINE, linewidth=1.2, zorder=1, alpha=0.85,
    )

    # WT dot
    ax.scatter(x_wt,  wt,  color=COLOR_WT,   s=70, zorder=3, linewidths=0)
    # Mutant dot
    ax.scatter(x_mut, mut, color=mut_color,  s=70, zorder=3, linewidths=0)

    # Pair name as column title
    ax.set_title(name, fontsize=9, fontweight="bold", pad=6, color="#222222")

    # X-axis
    ax.set_xticks([x_wt, x_mut])
    ax.set_xticklabels(["WT", "Mut"], fontsize=8.5, fontweight="bold")
    ax.set_xlim(-0.5, 1.5)
    ax.tick_params(axis="x", length=0)

    # Minimalist spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

# ── Y-AXIS: leftmost panel only ───────────────────────────────────────────────
axes[0].set_ylabel(YLABEL, fontsize=9, labelpad=8)
axes[0].tick_params(axis="y", labelsize=8)
axes[0].yaxis.set_ticks_position("left")
axes[0].spines["left"].set_visible(True)

for ax in axes[1:]:
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)

# ── LEGEND ────────────────────────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(color=COLOR_WT,     label="Wild type"),
    mpatches.Patch(color=COLOR_STRONG, label="Mutation — stronger (more negative)"),
    mpatches.Patch(color=COLOR_WEAK,   label="Mutation — weaker (less negative)"),
]
fig.legend(
    handles=legend_handles,
    frameon=False,
    fontsize=7.5,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.08),
    ncol=3,
)

# ── SAVE ──────────────────────────────────────────────────────────────────────
plt.tight_layout()
plt.savefig("binding_validation_slopegraph.pdf", dpi=300, bbox_inches="tight")
plt.show()
print("Saved: binding_validation_slopegraph.pdf")