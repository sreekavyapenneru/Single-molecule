import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d

# ========== USER INPUT ==========
sequence = """EQSLDQLVEEDKKFLDALNENFKDNPYSSGNCNDDCNDTCVDNQYLYDDITDDIMCAKAIIESDKALGQICGKKLSEKDLKFLKDEYLKHLSDEKYNLIDQWLCEKL"""
window_size = 5

# separate smoothing parameters
smooth_window_top = 4
smooth_window_toc = 5
gauss_sigma_top = 2
gauss_sigma_toc = 1
# =================================

# Residue groups
positive_residues = set(['R','K','H'])
negative_residues = set(['D','E'])
aromatic_residues = set(['F','W','Y'])
helix_breaker_residues = set(['G','P'])

# ---------- sliding-window raw values ----------
pos_vals, neg_vals, arom_vals, hb_vals = [], [], [], []
for i in range(len(sequence) - window_size + 1):
    window = sequence[i:i+window_size]
    pos_vals.append(sum(aa in positive_residues for aa in window))
    neg_vals.append(sum(aa in negative_residues for aa in window))
    arom_vals.append(sum(aa in aromatic_residues for aa in window))
    hb_vals.append(sum(aa in helix_breaker_residues for aa in window) >= 2)

pos_vals = np.array(pos_vals)
neg_vals = np.array(neg_vals)
arom_vals = np.array(arom_vals)
hb_vals = np.array(hb_vals, dtype=float)

# ---------- align to full seq length ----------
def align(values):
    pad_left = window_size // 2
    pad_right = len(sequence) - len(values) - pad_left
    return np.pad(values, (pad_left, pad_right), mode='edge')

pos_full = align(pos_vals)
neg_full = align(neg_vals)
arom_full = align(arom_vals)
hb_full = align(hb_vals)

x = np.arange(len(sequence))

# ---------- smoothing ----------
def smooth_box(data, w):
    if w <= 1:
        return data
    return np.convolve(data, np.ones(w)/w, mode='same')

def smooth_apply(data, w, sigma):
    out = smooth_box(data, w)
    if sigma > 0:
        out = gaussian_filter1d(out, sigma=sigma)
    return out

# smoothed top panel
pos_s = smooth_apply(pos_full, smooth_window_top, gauss_sigma_top)
neg_s = smooth_apply(neg_full, smooth_window_top, gauss_sigma_top)
arom_s = smooth_apply(arom_full, smooth_window_top, gauss_sigma_top)
hb_s = smooth_apply(hb_full, smooth_window_top, gauss_sigma_top)

# toc score
toc_raw = (3*arom_full + pos_full + 3*hb_full)/(1 + 3*neg_full)
toc_s = smooth_apply(toc_raw, smooth_window_toc, gauss_sigma_toc)

# GP positions
gp_positions = np.where(hb_full >= 1)[0]


# ============================================================
# ============   FIGURE 1 (TOP PANEL + SEQUENCE)   ============
# ============================================================
fig1, axes1 = plt.subplots(2, 1, figsize=(16, 5),
                           gridspec_kw={'height_ratios': [4, 1]},
                           sharex=True)

# ---- PANEL 1: Top feature panel ----
ax = axes1[0]
ax.plot(x, arom_s, color='tab:blue', linewidth=2, label='Aromatic')
ax.fill_between(x, arom_s, color='tab:blue', alpha=0.25)

ax.plot(x, pos_s, color='tab:green', linewidth=2, label='Positive')
ax.fill_between(x, pos_s, color='tab:green', alpha=0.25)

ax.plot(x, neg_s, color='tab:red', linewidth=2, label='Negative')
ax.fill_between(x, neg_s, color='tab:red', alpha=0.25)

ax.scatter(gp_positions, np.zeros(len(gp_positions)),
           color='black', s=60, zorder=5, label='GP-rich')

ax.set_ylim(-0.1, 2.5)
ax.set_ylabel("Value")
ax.set_title("Sequence Features of RCMLa")
ax.legend(loc="upper right")

# ---- PANEL 2: Sequence panel ----
ax3 = axes1[1]
ax3.set_ylim(-0.1, 1)
ax3.set_yticks([])

left = -0.4
right = len(sequence) - 0.6
ax3.fill_between([left, right], 0.15, 0.85,
                 color='white', edgecolor='black', linewidth=0.7)

for i, aa in enumerate(sequence):
    ax3.text(i, 0.5, aa, ha="center", va="center", fontsize=9)

ax3.set_xlim(left, right)
ax3.set_xlabel("Sequence Position")

plt.tight_layout()
plt.show()



# ============================================================
# ============   FIGURE 2 (TOC SCORE + SEQUENCE)   ============
# ============================================================
fig2, axes2 = plt.subplots(2, 1, figsize=(16, 2),
                           gridspec_kw={'height_ratios': [3, 1]},
                           sharex=True)

# ---- PANEL 1: TOC score ----
ax2 = axes2[0]
ax2.plot(x, toc_s, color='purple', linewidth=2, label='Toc34 Recognition Score')
ax2.fill_between(x, toc_s, color='purple', alpha=0.3)

ax2.set_ylim(-0.1, 5.5)
ax2.set_ylabel("Toc34 Recognition Score")
ax2.legend(loc='upper right')

# ---- PANEL 2: Sequence panel ----
ax4 = axes2[1]
ax4.set_ylim(0, 1)
ax4.set_yticks([])

ax4.fill_between([left, right], 0.15, 0.85,
                 color='white', edgecolor='black', linewidth=0.7)

for i, aa in enumerate(sequence):
    ax4.text(i, 0.5, aa, ha="center", va="center", fontsize=9)

ax4.set_xlim(left, right)
ax4.set_xlabel("Sequence Position")

plt.tight_layout()
plt.show()
