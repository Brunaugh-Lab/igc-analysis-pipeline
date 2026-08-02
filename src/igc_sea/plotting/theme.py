"""Journal-quality matplotlib theme for chromatography figures.

Specs:
  - Single column: 8.5 cm (3.35 in)
  - Double column: 17.4 cm (6.85 in)
  - Font: 8–10 pt sans-serif (Arial)
  - Export as PDF (vector) at 300+ dpi for raster elements
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Figure widths in inches
SINGLE_COL = 3.35  # 8.5 cm
DOUBLE_COL = 6.85  # 17.4 cm

# Consistent color palette across all modules
PALETTE = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # teal
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # pink
    "#999999",  # grey
]

# Component colors (for scatter plots vs composition)
COMPONENT_COLORS = {
    "x1": "#E69F00",
    "x2": "#56B4E9",
    "x3": "#009E73",
}


def apply_journal_theme(font_size: int = 9) -> None:
    """Apply publication-quality theme to all subsequent matplotlib figures.

    Parameters
    ----------
    font_size : int
        Base font size in points (default 9, acceptable range 8–10).
    """
    sns.set_theme(style="ticks", font_scale=1.0)

    mpl.rcParams.update({
        # Font
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": font_size,
        "axes.titlesize": font_size + 1,
        "axes.labelsize": font_size,
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,
        "legend.fontsize": font_size - 1,
        # Lines
        "lines.linewidth": 1.0,
        "lines.markersize": 4,
        # Axes
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        # Ticks
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.minor.visible": False,
        "ytick.minor.visible": False,
        # Grid
        "axes.grid": False,
        # Figure
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        # PDF export
        "pdf.fonttype": 42,  # TrueType — editable in Illustrator
        "ps.fonttype": 42,
    })


def get_palette(n: int | None = None) -> list[str]:
    """Return the standard color palette.

    Parameters
    ----------
    n : int, optional
        Number of colors to return. If None, returns the full palette.
    """
    if n is None:
        return PALETTE.copy()
    return PALETTE[:n]


def save_figure(fig: plt.Figure, path, width: float = DOUBLE_COL,
                height: float | None = None) -> None:
    """Save a figure as PDF at journal dimensions.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    path : str or Path
        Output file path (should end in .pdf).
    width : float
        Figure width in inches (default: double column).
    height : float, optional
        Figure height in inches. If None, uses current aspect ratio.
    """
    if height is None:
        w_old, h_old = fig.get_size_inches()
        height = h_old * (width / w_old)
    fig.set_size_inches(width, height)
    fig.savefig(path, format="pdf")
