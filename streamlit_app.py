"""Streamlit app for visualizing 2D PES files."""

from __future__ import annotations

from pathlib import Path
import tempfile

import numpy as np
import plotly.graph_objects as go
import streamlit as st

# Ensure src is on the path when running from repo root.
import sys

SRC_ROOT = Path(__file__).resolve().parent / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from python_scripts.dynamics_1d.constants import CONST
from python_scripts.dynamics_2d.io import read_pes_file_2d


st.set_page_config(page_title="PES 2D Viewer", layout="wide")


@st.cache_data(show_spinner=False)
def load_pes_from_bytes(
    data: bytes,
    position_units: str,
    energy_units: str,
    index_order: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load PES using the same reader as dynamics_2d."""
    with tempfile.NamedTemporaryFile(suffix=".dat", delete=True) as tmp:
        tmp.write(data)
        tmp.flush()
        return read_pes_file_2d(
            Path(tmp.name),
            position_units=position_units,
            energy_units=energy_units,
            index_order=index_order,
        )


def convert_positions_for_display(x: np.ndarray, units: str) -> np.ndarray:
    if units.lower() == "angstrom":
        return x * 1e10
    if units.lower() == "bohr":
        return x / CONST.bohr
    return x


def convert_energy_for_display(E: np.ndarray, units: str) -> np.ndarray:
    if units.lower() == "hartree":
        return E / CONST.hartree
    if units.lower() == "ev":
        return E / CONST.eV
    return E


st.title("2D Potential Energy Surface Viewer")

with st.sidebar:
    st.header("Load PES")
    source_mode = st.radio(
        "Source",
        ["Upload file", "Local path"],
        index=0,
    )
    uploaded = None
    local_path = None
    if source_mode == "Upload file":
        uploaded = st.file_uploader(
            "PES file (columns: x1 x2 E)",
            type=["dat", "txt", "pes"],
        )
    else:
        local_path = st.text_input(
            "Local file path",
            value="",
            placeholder="/path/to/pes.dat",
        )
    position_units = st.selectbox(
        "Position units in file",
        ["angstrom", "bohr"],
        index=0,
    )
    energy_units = st.selectbox(
        "Energy units in file",
        ["hartree", "ev"],
        index=0,
    )
    index_order = st.selectbox(
        "Index ordering",
        ["C", "F"],
        index=0,
        help="C: x2 varies fastest, F: x1 varies fastest",
    )

    st.header("Display options")
    zero_at_min = st.checkbox("Zero energy at minimum", value=True)
    show_range = st.checkbox("Zoom to energy range", value=False)
    stride = st.select_slider(
        "Downsample factor",
        options=[1, 2, 4, 6, 8, 10],
        value=2,
        help="Use larger values for very large grids.",
    )
    contour_levels = st.select_slider(
        "Contour levels",
        options=[10, 15, 20, 30, 40, 60],
        value=30,
    )
    colormap = st.selectbox(
        "Colormap",
        ["Viridis", "Cividis", "Plasma", "Magma", "Turbo"],
        index=0,
    )

if source_mode == "Upload file":
    if not uploaded:
        st.info("Upload a PES file to get started.")
        st.stop()
else:
    if not local_path:
        st.info("Enter a local PES file path to get started.")
        st.stop()
    if not Path(local_path).is_file():
        st.error("Local file path not found.")
        st.stop()

try:
    if source_mode == "Upload file":
        x1_si, x2_si, E_si = load_pes_from_bytes(
            uploaded.getvalue(),
            position_units=position_units,
            energy_units=energy_units,
            index_order=index_order,
        )
    else:
        x1_si, x2_si, E_si = read_pes_file_2d(
            Path(local_path),
            position_units=position_units,
            energy_units=energy_units,
            index_order=index_order,
        )
except Exception as exc:  # noqa: BLE001 - surface file errors to user
    st.error(f"Failed to load PES: {exc}")
    st.stop()

x1 = convert_positions_for_display(x1_si, position_units)
x2 = convert_positions_for_display(x2_si, position_units)
E = convert_energy_for_display(E_si, energy_units)

if zero_at_min:
    E = E - float(np.min(E))

# Energy range controls (for zooming into the PES bottom)
E_min = float(np.min(E))
E_max = float(np.max(E))
if show_range:
    E_low, E_high = st.slider(
        "Energy range",
        min_value=E_min,
        max_value=E_max,
        value=(E_min, E_min + 0.1 * (E_max - E_min)),
        format="%.6g",
    )
    E_range = (E_low, E_high)
else:
    E_range = (E_min, E_max)

# Downsample for interactive plotting
x1_plot = x1[::stride]
x2_plot = x2[::stride]
E_plot = E[::stride, ::stride]

x1_label = f"x1 ({position_units})"
x2_label = f"x2 ({position_units})"
E_label = f"Energy ({energy_units})"

info_cols = st.columns(4)
info_cols[0].metric("x1 points", f"{len(x1)}")
info_cols[1].metric("x2 points", f"{len(x2)}")
info_cols[2].metric("E min", f"{np.min(E):.6g}")
info_cols[3].metric("E max", f"{np.max(E):.6g}")

col1, col2 = st.columns(2)

with col1:
    st.subheader("3D Surface")
    surface = go.Surface(
        x=x1_plot,
        y=x2_plot,
        z=E_plot.T,
        colorscale=colormap,
        cmin=E_range[0],
        cmax=E_range[1],
        showscale=True,
    )
    fig_surface = go.Figure(data=[surface])
    fig_surface.update_layout(
        scene=dict(
            xaxis_title=x1_label,
            yaxis_title=x2_label,
            zaxis_title=E_label,
            zaxis=dict(range=[E_range[0], E_range[1]]),
        ),
        margin=dict(l=0, r=0, t=30, b=0),
        height=620,
    )
    st.plotly_chart(fig_surface, use_container_width=True)

with col2:
    st.subheader("Contour")
    contour = go.Contour(
        x=x1_plot,
        y=x2_plot,
        z=E_plot.T,
        colorscale=colormap,
        zmin=E_range[0],
        zmax=E_range[1],
        contours=dict(showlabels=True, start=E_range[0], end=E_range[1], size=(E_range[1] - E_range[0]) / contour_levels),
    )
    fig_contour = go.Figure(data=[contour])
    fig_contour.update_layout(
        xaxis_title=x1_label,
        yaxis_title=x2_label,
        margin=dict(l=0, r=0, t=30, b=0),
        height=620,
    )
    st.plotly_chart(fig_contour, use_container_width=True)
