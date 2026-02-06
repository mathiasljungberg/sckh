"""Streamlit app for visualizing 2D PES and dipole surface files."""

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
from python_scripts.dynamics_2d.io import read_surface_file_2d_raw


st.set_page_config(page_title="2D Surface Viewer", layout="wide")


@st.cache_data(show_spinner=False)
def load_surface_from_bytes(
    data: bytes,
    surface_type: str,
    position_units: str,
    energy_units: str,
    index_order: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Load surface data from bytes and convert coordinates/energy to SI."""
    with tempfile.NamedTemporaryFile(suffix=".dat", delete=True) as tmp:
        tmp.write(data)
        tmp.flush()
        return load_surface_from_path(
            Path(tmp.name),
            surface_type=surface_type,
            position_units=position_units,
            energy_units=energy_units,
            index_order=index_order,
        )


@st.cache_data(show_spinner=False)
def load_surface_from_path(
    path: Path,
    surface_type: str,
    position_units: str,
    energy_units: str,
    index_order: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Load surface data from path and convert to SI where relevant."""
    if surface_type == "PES":
        x1_raw, x2_raw, values, resolved_order = read_surface_file_2d_raw(
            path,
            value_columns=[2],
            index_order=index_order,
        )
    else:
        x1_raw, x2_raw, values, resolved_order = read_surface_file_2d_raw(
            path,
            value_columns=[2, 3, 4],
            index_order=index_order,
        )

    if position_units.lower() == "angstrom":
        x1_si = x1_raw * 1e-10
        x2_si = x2_raw * 1e-10
    elif position_units.lower() == "bohr":
        x1_si = x1_raw * CONST.bohr
        x2_si = x2_raw * CONST.bohr
    else:
        raise ValueError(f"Unknown position units: {position_units}")

    if surface_type == "PES":
        if energy_units.lower() == "hartree":
            values_si = values * CONST.hartree
        elif energy_units.lower() == "ev":
            values_si = values * CONST.eV
        else:
            raise ValueError(f"Unknown energy units: {energy_units}")
    else:
        # Dipole values are kept in file units (typically atomic units).
        values_si = values

    return x1_si, x2_si, values_si, resolved_order


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


def dipole_quantity(dipole_components: np.ndarray, mode: str) -> tuple[np.ndarray, str]:
    """Convert dipole component array to a scalar field for plotting."""
    if mode == "mu_x":
        return dipole_components[0], "Dipole (a.u.)"
    if mode == "mu_y":
        return dipole_components[1], "Dipole (a.u.)"
    if mode == "mu_z":
        return dipole_components[2], "Dipole (a.u.)"
    if mode == "|mu_x|":
        return np.abs(dipole_components[0]), "|Dipole| (a.u.)"
    if mode == "|mu_y|":
        return np.abs(dipole_components[1]), "|Dipole| (a.u.)"
    if mode == "|mu_z|":
        return np.abs(dipole_components[2]), "|Dipole| (a.u.)"
    if mode == "mu_x^2 + mu_y^2 + mu_z^2":
        return np.sum(dipole_components**2, axis=0), "Dipole^2 (a.u.^2)"
    raise ValueError(f"Unknown dipole display mode: {mode}")


def cut_line(
    z: np.ndarray,
    x1: np.ndarray,
    x2: np.ndarray,
    cut_direction: str,
    fixed_index: int,
) -> tuple[np.ndarray, np.ndarray, str, float]:
    """Extract a 1D cut from a 2D scalar surface."""
    if cut_direction == "x1 at fixed x2":
        return x1, z[:, fixed_index], "x1", float(x2[fixed_index])
    return x2, z[fixed_index, :], "x2", float(x1[fixed_index])


st.title("2D Surface Viewer")

with st.sidebar:
    st.header("Load Surface")
    surface_type = st.radio(
        "Surface type",
        ["PES", "Dipole"],
        index=0,
    )
    source_mode = st.radio(
        "Source",
        ["Upload file", "Local path"],
        index=0,
    )

    uploaded = None
    local_path = None
    if source_mode == "Upload file":
        uploaded = st.file_uploader(
            "Surface file",
            type=["dat", "txt", "pes"],
            help=(
                "PES columns: x1 x2 E. Dipole columns: x1 x2 mu_x mu_y mu_z."
            ),
        )
    else:
        local_path = st.text_input(
            "Local file path",
            value="",
            placeholder="/path/to/surface.dat",
        )

    position_units = st.selectbox(
        "Position units in file",
        ["angstrom", "bohr"],
        index=0,
    )

    if surface_type == "PES":
        energy_units = st.selectbox(
            "Energy units in file",
            ["hartree", "ev"],
            index=0,
        )
    else:
        energy_units = "hartree"

    index_order = st.selectbox(
        "Index ordering",
        ["auto", "C", "F"],
        index=0,
        help="Auto infers whether x2-fast (C) or x1-fast (F) ordering was used.",
    )

    st.header("Display options")

    if surface_type == "PES":
        zero_at_min = st.checkbox("Zero value at minimum", value=True)
        show_range = st.checkbox("Zoom to energy range", value=False)
        range_slider_label = "Energy range"
    else:
        zero_at_min = False
        show_range = st.checkbox("Zoom to value range", value=True)
        range_slider_label = "Value range"

    if surface_type == "Dipole":
        dipole_mode = st.selectbox(
            "Dipole quantity",
            [
                "mu_x",
                "mu_y",
                "mu_z",
                "|mu_x|",
                "|mu_y|",
                "|mu_z|",
                "mu_x^2 + mu_y^2 + mu_z^2",
            ],
            index=0,
        )
    else:
        dipole_mode = "mu_x"

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
        st.info("Upload a surface file to get started.")
        st.stop()
else:
    if not local_path:
        st.info("Enter a local surface file path to get started.")
        st.stop()
    if not Path(local_path).is_file():
        st.error("Local file path not found.")
        st.stop()

try:
    if source_mode == "Upload file":
        x1_si, x2_si, values_si, resolved_order = load_surface_from_bytes(
            uploaded.getvalue(),
            surface_type=surface_type,
            position_units=position_units,
            energy_units=energy_units,
            index_order=index_order,
        )
    else:
        x1_si, x2_si, values_si, resolved_order = load_surface_from_path(
            Path(local_path),
            surface_type=surface_type,
            position_units=position_units,
            energy_units=energy_units,
            index_order=index_order,
        )
except Exception as exc:  # noqa: BLE001 - surface file errors to user
    st.error(f"Failed to load surface: {exc}")
    st.stop()

x1 = convert_positions_for_display(x1_si, position_units)
x2 = convert_positions_for_display(x2_si, position_units)

if surface_type == "PES":
    z = convert_energy_for_display(values_si[0], energy_units)
    z_label = f"Energy ({energy_units})"
    if zero_at_min:
        z = z - float(np.min(z))
else:
    z, z_label = dipole_quantity(values_si, dipole_mode)

if index_order.lower() == "auto":
    st.caption(f"Detected index ordering: `{resolved_order}`")

z_min = float(np.min(z))
z_max = float(np.max(z))
if show_range:
    z_low, z_high = st.slider(
        range_slider_label,
        min_value=z_min,
        max_value=z_max,
        value=(z_min, z_max),
        format="%.6g",
    )
    z_range = (z_low, z_high)
else:
    z_range = (z_min, z_max)

x1_plot = x1[::stride]
x2_plot = x2[::stride]
z_plot = z[::stride, ::stride]

x1_label = f"x1 ({position_units})"
x2_label = f"x2 ({position_units})"

info_cols = st.columns(4)
info_cols[0].metric("x1 points", f"{len(x1)}")
info_cols[1].metric("x2 points", f"{len(x2)}")
info_cols[2].metric("Min", f"{np.min(z):.6g}")
info_cols[3].metric("Max", f"{np.max(z):.6g}")

contour_size = (z_range[1] - z_range[0]) / contour_levels
if contour_size <= 0:
    contour_size = 1.0

col1, col2 = st.columns(2)

with col1:
    st.subheader("3D Surface")
    surface = go.Surface(
        x=x1_plot,
        y=x2_plot,
        z=z_plot.T,
        colorscale=colormap,
        cmin=z_range[0],
        cmax=z_range[1],
        showscale=True,
    )
    fig_surface = go.Figure(data=[surface])
    fig_surface.update_layout(
        scene=dict(
            xaxis_title=x1_label,
            yaxis_title=x2_label,
            zaxis_title=z_label,
            zaxis=dict(range=[z_range[0], z_range[1]]),
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
        z=z_plot.T,
        colorscale=colormap,
        zmin=z_range[0],
        zmax=z_range[1],
        contours=dict(
            showlabels=True,
            start=z_range[0],
            end=z_range[1],
            size=contour_size,
        ),
    )
    fig_contour = go.Figure(data=[contour])
    fig_contour.update_layout(
        xaxis_title=x1_label,
        yaxis_title=x2_label,
        margin=dict(l=0, r=0, t=30, b=0),
        height=620,
    )
    st.plotly_chart(fig_contour, use_container_width=True)

st.subheader("Surface Cuts")
cut_cols = st.columns(3)
with cut_cols[0]:
    cut_direction = st.radio(
        "Cut direction",
        ["x1 at fixed x2", "x2 at fixed x1"],
        index=0,
    )

if cut_direction == "x1 at fixed x2":
    with cut_cols[1]:
        fixed_index = st.slider("x2 index", 0, len(x2) - 1, len(x2) // 2)
    fixed_value = x2[fixed_index]
    fixed_label = f"x2 = {fixed_value:.6g} {position_units}"
else:
    with cut_cols[1]:
        fixed_index = st.slider("x1 index", 0, len(x1) - 1, len(x1) // 2)
    fixed_value = x1[fixed_index]
    fixed_label = f"x1 = {fixed_value:.6g} {position_units}"

with cut_cols[2]:
    st.write("Fixed coordinate")
    st.write(fixed_label)

line_x, line_y, cut_axis, cut_fixed = cut_line(
    z=z,
    x1=x1,
    x2=x2,
    cut_direction=cut_direction,
    fixed_index=fixed_index,
)

if cut_axis == "x1":
    line_x_label = x1_label
    line_title = f"1D Cut: {z_label} at x2 = {cut_fixed:.6g} {position_units}"
else:
    line_x_label = x2_label
    line_title = f"1D Cut: {z_label} at x1 = {cut_fixed:.6g} {position_units}"

fig_cut = go.Figure(
    data=[
        go.Scatter(
            x=line_x,
            y=line_y,
            mode="lines",
            line=dict(width=2),
        )
    ]
)
fig_cut.update_layout(
    title=line_title,
    xaxis_title=line_x_label,
    yaxis_title=z_label,
    margin=dict(l=0, r=0, t=40, b=0),
    height=360,
)
st.plotly_chart(fig_cut, use_container_width=True)
