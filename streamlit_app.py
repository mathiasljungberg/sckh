"""Streamlit app for visualizing and comparing generic 2D surface files."""

from __future__ import annotations

from io import StringIO
from pathlib import Path
import sys
import tempfile

import numpy as np
import plotly.graph_objects as go
import streamlit as st

SRC_ROOT = Path(__file__).resolve().parent / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from python_scripts.dynamics_1d.constants import CONST
from python_scripts.dynamics_2d.io import read_surface_file_2d_raw


st.set_page_config(page_title="2D Surface Viewer", layout="wide")

VALID_SUFFIXES = {".dat", ".txt", ".pes"}
TRANSFORM_OPTIONS = ["raw", "abs", "sumsq"]


def loadtxt_fortran(filepath: Path) -> np.ndarray:
    """Load text file with support for Fortran D-format exponents."""
    text = filepath.read_text().replace("D", "E").replace("d", "e")
    return np.atleast_2d(np.loadtxt(StringIO(text)))


def convert_positions_to_si(
    x1_raw: np.ndarray,
    x2_raw: np.ndarray,
    position_units: str,
) -> tuple[np.ndarray, np.ndarray]:
    if position_units.lower() == "angstrom":
        return x1_raw * 1e-10, x2_raw * 1e-10
    if position_units.lower() == "bohr":
        return x1_raw * CONST.bohr, x2_raw * CONST.bohr
    raise ValueError(f"Unknown position units: {position_units}")


@st.cache_data(show_spinner=False)
def load_surface_from_path(
    path_str: str,
    position_units: str,
    index_order: str,
    mtime_ns: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str, list[str]]:
    """Load generic 2D surface data and convert coordinates to SI."""
    _ = mtime_ns  # cache key component to invalidate when file changes

    path = Path(path_str)
    data = loadtxt_fortran(path)
    if data.shape[1] < 3:
        raise ValueError(
            f"Surface file must have at least 3 columns (x1 x2 value), got {data.shape[1]}"
        )

    value_columns = list(range(2, data.shape[1]))
    x1_raw, x2_raw, values_raw, resolved_order = read_surface_file_2d_raw(
        path,
        value_columns=value_columns,
        index_order=index_order,
    )

    x1_si, x2_si = convert_positions_to_si(x1_raw, x2_raw, position_units)

    component_labels = [f"c{i}" for i in range(values_raw.shape[0])]
    return x1_si, x2_si, values_raw, resolved_order, component_labels


@st.cache_data(show_spinner=False)
def load_surface_from_bytes(
    file_bytes: bytes,
    position_units: str,
    index_order: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str, list[str]]:
    """Load generic 2D surface from uploaded bytes."""
    with tempfile.NamedTemporaryFile(mode="wb", suffix=".dat", delete=True) as tmp:
        tmp.write(file_bytes)
        tmp.flush()
        data = loadtxt_fortran(Path(tmp.name))
        if data.shape[1] < 3:
            raise ValueError(
                f"Surface file must have at least 3 columns (x1 x2 value), got {data.shape[1]}"
            )
        value_columns = list(range(2, data.shape[1]))
        x1_raw, x2_raw, values_raw, resolved_order = read_surface_file_2d_raw(
            Path(tmp.name),
            value_columns=value_columns,
            index_order=index_order,
        )
    x1_si, x2_si = convert_positions_to_si(x1_raw, x2_raw, position_units)
    component_labels = [f"c{i}" for i in range(values_raw.shape[0])]
    return x1_si, x2_si, values_raw, resolved_order, component_labels


def convert_positions_for_display(x: np.ndarray, units: str) -> np.ndarray:
    if units.lower() == "angstrom":
        return x * 1e10
    if units.lower() == "bohr":
        return x / CONST.bohr
    return x


def compute_display_field(
    values: np.ndarray,
    component_idx: int,
    transform: str,
) -> tuple[np.ndarray, str]:
    """Get scalar display field from component stack."""
    n_components = values.shape[0]
    component_idx = max(0, min(component_idx, n_components - 1))

    if transform == "raw":
        z = values[component_idx]
        label = f"c{component_idx}"
    elif transform == "abs":
        z = np.abs(values[component_idx])
        label = f"|c{component_idx}|"
    elif transform == "sumsq":
        z = np.sum(values**2, axis=0)
        label = "sum(c_i^2)"
    else:
        raise ValueError(f"Unknown transform: {transform}")

    return z, label


def cut_line(
    z: np.ndarray,
    x1: np.ndarray,
    x2: np.ndarray,
    cut_direction: str,
    fixed_value: float,
) -> tuple[np.ndarray, np.ndarray, str, float]:
    """Extract a 1D cut at the nearest coordinate to fixed_value."""
    if cut_direction == "x1 at fixed x2":
        idx = int(np.argmin(np.abs(x2 - fixed_value)))
        return x1, z[:, idx], "x1", float(x2[idx])

    idx = int(np.argmin(np.abs(x1 - fixed_value)))
    return x2, z[idx, :], "x2", float(x1[idx])


def ensure_state() -> None:
    if "active_entries" not in st.session_state:
        st.session_state.active_entries = []
    if "entry_counter" not in st.session_state:
        st.session_state.entry_counter = 0
    if "browser_dir" not in st.session_state:
        st.session_state.browser_dir = str(Path.cwd())
    if "selected_files_to_add" not in st.session_state:
        st.session_state.selected_files_to_add = []
    if "uploaded_files_store" not in st.session_state:
        st.session_state.uploaded_files_store = {}
    if "processed_upload_signatures" not in st.session_state:
        st.session_state.processed_upload_signatures = set()
    if "uploader_nonce" not in st.session_state:
        st.session_state.uploader_nonce = 0


def add_entry(path_str: str, component_idx: int = 0, transform: str = "raw") -> None:
    entry_id = f"entry_{st.session_state.entry_counter}"
    st.session_state.entry_counter += 1
    st.session_state.active_entries.append(
        {
            "id": entry_id,
            "source": "path",
            "path": path_str,
            "label": Path(path_str).name,
            "component_idx": component_idx,
            "transform": transform,
            "visible": True,
        }
    )


def add_uploaded_entry(filename: str, file_bytes: bytes) -> None:
    entry_id = f"entry_{st.session_state.entry_counter}"
    st.session_state.entry_counter += 1
    upload_key = f"{entry_id}:{filename}"
    st.session_state.uploaded_files_store[upload_key] = file_bytes
    st.session_state.active_entries.append(
        {
            "id": entry_id,
            "source": "upload",
            "upload_key": upload_key,
            "label": filename,
            "component_idx": 0,
            "transform": "raw",
            "visible": True,
        }
    )


def duplicate_entry(original: dict[str, object]) -> None:
    """Duplicate an existing active entry preserving its source."""
    source = str(original.get("source", "path"))
    entry_id = f"entry_{st.session_state.entry_counter}"
    st.session_state.entry_counter += 1

    duplicated = {
        "id": entry_id,
        "source": source,
        "label": str(original.get("label", "surface")),
        "component_idx": int(original.get("component_idx", 0)),
        "transform": str(original.get("transform", "raw")),
        "visible": True,
    }
    if source == "upload":
        duplicated["upload_key"] = str(original.get("upload_key", ""))
    else:
        duplicated["path"] = str(original.get("path", ""))

    st.session_state.active_entries.append(duplicated)


def browse_dir_listing(base_dir: Path) -> tuple[list[str], list[str]]:
    subdirs: list[str] = []
    files: list[str] = []

    if base_dir.parent != base_dir:
        subdirs.append("..")

    for item in sorted(base_dir.iterdir(), key=lambda p: p.name.lower()):
        if item.is_dir():
            subdirs.append(item.name)
        elif item.is_file() and item.suffix.lower() in VALID_SUFFIXES:
            files.append(item.name)

    return subdirs, files


def recursive_file_listing(base_dir: Path) -> list[str]:
    """List compatible files recursively relative to base_dir."""
    files: list[str] = []
    for item in sorted(base_dir.rglob("*"), key=lambda p: str(p).lower()):
        if item.is_file() and item.suffix.lower() in VALID_SUFFIXES:
            files.append(str(item.relative_to(base_dir)))
    return files


ensure_state()
st.title("2D Surface Viewer")

with st.sidebar:
    st.header("Files")

    st.subheader("File Picker")
    uploaded_files = st.file_uploader(
        "Choose files",
        type=["dat", "txt", "pes"],
        accept_multiple_files=True,
        help="Upload one or more surface files.",
        key=f"surface_uploader_{st.session_state.uploader_nonce}",
    )
    if not uploaded_files:
        st.session_state.processed_upload_signatures = set()
    else:
        added = 0
        for uploaded in uploaded_files:
            signature = f"{uploaded.name}:{uploaded.size}"
            if signature in st.session_state.processed_upload_signatures:
                continue
            add_uploaded_entry(uploaded.name, uploaded.getvalue())
            st.session_state.processed_upload_signatures.add(signature)
            added += 1
        if added > 0:
            # Reset uploader widget so selected files are not shown in sidebar anymore.
            st.session_state.uploader_nonce += 1
            st.rerun()

    with st.expander("Advanced local path tools", expanded=False):
        nav_cols = st.columns(2)
        if nav_cols[0].button("Up one level", use_container_width=True):
            st.session_state.browser_dir = str(
                Path(st.session_state.browser_dir).expanduser().resolve().parent
            )
            st.rerun()
        if nav_cols[1].button("Home (cwd)", use_container_width=True):
            st.session_state.browser_dir = str(Path.cwd())
            st.rerun()

        browser_dir_input = st.text_input(
            "Directory",
            value=st.session_state.browser_dir,
            key="browser_dir_input",
        )

        browser_dir = Path(browser_dir_input).expanduser().resolve()
        if browser_dir.is_dir():
            st.session_state.browser_dir = str(browser_dir)
            subdirs, files = browse_dir_listing(browser_dir)
            direct_subdirs = [d for d in subdirs if d != ".."]

            selected_subdir = st.selectbox(
                "Subfolder",
                options=["(stay in current directory)"] + direct_subdirs,
                index=0,
                key="selected_subdir",
            )
            if st.button("Open selected subfolder", key="open_folder", use_container_width=True):
                if selected_subdir != "(stay in current directory)":
                    st.session_state.browser_dir = str((browser_dir / selected_subdir).resolve())
                    st.rerun()
            st.caption(f"Current: `{browser_dir}`")

            recursive_search = st.checkbox("Search files recursively", value=False)
            file_filter = st.text_input(
                "Filter files",
                value="",
                placeholder="Type part of a filename",
            )

            if recursive_search:
                file_options = recursive_file_listing(browser_dir)
            else:
                file_options = files

            if file_filter.strip():
                needle = file_filter.strip().lower()
                file_options = [f for f in file_options if needle in f.lower()]

            selected_local_files = st.multiselect(
                "Local files",
                options=file_options,
                key="selected_files_to_add",
                help="Pick one or more files and click Add local files.",
            )

            if st.button("Add local files", use_container_width=True):
                for file_name in selected_local_files:
                    path_to_add = (browser_dir / file_name).resolve()
                    add_entry(str(path_to_add))
                st.session_state.selected_files_to_add = []
                st.rerun()
        else:
            st.warning("Directory does not exist.")

        manual_path = st.text_input(
            "Manual path",
            value="",
            placeholder="/path/to/surface.dat",
        )
        if st.button("Add manual path", use_container_width=True):
            manual = Path(manual_path).expanduser().resolve()
            if not manual.is_file():
                st.error("Manual path not found.")
            elif manual.suffix.lower() not in VALID_SUFFIXES:
                st.error("File suffix must be .dat, .txt, or .pes")
            else:
                add_entry(str(manual))
            st.rerun()

    with st.expander("Grid and Display options", expanded=False):
        position_units = st.selectbox("Position units in file", ["angstrom", "bohr"], index=0)
        index_order = st.selectbox(
            "Index ordering",
            ["auto", "C", "F"],
            index=0,
            help="Auto infers whether x2-fast (C) or x1-fast (F) ordering was used.",
        )
        show_range = st.checkbox("Zoom to value range", value=True)
        stride = st.select_slider("Downsample factor", options=[1, 2, 4, 6, 8, 10], value=2)
        contour_levels = st.select_slider("Contour levels", options=[10, 15, 20, 30, 40, 60], value=30)
        colormap = st.selectbox("Colormap", ["Viridis", "Cividis", "Plasma", "Magma", "Turbo"], index=0)

if not st.session_state.active_entries:
    st.info("Add at least one surface file from the sidebar.")
    st.stop()

# Load distinct paths once per rerun.
datasets: dict[str, dict[str, object]] = {}
for entry in st.session_state.active_entries:
    source = str(entry.get("source", "path"))
    if source == "upload":
        upload_key = str(entry.get("upload_key", ""))
        data_key = f"upload:{upload_key}"
    else:
        path_str = str(entry.get("path", ""))
        data_key = f"path:{path_str}"

    if data_key in datasets:
        continue

    try:
        if source == "upload":
            raw_bytes = st.session_state.uploaded_files_store.get(upload_key)
            if raw_bytes is None:
                datasets[data_key] = {"error": "Uploaded content is no longer available."}
                continue
            x1_si, x2_si, values, resolved_order, component_labels = load_surface_from_bytes(
                raw_bytes,
                position_units=position_units,
                index_order=index_order,
            )
            datasets[data_key] = {
                "x1": convert_positions_for_display(x1_si, position_units),
                "x2": convert_positions_for_display(x2_si, position_units),
                "values": values,
                "resolved_order": resolved_order,
                "component_labels": component_labels,
            }
            continue

        path = Path(path_str)
        if not path.is_file():
            datasets[data_key] = {"error": "File not found."}
            continue

        mtime_ns = path.stat().st_mtime_ns
        x1_si, x2_si, values, resolved_order, component_labels = load_surface_from_path(
            path_str,
            position_units=position_units,
            index_order=index_order,
            mtime_ns=mtime_ns,
        )
        datasets[data_key] = {
            "x1": convert_positions_for_display(x1_si, position_units),
            "x2": convert_positions_for_display(x2_si, position_units),
            "values": values,
            "resolved_order": resolved_order,
            "component_labels": component_labels,
            "mtime_ns": mtime_ns,
        }
    except Exception as exc:  # noqa: BLE001
        datasets[data_key] = {"error": str(exc)}

valid_component_counts = [
    int(d["values"].shape[0])  # type: ignore[index]
    for d in datasets.values()
    if "error" not in d
]
max_component_count = max(valid_component_counts) if valid_component_counts else 1

with st.sidebar:
    st.header("Global Override")
    use_global_override = st.checkbox(
        "Use global component/transform",
        value=False,
        help="Overrides per-file selectors for plotting.",
    )
    global_transform = st.selectbox(
        "Global transform",
        options=TRANSFORM_OPTIONS,
        index=0,
        format_func=lambda t: {
            "raw": "Raw",
            "abs": "Absolute value",
            "sumsq": "Sum of squares",
        }[t],
    )
    if global_transform in {"raw", "abs"}:
        global_component_idx = st.selectbox(
            "Global component",
            options=list(range(max_component_count)),
            format_func=lambda i: f"c{i}",
            index=0,
        )
    else:
        global_component_idx = 0

st.subheader("Active Files")
if use_global_override:
    st.caption("Global override is active: duplicate entries are collapsed to one plot per file.")
remove_id: str | None = None
duplicate_id: str | None = None

for idx, entry in enumerate(st.session_state.active_entries, start=1):
    entry_id = entry["id"]
    source = str(entry.get("source", "path"))
    path_str = str(entry.get("path", ""))
    upload_key = str(entry.get("upload_key", ""))
    data_key = f"upload:{upload_key}" if source == "upload" else f"path:{path_str}"
    path_name = str(entry.get("label", Path(path_str).name))
    dataset = datasets.get(data_key, {"error": "Missing dataset"})

    with st.container(border=True):
        head_cols = st.columns([6, 2, 2])
        head_cols[0].markdown(f"**{idx}. {path_name}**")
        if head_cols[1].button("Duplicate", key=f"dup_{entry_id}", use_container_width=True):
            duplicate_id = entry_id
        if head_cols[2].button("Remove", key=f"rm_{entry_id}", use_container_width=True):
            remove_id = entry_id

        if source == "upload":
            st.caption("Source: uploaded content")
        else:
            st.caption(path_str)

        if "error" in dataset:
            st.error(f"Load error: {dataset['error']}")
            entry["visible"] = False
            continue

        comp_labels: list[str] = dataset["component_labels"]  # type: ignore[assignment]

        settings_cols = st.columns(4)
        entry["visible"] = settings_cols[0].checkbox(
            "Active",
            value=entry.get("visible", True),
            key=f"visible_{entry_id}",
        )

        max_idx = len(comp_labels) - 1
        current_component = int(entry.get("component_idx", 0))
        current_component = max(0, min(current_component, max_idx))

        selected_component = settings_cols[1].selectbox(
            "Component",
            options=list(range(len(comp_labels))),
            format_func=lambda i: comp_labels[i],
            index=current_component,
            key=f"component_{entry_id}",
        )
        entry["component_idx"] = int(selected_component)

        current_transform = entry.get("transform", "raw")
        if current_transform not in TRANSFORM_OPTIONS:
            current_transform = "raw"

        entry["transform"] = settings_cols[2].selectbox(
            "Transform",
            options=TRANSFORM_OPTIONS,
            index=TRANSFORM_OPTIONS.index(current_transform),
            format_func=lambda t: {
                "raw": "Raw",
                "abs": "Absolute value",
                "sumsq": "Sum of squares",
            }[t],
            key=f"transform_{entry_id}",
        )

        settings_cols[3].write("Detected order")
        settings_cols[3].code(str(dataset["resolved_order"]))

if remove_id is not None:
    st.session_state.active_entries = [
        e for e in st.session_state.active_entries if e["id"] != remove_id
    ]
    st.rerun()

if duplicate_id is not None:
    original = next(
        (e for e in st.session_state.active_entries if e["id"] == duplicate_id),
        None,
    )
    if original is not None:
        duplicate_entry(original)
    st.rerun()

visible_fields: list[dict[str, object]] = []
seen_sources: set[str] = set()
for entry in st.session_state.active_entries:
    if not entry.get("visible", True):
        continue

    source = str(entry.get("source", "path"))
    path_str = str(entry.get("path", ""))
    upload_key = str(entry.get("upload_key", ""))
    data_key = f"upload:{upload_key}" if source == "upload" else f"path:{path_str}"
    source_key = data_key
    if use_global_override and source_key in seen_sources:
        continue
    dataset = datasets.get(data_key)
    if not dataset or "error" in dataset:
        continue

    values: np.ndarray = dataset["values"]  # type: ignore[assignment]
    if use_global_override:
        component_idx = int(global_component_idx)
        transform = str(global_transform)
    else:
        component_idx = int(entry.get("component_idx", 0))
        transform = str(entry.get("transform", "raw"))
    z, z_label = compute_display_field(values, component_idx, transform)

    filename = str(entry.get("label", Path(path_str).name))
    if transform == "sumsq":
        name = f"{filename} | sum(c_i^2)"
    elif transform == "abs":
        name = f"{filename} | |c{component_idx}|"
    else:
        name = f"{filename} | c{component_idx}"

    visible_fields.append(
        {
            "name": name,
            "z_label": z_label,
            "x1": dataset["x1"],
            "x2": dataset["x2"],
            "z": z,
            "path": path_str if source == "path" else f"upload:{filename}",
        }
    )
    if use_global_override:
        seen_sources.add(source_key)

if not visible_fields:
    st.warning("No active, valid entries to plot.")
    st.stop()

x1_label = f"x1 ({position_units})"
x2_label = f"x2 ({position_units})"

z_all_min = min(float(np.min(field["z"])) for field in visible_fields)
z_all_max = max(float(np.max(field["z"])) for field in visible_fields)

info_cols = st.columns(4)
info_cols[0].metric("Active entries", f"{len(visible_fields)}")
info_cols[1].metric("Distinct files", f"{len({f['path'] for f in visible_fields})}")
info_cols[2].metric("Global min", f"{z_all_min:.6g}")
info_cols[3].metric("Global max", f"{z_all_max:.6g}")

if show_range:
    z_low, z_high = st.slider(
        "Value range",
        min_value=z_all_min,
        max_value=z_all_max,
        value=(z_all_min, z_all_max),
        format="%.6g",
    )
    z_range = (z_low, z_high)
else:
    z_range = (z_all_min, z_all_max)

contour_size = (z_range[1] - z_range[0]) / contour_levels
if contour_size <= 0:
    contour_size = 1.0

col1, col2 = st.columns(2)

with col1:
    st.subheader("3D Surface")
    fig_surface = go.Figure()
    for idx, field in enumerate(visible_fields):
        x1 = field["x1"]
        x2 = field["x2"]
        z = field["z"]
        x1_plot = x1[::stride]
        x2_plot = x2[::stride]
        z_plot = z[::stride, ::stride]

        fig_surface.add_trace(
            go.Surface(
                x=x1_plot,
                y=x2_plot,
                z=z_plot.T,
                colorscale=colormap,
                cmin=z_range[0],
                cmax=z_range[1],
                opacity=0.62,
                showscale=idx == 0,
                name=str(field["name"]),
            )
        )

    fig_surface.update_layout(
        scene=dict(
            xaxis_title=x1_label,
            yaxis_title=x2_label,
            zaxis_title="Value",
            zaxis=dict(range=[z_range[0], z_range[1]]),
        ),
        margin=dict(l=0, r=0, t=30, b=0),
        height=620,
    )
    st.plotly_chart(fig_surface, use_container_width=True)

with col2:
    st.subheader("Contour")
    fig_contour = go.Figure()
    for field in visible_fields:
        x1 = field["x1"]
        x2 = field["x2"]
        z = field["z"]
        x1_plot = x1[::stride]
        x2_plot = x2[::stride]
        z_plot = z[::stride, ::stride]

        fig_contour.add_trace(
            go.Contour(
                x=x1_plot,
                y=x2_plot,
                z=z_plot.T,
                colorscale=colormap,
                zmin=z_range[0],
                zmax=z_range[1],
                contours=dict(
                    showlabels=False,
                    start=z_range[0],
                    end=z_range[1],
                    size=contour_size,
                ),
                name=str(field["name"]),
                opacity=0.85,
            )
        )

    fig_contour.update_layout(
        xaxis_title=x1_label,
        yaxis_title=x2_label,
        margin=dict(l=0, r=0, t=30, b=0),
        height=620,
    )
    st.plotly_chart(fig_contour, use_container_width=True)

st.subheader("Surface Cuts")
reference = visible_fields[0]
ref_x1: np.ndarray = reference["x1"]  # type: ignore[assignment]
ref_x2: np.ndarray = reference["x2"]  # type: ignore[assignment]

cut_cols = st.columns(3)
with cut_cols[0]:
    cut_direction = st.radio("Cut direction", ["x1 at fixed x2", "x2 at fixed x1"], index=0)

if cut_direction == "x1 at fixed x2":
    with cut_cols[1]:
        ref_idx = st.slider("x2 index (reference)", 0, len(ref_x2) - 1, len(ref_x2) // 2)
    fixed_value = float(ref_x2[ref_idx])
    fixed_label = f"x2 ~ {fixed_value:.6g} {position_units}"
else:
    with cut_cols[1]:
        ref_idx = st.slider("x1 index (reference)", 0, len(ref_x1) - 1, len(ref_x1) // 2)
    fixed_value = float(ref_x1[ref_idx])
    fixed_label = f"x1 ~ {fixed_value:.6g} {position_units}"

with cut_cols[2]:
    st.write("Fixed coordinate")
    st.write(fixed_label)

fig_cut = go.Figure()
for field in visible_fields:
    x1 = field["x1"]
    x2 = field["x2"]
    z = field["z"]

    line_x, line_y, cut_axis, cut_fixed = cut_line(
        z=z,
        x1=x1,
        x2=x2,
        cut_direction=cut_direction,
        fixed_value=fixed_value,
    )

    fig_cut.add_trace(
        go.Scatter(
            x=line_x,
            y=line_y,
            mode="lines",
            line=dict(width=2),
            name=str(field["name"]),
            hovertemplate=(
                "%{x:.6g}, %{y:.6g}<br>"
                f"{cut_axis}-cut at {cut_fixed:.6g} {position_units}"
                "<extra></extra>"
            ),
        )
    )

line_x_label = x1_label if cut_direction == "x1 at fixed x2" else x2_label
fig_cut.update_layout(
    title=f"1D Cuts at {fixed_label}",
    xaxis_title=line_x_label,
    yaxis_title="Value",
    margin=dict(l=0, r=0, t=40, b=0),
    height=380,
)
st.plotly_chart(fig_cut, use_container_width=True)
