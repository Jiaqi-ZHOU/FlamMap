from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib
from matplotlib import font_manager

matplotlib.use("Agg")


def _has_font(font_name: str) -> bool:
    try:
        font_path = font_manager.findfont(font_name, fallback_to_default=False)
    except ValueError:
        return False
    return bool(font_path)


if _has_font("Arial"):
    matplotlib.rcParams["font.family"] = "Arial"
    matplotlib.rcParams["font.sans-serif"] = ["Arial"]
    matplotlib.rcParams["mathtext.fontset"] = "custom"
    matplotlib.rcParams["mathtext.rm"] = "Arial"
    matplotlib.rcParams["mathtext.it"] = "Arial:italic"
    matplotlib.rcParams["mathtext.bf"] = "Arial:bold"

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
import numpy as np


TRIANGLE_HEIGHT = math.sqrt(3.0) / 2.0
AIR_BARY = np.array([0.21, 0.79, 0.0], dtype=float)
FUEL_BARY = np.array([0.0, 0.0, 1.0], dtype=float)


def bary_to_cart(o2: float, n2: float, fuel: float) -> tuple[float, float]:
    x = n2 + 0.5 * fuel
    y = TRIANGLE_HEIGHT * fuel
    return x, y


def cart_to_bary(x: float, y: float) -> np.ndarray:
    fuel = y / TRIANGLE_HEIGHT
    n2 = x - 0.5 * fuel
    o2 = 1.0 - n2 - fuel
    return np.array([o2, n2, fuel], dtype=float)


def load_dat_file(
    dat_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(dat_file, comments="#")
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]


def _formula_to_mathtext(formula: str) -> str:
    parts = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
    if not parts:
        return formula
    tokens = []
    for element, count in parts:
        if count:
            tokens.append(rf"\mathrm{{{element}}}_{{{count}}}")
        else:
            tokens.append(rf"\mathrm{{{element}}}")
    return f"${''.join(tokens)}$"


def _formula_counts(formula: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[element] = counts.get(element, 0) + (int(count) if count else 1)
    return counts


def _draw_grid(ax, tick_values: np.ndarray) -> None:
    for value in tick_values[1:-1]:
        p1 = bary_to_cart(1.0 - value, 0.0, value)
        p2 = bary_to_cart(0.0, 1.0 - value, value)
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            color="#7a7a7a",
            linewidth=0.6,
            alpha=0.9,
            zorder=2,
        )

        p1 = bary_to_cart(value, 1.0 - value, 0.0)
        p2 = bary_to_cart(value, 0.0, 1.0 - value)
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            color="#7a7a7a",
            linewidth=0.6,
            alpha=0.9,
            zorder=2,
        )

        p1 = bary_to_cart(1.0 - value, value, 0.0)
        p2 = bary_to_cart(0.0, value, 1.0 - value)
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            color="#7a7a7a",
            linewidth=0.6,
            alpha=0.9,
            zorder=2,
        )


def _draw_axis_ticks(ax, tick_values: np.ndarray) -> None:
    tick_fontsize = 7
    tick_len = 0.018

    left_tick = np.array([-1.0, 0.0])
    for value in tick_values:
        fuel = value
        x, y = bary_to_cart(1.0 - fuel, 0.0, fuel)
        end = np.array([x, y]) + tick_len * left_tick
        label_pos = np.array([x, y]) + 0.042 * left_tick
        ax.plot([x, end[0]], [y, end[1]], color="black", linewidth=0.6)
        ax.text(
            label_pos[0],
            label_pos[1],
            f"{value:.1f}",
            ha="right",
            va="center",
            fontsize=tick_fontsize,
        )

    right_angle = 60
    right_tick = np.array(
        [math.cos(math.radians(right_angle)), math.sin(math.radians(right_angle))]
    )
    for value in tick_values:
        fuel = value
        x, y = bary_to_cart(0.0, 1.0 - fuel, fuel)
        end = np.array([x, y]) + tick_len * right_tick
        label_pos = np.array([x, y]) + 0.05 * right_tick
        ax.plot([x, end[0]], [y, end[1]], color="black", linewidth=0.6)
        ax.text(
            label_pos[0],
            label_pos[1],
            f"{1.0 - value:.1f}",
            ha="left",
            va="center",
            rotation=right_angle,
            rotation_mode="anchor",
            fontsize=tick_fontsize,
        )

    bottom_tick = np.array([0.5, -math.sqrt(3.0) / 2.0])
    for value in tick_values:
        n2 = value
        x, y = bary_to_cart(1.0 - n2, n2, 0.0)
        end = np.array([x, y]) + tick_len * bottom_tick
        label_pos = np.array([x, y]) + 0.07 * bottom_tick
        ax.plot([x, end[0]], [y, end[1]], color="black", linewidth=0.6)
        ax.text(
            label_pos[0],
            label_pos[1],
            f"{1.0 - value:.1f}",
            ha="right",
            va="top",
            rotation=-60,
            rotation_mode="anchor",
            fontsize=tick_fontsize,
        )


def _draw_reference_lines(ax, formula: str) -> None:
    air_point = np.array(bary_to_cart(*AIR_BARY))
    fuel_point = np.array(bary_to_cart(*FUEL_BARY))
    ax.plot(
        [fuel_point[0], air_point[0]],
        [fuel_point[1], air_point[1]],
        color="red",
        linewidth=1.2,
        zorder=4,
    )

    air_vector = air_point - fuel_point
    air_angle = math.degrees(math.atan2(air_vector[1], air_vector[0]))
    air_label = 0.52 * fuel_point + 0.48 * air_point + np.array([0.025, 0.02])
    ax.text(
        air_label[0],
        air_label[1],
        "Air line",
        color="red",
        fontsize=9,
        rotation=air_angle,
        rotation_mode="anchor",
        ha="left",
        va="center",
    )

    counts = _formula_counts(formula)
    stoich_o2 = (
        counts.get("C", 0) + 0.25 * counts.get("H", 0) - 0.5 * counts.get("O", 0)
    )
    if stoich_o2 <= 0:
        return

    fuel_fraction = 1.0 / (stoich_o2 + 1.0)
    x0 = np.array(bary_to_cart(0.0, 1.0, 0.0))
    x1 = np.array(bary_to_cart(1.0 - fuel_fraction, 0.0, fuel_fraction))
    ax.plot([x0[0], x1[0]], [x0[1], x1[1]], color="blue", linewidth=1.0, zorder=4)

    stoich_vector = x1 - x0
    stoich_angle = math.degrees(math.atan2(stoich_vector[1], stoich_vector[0]))
    if stoich_angle > 90.0 or stoich_angle < -90.0:
        stoich_angle += 180.0
    label_point = 0.5 * x0 + 0.5 * x1 + np.array([-0.03, 0.04])
    ax.text(
        label_point[0],
        label_point[1],
        "Stoichiometric line",
        color="blue",
        fontsize=9,
        rotation=stoich_angle,
        rotation_mode="anchor",
        ha="right",
        va="center",
    )


def _line_intersection(p1, p2, q1, q2, eps: float = 1e-9):
    p = np.asarray(p1, dtype=float)
    r = np.asarray(p2, dtype=float) - p
    q = np.asarray(q1, dtype=float)
    s = np.asarray(q2, dtype=float) - q

    rxs = r[0] * s[1] - r[1] * s[0]
    q_p = q - p
    qpxr = q_p[0] * r[1] - q_p[1] * r[0]

    if abs(rxs) < eps:
        return None

    t = (q_p[0] * s[1] - q_p[1] * s[0]) / rxs
    u = qpxr / rxs
    if -eps <= t <= 1.0 + eps and -eps <= u <= 1.0 + eps:
        return p + t * r
    return None


def _deduplicate_points(
    points: list[np.ndarray], tol: float = 1e-5
) -> list[np.ndarray]:
    unique: list[np.ndarray] = []
    for point in points:
        if not any(np.linalg.norm(point - existing) < tol for existing in unique):
            unique.append(point)
    return unique


def extract_contour_segments(dat_file: str | Path, threshold_temperature: float):
    o2, n2, fuel, temperatures = load_dat_file(dat_file)
    x, y = bary_to_cart(o2, n2, fuel)
    triangulation = mtri.Triangulation(x, y)

    fig, ax = plt.subplots()
    contour = ax.tricontour(triangulation, temperatures, levels=[threshold_temperature])
    segments = [
        np.asarray(seg, dtype=float) for seg in contour.allsegs[0] if len(seg) >= 2
    ]
    plt.close(fig)
    return segments


def compute_flammability_limits(
    dat_file: str | Path, threshold_temperature: float
) -> tuple[float, float, list[np.ndarray]]:
    segments = extract_contour_segments(dat_file, threshold_temperature)
    if not segments:
        return math.nan, math.nan, segments

    air_point = np.array(bary_to_cart(*AIR_BARY), dtype=float)
    fuel_point = np.array(bary_to_cart(*FUEL_BARY), dtype=float)

    intersections: list[np.ndarray] = []
    for segment in segments:
        for idx in range(len(segment) - 1):
            intersection = _line_intersection(
                segment[idx], segment[idx + 1], air_point, fuel_point
            )
            if intersection is not None:
                intersections.append(intersection)

    intersections = _deduplicate_points(intersections)
    if not intersections:
        return math.nan, math.nan, segments

    bary_points = [cart_to_bary(point[0], point[1]) for point in intersections]
    fuel_percents = sorted(float(point[2] * 100.0) for point in bary_points)
    lfl = fuel_percents[0]
    ufl = 100.0 if len(fuel_percents) == 1 else fuel_percents[-1]
    return lfl, ufl, segments


def plot_phase_diagram(
    dat_file: str | Path,
    pdf_file: str | Path,
    *,
    formula: str,
    threshold_temperature: float,
    lfl_percent: float,
    ufl_percent: float,
) -> None:
    o2, n2, fuel, temperatures = load_dat_file(dat_file)
    x, y = bary_to_cart(o2, n2, fuel)
    triangulation = mtri.Triangulation(x, y)
    segments = extract_contour_segments(dat_file, threshold_temperature)
    color_norm = mcolors.Normalize(vmin=0.0, vmax=3600.0, clip=True)
    cmap = plt.get_cmap("jet")

    fig, ax = plt.subplots(figsize=(6.0, 6.2), dpi=300)
    triangle_patch = Polygon(
        [
            bary_to_cart(1.0, 0.0, 0.0),
            bary_to_cart(0.0, 1.0, 0.0),
            bary_to_cart(0.0, 0.0, 1.0),
        ],
        closed=True,
        facecolor=cmap(color_norm(3600.0)),
        edgecolor="none",
        zorder=0,
    )
    ax.add_patch(triangle_patch)
    filled = ax.tricontourf(
        triangulation,
        temperatures,
        levels=np.linspace(0, 3600, 25),
        cmap=cmap,
        norm=color_norm,
    )
    tick_values = np.round(np.linspace(0.0, 1.0, 11), 1)
    _draw_grid(ax, tick_values)
    ax.tricontour(
        triangulation,
        temperatures,
        levels=[threshold_temperature],
        colors="black",
        linewidths=1.0,
    )

    vertices = np.array(
        [
            bary_to_cart(1.0, 0.0, 0.0),
            bary_to_cart(0.0, 1.0, 0.0),
            bary_to_cart(0.0, 0.0, 1.0),
            bary_to_cart(1.0, 0.0, 0.0),
        ]
    )
    ax.plot(vertices[:, 0], vertices[:, 1], color="black", linewidth=0.8)

    _draw_reference_lines(ax, formula)

    for seg in segments:
        ax.plot(seg[:, 0], seg[:, 1], color="black", linewidth=0.8)

    for value in [lfl_percent, ufl_percent]:
        if math.isnan(value):
            continue
        fuel_frac = value / 100.0
        point = AIR_BARY + fuel_frac * (FUEL_BARY - AIR_BARY)
        px, py = bary_to_cart(*point)
        ax.scatter(px, py, color="purple", s=18, marker="x", zorder=5)

    _draw_axis_ticks(ax, tick_values)

    ax.text(
        *bary_to_cart(0.45, 0.55, -0.10),
        r"$\mathrm{O}_{2}$",
        ha="center",
        va="top",
        fontsize=10,
    )
    ax.text(
        *bary_to_cart(-0.15, 0.575, 0.575),
        r"$\mathrm{N}_{2}$",
        ha="left",
        va="center",
        fontsize=10,
        rotation=-60,
        rotation_mode="anchor",
    )
    ax.text(
        *bary_to_cart(0.575, -0.15, 0.575),
        _formula_to_mathtext(formula),
        ha="right",
        va="center",
        fontsize=11,
        rotation=60,
        rotation_mode="anchor",
    )

    ax.set_aspect("equal")
    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.12, TRIANGLE_HEIGHT + 0.16)
    ax.axis("off")

    fig.subplots_adjust(left=0.08, right=0.92, top=0.96, bottom=0.18)
    ax_pos = ax.get_position()
    cbar_width = ax_pos.width * 0.72
    cbar_x0 = ax_pos.x0 + 0.5 * (ax_pos.width - cbar_width)
    cax = fig.add_axes([cbar_x0, ax_pos.y0 - 0.06, cbar_width, 0.025])
    cbar = fig.colorbar(filled, cax=cax, orientation="horizontal")
    cbar.set_label("Calculated adiabatic flame temperature (K)")
    cbar.set_ticks([0, 1200, 2400, 3600])
    cbar.ax.tick_params(labelsize=8)

    pdf_file = Path(pdf_file)
    pdf_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_file, bbox_inches="tight")
    plt.close(fig)
