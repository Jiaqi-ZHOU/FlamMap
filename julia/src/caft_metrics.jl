module CAFTMetrics

# =========================
# Imports
# =========================
using LinearAlgebra
using Printf
using DelimitedFiles
using Makie
using CairoMakie
using TernaryDiagrams

using CSV, DataFrames
using OrderedCollections

# =========================
# Constants + small utils
# =========================
const airline_oxygen = 0.21
const airline_nitrogen = 0.79

round3(x) = round(Float64(x); digits=3)
labelfontsize = 9.333
tickfontsize = 8 - 1
export write_results_csv, write_results_csv_full, plot_PD, compute_PD_metrics

# =========================
# CSV helpers
# =========================
function find_Hf_298K(jsonPrefix, Hf_csv::AbstractString)
    println("Finding Hf_298K for jsonPrefix: $(jsonPrefix) in $(Hf_csv)...")
    df = CSV.read(Hf_csv, DataFrame)
    row = filter(row -> row.jsonPrefix == jsonPrefix, df)
    Hf = row[1, :Hf_298K_kJ]
    return Printf.@sprintf("%.1f", Hf)
end

"""
Create a new CSV that contains jsonPrefix/formula + LFL/UFL at multiple `levels`.
- Writes LFL/UFL in **percentage** (×100), rounded to 3 decimals.
"""
function write_results_csv(dat_dir::AbstractString, out_csv::AbstractString; levels)
    files = sort(filter(f -> endswith(f, ".dat"), readdir(dat_dir)))
    total = length(files)
    total == 0 && error("No .dat files found in $dat_dir")

    # output table schema
    df = DataFrame(jsonPrefix=String[], formula=String[])
    for level in levels
        df[!, Symbol("LFL_$(level)")] = Float64[]
        df[!, Symbol("UFL_$(level)")] = Float64[]
    end

    # collect into dict (avoid duplicates)
    rows = Dict{String,Dict{Symbol,Any}}()

    for (i, filename) in enumerate(files)
        fullpath = joinpath(dat_dir, filename)
        println("[$i/$total] Processing: $filename")

        for level in levels
            try
                res = compute_PD_metrics(fullpath, level)
                jp = String(res.jsonPrefix)

                if !haskey(rows, jp)
                    rows[jp] = Dict{Symbol,Any}(
                        :jsonPrefix => jp,
                        :formula => String(res.formula),
                    )
                end

                # store as percentage
                rows[jp][Symbol("LFL_$(level)")] = round3(res.lfl[3] * 100)
                rows[jp][Symbol("UFL_$(level)")] = round3(res.ufl[3] * 100)

            catch err
                @warn "Failed $filename at $(level)K" exception = (err, catch_backtrace())
            end
        end
    end

    # assemble dataframe
    for (_, d) in rows
        row = Dict{Symbol,Any}(
            :jsonPrefix => d[:jsonPrefix],
            :formula => d[:formula],
        )
        for level in levels
            row[Symbol("LFL_$(level)")] = get(d, Symbol("LFL_$(level)"), missing)
            row[Symbol("UFL_$(level)")] = get(d, Symbol("UFL_$(level)"), missing)
        end
        push!(df, row)
    end

    sort!(df, :formula)
    CSV.write(out_csv, df)
    println("Wrote: $(abspath(out_csv))")
    return df
end


function write_results_csv_full(dat_dir::AbstractString, out_csv::AbstractString, Hf_csv::AbstractString;
    levels=(1600,))
    # Need to list the dir to know how many files we have for progress tracking
    # for the number of csvfile row. 
    files = sort(filter(f -> endswith(f, ".dat"), readdir(dat_dir)))
    total = length(files)
    total == 0 && error("No .dat files found in $dat_dir")

    # output schema
    df = DataFrame(jsonPrefix=String[], formula=String[], Hf_298K_kJ=Union{Missing,Float64}[])
    for level in levels
        df[!, Symbol("LFL_$(level)")] = Float64[]
        df[!, Symbol("UFL_$(level)")] = Float64[]
    end

    rows = Dict{String,Dict{Symbol,Any}}()

    for (i, filename) in enumerate(files)
        fullpath = joinpath(dat_dir, filename)
        println("[$i/$total] Processing: $filename")

        for level in levels
            try
                res = compute_PD_metrics(fullpath, level)
                jp = String(res.jsonPrefix)

                if !haskey(rows, jp)
                    rows[jp] = Dict{Symbol,Any}(
                        :jsonPrefix => jp,
                        :formula => String(res.formula),
                        :Hf_298K_kJ => missing,  # fill below
                    )
                    # Hf is independent of level → only lookup once
                    try
                        hf_str = find_Hf_298K(jp, Hf_csv)      # returns like "-74.5"
                        rows[jp][:Hf_298K_kJ] = parse(Float64, hf_str)
                    catch
                        rows[jp][:Hf_298K_kJ] = missing
                    end
                end

                rows[jp][Symbol("LFL_$(level)")] = round3(res.lfl[3] * 100)
                rows[jp][Symbol("UFL_$(level)")] = round3(res.ufl[3] * 100)

            catch err
                @warn "Failed $filename at $(level)K" exception = (err, catch_backtrace())
            end
        end
    end

    # assemble dataframe
    for (_, d) in rows
        row = Dict{Symbol,Any}(
            :jsonPrefix => d[:jsonPrefix],
            :formula => d[:formula],
            :Hf_298K_kJ => d[:Hf_298K_kJ],
        )
        for level in levels
            row[Symbol("LFL_$(level)")] = get(d, Symbol("LFL_$(level)"), missing)
            row[Symbol("UFL_$(level)")] = get(d, Symbol("UFL_$(level)"), missing)
        end
        push!(df, row)
    end

    sort!(df, :formula)
    CSV.write(out_csv, df)
    println("Wrote: $(abspath(out_csv))")
    return df
end

# =========================
# Core computation (no plotting)
# =========================

"""
Read PD .dat file:
columns: a1 a2 a3 mus
"""
function load_PD_data(filename::AbstractString)
    data = readdlm(filename)
    a1 = data[:, 1]
    a2 = data[:, 2]
    a3 = data[:, 3]
    mus = data[:, 4]
    return a1, a2, a3, mus
end

"""
Compute contour lines in Cartesian coordinates for a given temperature `level`.
Returns Vector{Vector{Point2}}.
"""
function get_contour(as, bs, cs, ds, level)
    xs = Float64[]
    ys = Float64[]
    ws = Float64[]
    for (a, b, c, d) in zip(as, bs, cs, ds)
        carts = TernaryDiagrams.R * [a, b, c]
        push!(xs, carts[2])
        push!(ys, carts[3])
        push!(ws, d)
    end

    scaled_coords = TernaryDiagrams.delaunay_scale.([
        TernaryDiagrams.gp.Point2D.(x, y) for (x, y) in zip(xs, ys)
    ])

    level_edges, _ = TernaryDiagrams.contour_triangle(scaled_coords, [level], ws, 1)
    edges = TernaryDiagrams.split_edges(level_edges[1])

    lines = map(edges) do curve
        [Point2(TernaryDiagrams.unpack(TernaryDiagrams.delaunay_unscale(v))...) for v in curve]
    end
    return lines
end


"""
Compute LFL/UFL intersections with airline.
Returns (lfl, ufl) in barycentric coordinates.
"""
function get_fl(contour; thr_dist=5e-2)
    p_air = TernaryDiagrams.from_bary_to_cart(airline_oxygen, airline_nitrogen, 1 - airline_oxygen - airline_nitrogen)
    p0 = TernaryDiagrams.from_bary_to_cart(0, 0, 1)

    yair(x) = p_air[2] + (x - p_air[1]) * (p_air[2] - p0[2]) / (p_air[1] - p0[1])

    xs = Float64[]
    ys = Float64[]
    for line in contour
        for p in line
            push!(xs, p[1])
            push!(ys, p[2])
        end
    end

    ys_air = yair.(xs)
    diff = abs.(ys_air .- ys)
    idx = sortperm(diff)

    i1 = idx[1]
    e1 = [xs[i1], ys[i1]]
    i2 = idx[2]
    e2 = [xs[i2], ys[i2]]

    dist = norm(e1 - e2)
    j = 2
    while (dist < thr_dist) && (j < length(xs))
        j += 1
        i2 = idx[j]
        e2 = [xs[i2], ys[i2]]
        dist = norm(e1 - e2)
    end

    # ensure e1 is lower limit
    if e1[2] > e2[2]
        e1, e2 = e2, e1
    end

    lfl = TernaryDiagrams.from_cart_to_bary(e1...)
    ufl = TernaryDiagrams.from_cart_to_bary(e2...)
    return lfl, ufl
end

"""
Compute metrics from `.dat` for a given threshold temperature `level`.
"""
function compute_PD_metrics(dat_filename::AbstractString, level::Real)
    lines = readlines(dat_filename)
    lines = filter(!startswith("#"), lines)
    data = [parse.(Float64, split(line)) for line in lines]
    data = reduce(vcat, permutedims.(data))

    a1 = data[:, 1]
    a2 = data[:, 2]
    a3 = data[:, 3]
    mus = data[:, 4]
    max_mus = maximum(filter(!isnan, mus))

    base = splitext(basename(dat_filename))[1]
    formula, jsonPrefix = split(base, "_")[1:2]

    if max_mus < level
        @warn("No data points above the level $level in $(dat_filename). Skipping.")
        lfl = [0, 0, 1]
        ufl = [1, 0, 0]
        t = nothing
    else
        t = get_contour(a1, a2, a3, mus, level)

        if any(occursin(x, dat_filename) for x in ["C5H8", "C5H10", "C5H12", "C4H10O", "C4H10",
            "C4H11N", "C4H8_6cec", "C4H9N", "CH4O3_6062"])
            thr_dist = 1e-2
        else
            thr_dist = 5e-2
        end

        lfl, ufl = get_fl(t; thr_dist)

        if ufl[3] - lfl[3] < 1e-3
            lfl, ufl = get_fl(t; thr_dist=1e-3)
            ufl = [0, 0, 1]
        elseif ufl[1] < 1e-3
            ufl = [0, 0, 1]
        end
    end

    return (; base, formula, jsonPrefix, a1, a2, a3, mus, max_mus, t, lfl, ufl)
end

# =========================
# Plotting helpers
# =========================

function count_elements(formula::AbstractString)
    element_counts = OrderedDict{String,Int}()
    regex = r"([A-Z][a-z]?)(\d*)"
    for match in eachmatch(regex, formula)
        element = match[1]
        count = isempty(match[2]) ? 1 : parse(Int, match[2])
        element_counts[element] = get(element_counts, element, 0) + count
    end
    return element_counts
end

function rich_formula(formula::AbstractString)
    unicode_subscripts = Dict(
        '0' => "₀", '1' => "₁", '2' => "₂", '3' => "₃", '4' => "₄",
        '5' => "₅", '6' => "₆", '7' => "₇", '8' => "₈", '9' => "₉",
    )
    elem_counts = count_elements(formula)
    form = ""
    for (e, n) in elem_counts
        if n == 1
            form *= e
        else
            sub = join(unicode_subscripts[c] for c in string(n))
            form *= e * sub
        end
    end
    return form
end

function rich_formula_rt(formula::AbstractString)
    elem_counts = count_elements(formula)
    parts = Any[]
    for (e, n) in elem_counts
        push!(parts, e)
        if n != 1
            push!(parts, subscript(string(n)))
        end
    end
    return rich(parts...)
end


function plot_airline(ax)
    ternarylines!(ax,
        [0.0, airline_oxygen],
        [0.0, airline_nitrogen],
        [1.0, 0.0];
        color=:red, linewidth=0.8, linestyle=:solid,
    )

    p = [0.1, 0.4, 0.5]
    px, py = TernaryDiagrams.from_bary_to_cart(p...)
    s0x, s0y = TernaryDiagrams.from_bary_to_cart(0, 0, 1)
    s1x, s1y = TernaryDiagrams.from_bary_to_cart(airline_oxygen, airline_nitrogen, 1 - airline_oxygen - airline_nitrogen)
    a = atan(s0y - s1y, s0x - s1x)



    return ax, px, py, a
end

function plot_stoichiometry(ax, element_counts)
    stoichiometric_O2 = element_counts["H"] * 0.25 + element_counts["C"] * 1.0
    if "O" in keys(element_counts)
        stoichiometric_O2 -= element_counts["O"] * 0.5
    end
    if stoichiometric_O2 <= 0
        @warn("Negative stoichiometric O2 requirement. Skipping stoichiometric line.")
        return
    end

    fuel_composition = 1 / (stoichiometric_O2 + 1)
    x0 = [0.0, 1.0, 0.0]
    x1 = [1 - fuel_composition, 0.0, fuel_composition]

    ternarylines!(ax, [x0[1], x1[1]], [x0[2], x1[2]], [x0[3], x1[3]];
        color=:blue, linewidth=0.8, linestyle=:solid,
    )

    p = x1 + [-0.05, 0.05, 0]
    px, py = TernaryDiagrams.from_bary_to_cart(p...)
    s0x, s0y = TernaryDiagrams.from_bary_to_cart(x1...)
    s1x, s1y = TernaryDiagrams.from_bary_to_cart(x0...)
    a = atan(s0y - s1y, s0x - s1x)

    println("Stoichiometric ", a)

    if a > 2.5
        text!(ax, px, py; text="Stoichiometric line", rotation=a - π, color=:blue, fontsize=tickfontsize)
    else
        text!(ax, px - 0.08, py - 0.065; text="Stoichiometric line", rotation=a - π, color=:blue, fontsize=tickfontsize)
    end
end

function plot_contour(ax, contour)
    color = :black
    map(contour) do l
        lines!(l; color, linewidth=0.8)
    end
end

function plot_lfl(ax, lel)
    xs, ys, zs = lel
    ternaryscatter!(ax, xs, ys, zs; color=:purple, marker=:xcross, markersize=7)
end

function plot_fl(ax, lel, uel)
    xs, ys, zs = zip(lel, uel)
    ternaryscatter!(ax, xs, ys, zs; color=:purple, marker=:xcross, markersize=7)
end

function set_arial_theme!()
    set_theme!(Theme(fonts=(; regular="Arial", bold="Arial Bold", italic="Arial Italic")))
    return nothing
end

chem_rich(sym::AbstractString, sub::AbstractString) = rich(sym, subscript(sub))

function plot_PD_scatter(fig, ax, fuel_name, a1, a2, a3, mus, level, withcolorbar, color=:black)
    Tmin, Tmax = 0, 3600
    colormap = :jet

    xlims!(ax, -0.05, 1.06)
    ylims!(ax, -0.15, 0.95)

    hidedecorations!(ax)
    hidespines!(ax)

    p = ternaryscatter!(ax, a1, a2, a3;
        color=mus,
        colormap,
        colorrange=[Tmin, Tmax],
        clipaxes=true,
    )

    if withcolorbar == true
        Colorbar(fig[2, 1], p;
            vertical=false,
            width=Relative(4 / 5),
            ticks=[Tmin, 1200, 2400, 3600],
            height=4,
            ticklabelsize=labelfontsize - 1.7,
            flipaxis=false,
            spinewidth=0.5,
            ticksize=2,
            tickwidth=0.5,
        )
    end

    ternaryaxis!(ax;
        labelx_arrow="",
        labely_arrow=chem_rich("N", "2"),
        labelz_arrow=fuel_name,
        hide_vertex_labels=true,
        hide_arrows=true,
        tick_fontsize=tickfontsize,
        label_edge_vertical_adjustment=0.16,
        arrow_label_fontsize=labelfontsize,
        arrow_label_rotation_adjustment=1,
    )

    Label(fig[1, 1, Bottom()], chem_rich("O", "2");
        fontsize=labelfontsize,
        tellheight=false,
        padding=(0, 0, 29, 0)
    )

    levels_line = [level]
    ternarycontour!(ax, a1, a2, a3, mus;
        levels=levels_line,
        color,
        linewidth=0.8,
        colormap,
    )
    if withcolorbar == true
        rowsize!(fig.layout, 2, Fixed(12))
        rowgap!(fig.layout, 1, 2)
    end
    return nothing
end

# =========================
# Main plot function
# =========================
"""
Plot phase diagram and annotate Hf.
NOTE: This function only plots and returns fig. It does NOT write any CSV.
"""
function plot_PD(Hf_csv, dat_filename, level, withcolorbar, output_pdf)

    # withcolorbar = true

    if withcolorbar == false
        fact = 1.22
        fig = Figure(size=(171 * fact, 188 * fact))
    else
        fig = Figure(size=(171 * 1.27, 215 * 1.27))
    end

    # ax = Axis(fig[1, 1], aspect=AxisAspect(1))
    ax = Axis(fig[1, 1], aspect=DataAspect())
    # ax.alignmode = Inside()

    res = compute_PD_metrics(dat_filename, level)

    fuel_name = rich_formula_rt(res.formula)
    plot_PD_scatter(fig, ax, fuel_name, res.a1, res.a2, res.a3, res.mus, level, withcolorbar)

    ax, px, py, a = plot_airline(ax)


    element_counts = count_elements(res.formula)

    if res.t !== nothing
        plot_contour(ax, res.t)
        if any(occursin(x, dat_filename) for x in ["CHNO3_0b6c"])
            @warn("Not plot stoichiometry.")
        else
            plot_stoichiometry(ax, element_counts)
        end
        if res.ufl[3] - res.lfl[3] < 1e-3 || res.ufl[1] < 1e-3
            plot_lfl(ax, res.lfl)
        else
            plot_fl(ax, res.lfl, res.ufl)
        end
    end

    text!(ax, px, py; text="Air line", rotation=a - π, color=:red, fontsize=tickfontsize)


    Hf_text = find_Hf_298K(res.jsonPrefix, Hf_csv)

    text_loc = [-0.2, 0.35, 0.82]
    px, py = TernaryDiagrams.from_bary_to_cart(text_loc...)
    lfl = copy(res.lfl)
    ufl = copy(res.ufl)

    if abs(lfl[3]) < 1e-6
        lfl[3] = 0.0
    end
    if abs(ufl[3]) < 1e-6
        ufl[3] = 0.0
    end

    lel_text = @sprintf("%.1f", lfl[3] * 100)
    uel_text = @sprintf("%.1f", ufl[3] * 100)

    if any(occursin(x, dat_filename) for x in ["CH2O3_08bda0f3"])
        text_str = [
            rich(
                "Δ",
                rich("H"; font=:italic),
                rich(
                    subscript("f", offset=(0.1, 0)),
                    superscript("°", offset=(-0.1, 0));
                ),
                " = $(Hf_text)"
            ),
            rich("LFL = N/A"),
            rich("UFL = N/A"),
        ]
    else
        text_str = Union{String,Makie.RichText}[
            rich(
                "Δ",
                rich("H"; font=:italic),
                rich(
                    subscript("f", offset=(0.1, 0)),
                    superscript("°", offset=(-0.1, 0));
                ),
                " = $(Hf_text)"
            ),
            "LFL = $(lel_text)%",
            "UFL = $(uel_text)%"
        ]
    end

    text!(ax,
        [px - 0.02 for _ in 1:length(text_str)],
        [py + 0.18 - i * 0.07 for i in 0:length(text_str)-1];
        text=text_str,
        fontsize=labelfontsize - 1,
        font="Helvetica",
    )
    xlims!(ax, -0.12, 1.1)
    ylims!(ax, -0.18, 1.04)
    save(output_pdf, fig)
    println("Saved plot to: $(abspath(output_pdf))")

    return fig

end

end # module CAFTMetrics
