# -*- coding: utf-8 -*-

"""
    CVExtraction plotting

Visualization helpers for CV parameter extraction results.
"""

using CairoMakie

"""
    plot_cv_analysis(result::CVExtractionResult; title::String = "CV Analysis") -> Figure

Plot the original (raw) and corrected mean cycles side by side.
"""
function plot_cv_analysis(result::CVExtractionResult; title::String = "CV Analysis")::Figure
    fig = Figure(size = (1200, 500))

    # Left: original CV
    ax_raw = Axis(
        fig[1, 1],
        title = "Original CV",
        xlabel = "U in V vs reference electrode",
        ylabel = "j in A cm⁻²",
    )
    raw = result.raw_mean_cycle
    lines!(ax_raw, raw.U, raw.j; color = :black, linewidth = 1.2, label = "Mean cycle")

    # Right: corrected CV
    ax_corr = Axis(
        fig[1, 2],
        title = "Corrected CV",
        xlabel = "U in V vs reference electrode",
        ylabel = "j in A cm⁻²",
    )
    mc = result.mean_cycle
    an = result.anodic
    cat = result.cathodic
    lines!(ax_corr, mc.U, mc.j; color = :black, linewidth = 1.2, label = "Mean cycle")
    lines!(ax_corr, an.U, an.j; color = :blue, linewidth = 0.8, label = "Anodic")
    lines!(ax_corr, cat.U, cat.j; color = :red, linewidth = 0.8, label = "Cathodic")
    axislegend(ax_corr; position = :rb)

    # Link axes for easy comparison
    linkxaxes!(ax_raw, ax_corr)
    linkyaxes!(ax_raw, ax_corr)

    return fig
end
