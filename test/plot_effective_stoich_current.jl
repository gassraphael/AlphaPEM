# -*- coding: utf-8 -*-

"""
Plots the effective current density `effective_stoich_current` used in `desired_flows`
(src/alphapem/core/models/fuelcell/velocity.jl) against the fuel cell current density i_fc_cell,
to visualize the minimum-flow floor at i_min_stoich and the smoothing at i_fc_cell = 0 and
i_fc_cell = i_min_stoich.
"""

using AlphaPEM
using GLMakie
using CairoMakie

const effective_stoich_current = AlphaPEM.Core.Models.effective_stoich_current

i_min_stoich = 0.5  # A.cm-2, same default value as OperatingConditions.i_min_stoich.

i_fc_cell_range = 0.0:1e-4:0.7  # A.cm-2
i_eff_range = effective_stoich_current.(i_fc_cell_range, i_min_stoich)

fig = Figure(size=(800, 500))
ax = Axis(fig[1, 1];
          xlabel="i_fc_cell (A.cm⁻²)",
          ylabel="i_eff (A.cm⁻²)",
          title="Effective current density used for the desired flows")

lines!(ax, i_fc_cell_range, i_eff_range; color=:royalblue, linewidth=2.8, label="i_eff")
lines!(ax, i_fc_cell_range, collect(i_fc_cell_range); color=:gray, linestyle=:dash, linewidth=1.5, label="y = i_fc_cell")
hlines!(ax, [i_min_stoich]; color=:red, linestyle=:dot, linewidth=1.5, label="i_min_stoich")

axislegend(ax; position=:rb)

GLMakie.activate!()
display(fig)
readline()
