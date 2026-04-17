# plot.jl
#
# Native CairoMakie plotting recipes for AlphaPEM (progressive migration stage).


using .PlotHelpers: _publication_colors,
                    _finalize_axis!,
                    _label_with_rmse,
                    _polarization_legend_base,
                    _experimental_marker,
                    _set_polarization_axis_limits!,
                    _set_polarization_fixed_ticks!,
                    _format_fixed,
                    _compact_tick_labels,
                    _nice_tick_step,
                    _colorbar_ticks_auto,
                    _rounded_major_ticks,
                    _set_dense_ticks!,
                    lsub

# ═══════════════════════════════════════════════════════════════════════════════
#  Internal States - Temporal Plots
# ═══════════════════════════════════════════════════════════════════════════════

"""Plot current density through the gas-channel direction as a function of time."""
function plot_ifc_1D_temporal(outputs::SimulationOutputs,
                              fc::AbstractFuelCell,
                              cd::AbstractCurrent,
                              cfg::SimulationConfig,
                              ax)
    palette = _publication_colors()
    t = masked_time_history(outputs, cd, cfg)
    nb_gc = fc.numerical_parameters.nb_gc

    i_fc_cell_t = current(cd, t) ./ 1e4
    lines!(ax, t, i_fc_cell_t; color=:black, linewidth=2.8, label=lsub("i", "fc,cell"))

    for i in 1:nb_gc
        i_fc_t = extract_masked_derived_gc_series(outputs, i, cd, cfg, x -> x.i_fc) ./ 1e4
        lines!(ax, t, i_fc_t; color=palette[mod1(i, length(palette))], label=lsub("i", "fc,$(i)"))
    end
    i_fc_series = [extract_masked_derived_gc_series(outputs, i, cd, cfg, x -> x.i_fc) ./ 1e4 for i in 1:nb_gc]
    _set_dense_ticks!(ax, t, vcat([i_fc_cell_t], i_fc_series))
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                    legend=true,
                    legend_position=:lb)
    return nothing
end

"""Plot vapour concentration through the MEA middle channel as a function of time."""
function plot_C_v_1D_temporal(outputs::SimulationOutputs,
                              fc::AbstractFuelCell,
                              cd::AbstractCurrent,
                              cfg::SimulationConfig,
                              ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(fc)
    nb_mpl_mid = middle_mpl_index(fc)
    t = masked_time_history(outputs, cd, cfg)
    T_ccl = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ccl.T)
    C_v_sat_ccl = [C_v_sat(T) for T in T_ccl]

    series = [
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agc.C_v), lsub("C", "v,agc")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agdl[nb_gdl_mid].C_v), lsub("C", "v,agdl")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ampl[nb_mpl_mid].C_v), lsub("C", "v,ampl")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.acl.C_v), lsub("C", "v,acl")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ccl.C_v), lsub("C", "v,ccl")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cmpl[nb_mpl_mid].C_v), lsub("C", "v,cmpl")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgdl[nb_gdl_mid].C_v), lsub("C", "v,cgdl")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgc.C_v), lsub("C", "v,cgc")),
    ]

    for (i, (y, lbl)) in enumerate(series)
        lines!(ax, t, y; color=palette[mod1(i, length(palette))], label=lbl)
    end
    # Reference saturation curve: bold black line
    lines!(ax, t, C_v_sat_ccl; color=:black, linewidth=2.8,
           label=rich(lsub("C", "v,sat,ccl"); font=:bold))
    _set_dense_ticks!(ax, t, vcat([s[1] for s in series], [C_v_sat_ccl]))
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Vapor concentration ", lsub("C", "v"), " (mol·m⁻³)"),
                    legend=true,
                    legend_position=:rt)
    return nothing
end

"""Plot membrane/CL water content as a function of time."""
function plot_lambda_1D_temporal(outputs::SimulationOutputs,
                                 fc::AbstractFuelCell,
                                 cd::AbstractCurrent,
                                 cfg::SimulationConfig,
                                 ax)
    palette = _publication_colors()
    if cfg.display_timing == :postrun && cfg.type_current isa Union{PolarizationParams, PolarizationCalibrationParams}
        t_hist = time_history(outputs)
        sample_indices = polarisation_sampling_indices(outputs, cd)
        ifc_discretized = [current(cd, t_hist[idx]) / 1e4 for idx in sample_indices]

        lambda_acl_t = extract_mid_mea_series(outputs, fc, mea -> mea.acl.lambda)[sample_indices]
        lambda_mem_t = extract_mid_mea_series(outputs, fc, mea -> mea.mem.lambda)[sample_indices]
        lambda_ccl_t = extract_mid_mea_series(outputs, fc, mea -> mea.ccl.lambda)[sample_indices]

        scatter!(ax, ifc_discretized, lambda_acl_t; color=palette[3], markersize=8, label=lsub("λ", "acl"))
        scatter!(ax, ifc_discretized, lambda_mem_t; color=palette[4], markersize=8, label=lsub("λ", "mem"))
        scatter!(ax, ifc_discretized, lambda_ccl_t; color=palette[5], markersize=8, label=lsub("λ", "ccl"))

        _set_dense_ticks!(ax, ifc_discretized, [lambda_acl_t, lambda_mem_t, lambda_ccl_t])
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Water content ", lsub("λ", ""), " (-)"),
                        legend=true,
                        legend_position=:lb)
    else
        t = masked_time_history(outputs, cd, cfg)
        lines!(ax, t, extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.acl.lambda);
               color=palette[3], label=lsub("λ", "acl"))
        lines!(ax, t, extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.mem.lambda);
               color=palette[4], label=lsub("λ", "mem"))
        lines!(ax, t, extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ccl.lambda);
               color=palette[5], label=lsub("λ", "ccl"))
        lambda_acl_t = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.acl.lambda)
        lambda_mem_t = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.mem.lambda)
        lambda_ccl_t = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ccl.lambda)
        _set_dense_ticks!(ax, t, [lambda_acl_t, lambda_mem_t, lambda_ccl_t])
        _finalize_axis!(ax;
                        xlabel=rich("Time ", lsub("t", ""), " (s)"),
                        ylabel=rich("Water content ", lsub("λ", ""), " (-)"),
                        legend=true,
                        legend_position=:lb)
    end
    return nothing
end

"""Plot liquid water saturation as a function of time."""
function plot_s_1D_temporal(outputs::SimulationOutputs,
                            fc::AbstractFuelCell,
                            cd::AbstractCurrent,
                            cfg::SimulationConfig,
                            ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(fc)
    nb_mpl_mid = middle_mpl_index(fc)
    if cfg.display_timing == :postrun && cfg.type_current isa Union{PolarizationParams, PolarizationCalibrationParams}
        t_hist = time_history(outputs)
        sample_indices = polarisation_sampling_indices(outputs, cd)
        ifc_discretized = [current(cd, t_hist[idx]) / 1e4 for idx in sample_indices]

        series = [
            (extract_mid_mea_series(outputs, fc, mea -> mea.agc.s)[sample_indices], lsub("s", "agc")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.agdl[nb_gdl_mid].s)[sample_indices], lsub("s", "agdl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.ampl[nb_mpl_mid].s)[sample_indices], lsub("s", "ampl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.acl.s)[sample_indices], lsub("s", "acl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.ccl.s)[sample_indices], lsub("s", "ccl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.cmpl[nb_mpl_mid].s)[sample_indices], lsub("s", "cmpl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.cgdl[nb_gdl_mid].s)[sample_indices], lsub("s", "cgdl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.cgc.s)[sample_indices], lsub("s", "cgc")),
        ]

        for (i, (y, lbl)) in enumerate(series)
            scatter!(ax, ifc_discretized, y; color=palette[mod1(i, length(palette))], markersize=8, label=lbl)
        end
        _set_dense_ticks!(ax, ifc_discretized, [s[1] for s in series])
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Liquid saturation ", lsub("s", ""), " (-)"),
                        legend=true,
                        legend_position=:rt)
    else
        t = masked_time_history(outputs, cd, cfg)

        series = [
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agc.s), lsub("s", "agc")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agdl[nb_gdl_mid].s), lsub("s", "agdl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ampl[nb_mpl_mid].s), lsub("s", "ampl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.acl.s), lsub("s", "acl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ccl.s), lsub("s", "ccl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cmpl[nb_mpl_mid].s), lsub("s", "cmpl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgdl[nb_gdl_mid].s), lsub("s", "cgdl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgc.s), lsub("s", "cgc")),
        ]

        for (i, (y, lbl)) in enumerate(series)
            lines!(ax, t, y; color=palette[mod1(i, length(palette))], label=lbl)
        end
        _set_dense_ticks!(ax, t, [s[1] for s in series])
        _finalize_axis!(ax;
                        xlabel=rich("Time ", lsub("t", ""), " (s)"),
                        ylabel=rich("Liquid saturation ", lsub("s", ""), " (-)"),
                        legend=true,
                        legend_position=:rt)
    end
    return nothing
end

"""Plot hydrogen concentration as a function of time."""
function plot_C_H2_1D_temporal(outputs::SimulationOutputs,
                               fc::AbstractFuelCell,
                               cd::AbstractCurrent,
                               cfg::SimulationConfig,
                               ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(fc)
    nb_mpl_mid = middle_mpl_index(fc)
    t = masked_time_history(outputs, cd, cfg)

    lines!(ax, t, extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agc.C_H2);
           color=palette[1], label=lsub("C", "H₂,agc"))
    lines!(ax, t, extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agdl[nb_gdl_mid].C_H2);
           color=palette[2], label=lsub("C", "H₂,agdl"))
    lines!(ax, t, extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ampl[nb_mpl_mid].C_H2);
           color=palette[3], label=lsub("C", "H₂,ampl"))
    lines!(ax, t, extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.acl.C_H2);
           color=palette[4], label=lsub("C", "H₂,acl"))

    C_H2_agc_t = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agc.C_H2)
    C_H2_agdl_t = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agdl[nb_gdl_mid].C_H2)
    C_H2_ampl_t = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ampl[nb_mpl_mid].C_H2)
    C_H2_acl_t = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.acl.C_H2)
    _set_dense_ticks!(ax, t, [C_H2_agc_t, C_H2_agdl_t, C_H2_ampl_t, C_H2_acl_t])
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Hydrogen concentration ", lsub("C", "H₂"), " (mol·m⁻³)"),
                    legend=true,
                    legend_position=:rt)
    return nothing
end

"""Plot oxygen concentration as a function of time."""
function plot_C_O2_1D_temporal(outputs::SimulationOutputs,
                               fc::AbstractFuelCell,
                               cd::AbstractCurrent,
                               cfg::SimulationConfig,
                               ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(fc)
    nb_mpl_mid = middle_mpl_index(fc)
    t = masked_time_history(outputs, cd, cfg)

    lines!(ax, t, extract_masked_mid_derived_gc_series(outputs, fc, cd, cfg, x -> x.C_O2_Pt);
           color=palette[10], label=lsub("C", "O₂,Pt"))
    lines!(ax, t, extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ccl.C_O2);
           color=palette[4], label=lsub("C", "O₂,ccl"))
    lines!(ax, t, extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cmpl[nb_mpl_mid].C_O2);
           color=palette[5], label=lsub("C", "O₂,cmpl"))
    lines!(ax, t, extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgdl[nb_gdl_mid].C_O2);
           color=palette[6], label=lsub("C", "O₂,cgdl"))
    lines!(ax, t, extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgc.C_O2);
           color=palette[7], label=lsub("C", "O₂,cgc"))

    C_O2_Pt_t = extract_masked_mid_derived_gc_series(outputs, fc, cd, cfg, x -> x.C_O2_Pt)
    C_O2_ccl_t = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ccl.C_O2)
    C_O2_cmpl_t = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cmpl[nb_mpl_mid].C_O2)
    C_O2_cgdl_t = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgdl[nb_gdl_mid].C_O2)
    C_O2_cgc_t = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgc.C_O2)
    _set_dense_ticks!(ax, t, [C_O2_Pt_t, C_O2_ccl_t, C_O2_cmpl_t, C_O2_cgdl_t, C_O2_cgc_t])
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Oxygen concentration ", lsub("C", "O₂"), " (mol·m⁻³)"),
                    legend=true,
                    legend_position=:lb)
    return nothing
end

"""Plot temperature through the MEA as a function of time."""
function plot_T_1D_temporal(outputs::SimulationOutputs,
                            fc::AbstractFuelCell,
                            cd::AbstractCurrent,
                            cfg::SimulationConfig,
                            ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(fc)
    nb_mpl_mid = middle_mpl_index(fc)
    T_des = fc.operating_conditions.T_des - 273.15
    if cfg.display_timing == :postrun && cfg.type_current isa Union{PolarizationParams, PolarizationCalibrationParams}
        t_hist = time_history(outputs)
        sample_indices = polarisation_sampling_indices(outputs, cd)
        ifc_discretized = [current(cd, t_hist[idx]) / 1e4 for idx in sample_indices]

        series = [
            (extract_mid_mea_series(outputs, fc, mea -> mea.agc.T)[sample_indices] .- 273.15, lsub("T", "agc")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.agdl[nb_gdl_mid].T)[sample_indices] .- 273.15, lsub("T", "agdl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.ampl[nb_mpl_mid].T)[sample_indices] .- 273.15, lsub("T", "ampl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.acl.T)[sample_indices] .- 273.15, lsub("T", "acl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.mem.T)[sample_indices] .- 273.15, lsub("T", "mem")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.ccl.T)[sample_indices] .- 273.15, lsub("T", "ccl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.cmpl[nb_mpl_mid].T)[sample_indices] .- 273.15, lsub("T", "cmpl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.cgdl[nb_gdl_mid].T)[sample_indices] .- 273.15, lsub("T", "cgdl")),
            (extract_mid_mea_series(outputs, fc, mea -> mea.cgc.T)[sample_indices] .- 273.15, lsub("T", "cgc")),
        ]

        for (i, (y, lbl)) in enumerate(series)
            scatter!(ax, ifc_discretized, y; color=palette[mod1(i, length(palette))], markersize=8, label=lbl)
        end
        lines!(ax, [minimum(ifc_discretized), maximum(ifc_discretized)], [T_des, T_des];
               color=:black, linestyle=:dash, label=lsub("T", "des"))

        _set_dense_ticks!(ax, ifc_discretized, [s[1] for s in series])
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Temperature ", lsub("T", ""), " (°C)"),
                        legend=true,
                        legend_position=:lt)
    else
        t = masked_time_history(outputs, cd, cfg)

        series = [
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agc.T) .- 273.15, lsub("T", "agc")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agdl[nb_gdl_mid].T) .- 273.15, lsub("T", "agdl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ampl[nb_mpl_mid].T) .- 273.15, lsub("T", "ampl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.acl.T) .- 273.15, lsub("T", "acl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.mem.T) .- 273.15, lsub("T", "mem")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.ccl.T) .- 273.15, lsub("T", "ccl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cmpl[nb_mpl_mid].T) .- 273.15, lsub("T", "cmpl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgdl[nb_gdl_mid].T) .- 273.15, lsub("T", "cgdl")),
            (extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgc.T) .- 273.15, lsub("T", "cgc")),
        ]

        for (i, (y, lbl)) in enumerate(series)
            lines!(ax, t, y; color=palette[mod1(i, length(palette))], label=lbl)
        end
        lines!(ax, t, fill(T_des, length(t)); color=:black, linestyle=:dash, label=lsub("T", "des"))

        _set_dense_ticks!(ax, t, [s[1] for s in series])
        _finalize_axis!(ax;
                        xlabel=rich("Time ", lsub("t", ""), " (s)"),
                        ylabel=rich("Temperature ", lsub("T", ""), " (°C)"),
                        legend=true,
                        legend_position=:lt)
    end
    return nothing
end

"""Plot cell voltage as a function of time."""
function plot_Ucell(outputs::SimulationOutputs,
                    cd::AbstractCurrent,
                    cfg::SimulationConfig,
                    ax)
    t = masked_time_history(outputs, cd, cfg)
    Ucell_t = extract_masked_derived_series(outputs, cd, cfg, x -> x.Ucell)
    lines!(ax, t, Ucell_t; color=:black, label=lsub("U", "cell"))
    _set_dense_ticks!(ax, t, [Ucell_t])
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Cell voltage ", lsub("U", "cell"), " (V)"),
                    legend=true,
                    legend_position=:lt)
    return nothing
end

"""Plot nitrogen concentration at AGC and CGC for each GC node as a function of time."""
function plot_C_N2_1D_temporal(outputs::SimulationOutputs,
                               fc::AbstractFuelCell,
                               cd::AbstractCurrent,
                               cfg::SimulationConfig,
                               ax)
    palette = _publication_colors()
    nb_gc = fc.numerical_parameters.nb_gc
    t = masked_time_history(outputs, cd, cfg)

    C_N2_series = Vector{Vector{Float64}}()
    for i in 1:nb_gc
        C_N2_agc_i = extract_masked_mea_series(outputs, i, cd, cfg, mea -> mea.agc.C_N2)
        C_N2_cgc_i = extract_masked_mea_series(outputs, i, cd, cfg, mea -> mea.cgc.C_N2)
        push!(C_N2_series, C_N2_agc_i, C_N2_cgc_i)
        lines!(ax, t, C_N2_agc_i; color=palette[mod1(i, length(palette))],
               linestyle=:solid, label=lsub("C", "N₂,agc,$(i)"))
        lines!(ax, t, C_N2_cgc_i; color=palette[mod1(i, length(palette))],
               linestyle=:dash,  label=lsub("C", "N₂,cgc,$(i)"))
    end
    _set_dense_ticks!(ax, t, C_N2_series)
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Nitrogen concentration ", lsub("C", "N₂"), " (mol·m⁻³)"),
                    legend=true,
                    legend_position=:rt)
    return nothing
end

"""Plot relative humidity through the anode layers at the middle GC node as a function of time."""
function plot_Phi_a_1D_temporal(outputs::SimulationOutputs,
                                fc::AbstractFuelCell,
                                cd::AbstractCurrent,
                                cfg::SimulationConfig,
                                ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(fc)
    nb_mpl_mid = middle_mpl_index(fc)
    t = masked_time_history(outputs, cd, cfg)

    series = [
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg,
             mea -> mea.agc.C_v / C_v_sat(mea.agc.T)),                                    lsub("Φ", "a,agc")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg,
             mea -> mea.agdl[nb_gdl_mid].C_v / C_v_sat(mea.agdl[nb_gdl_mid].T)),          lsub("Φ", "a,agdl")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg,
             mea -> mea.ampl[nb_mpl_mid].C_v / C_v_sat(mea.ampl[nb_mpl_mid].T)),          lsub("Φ", "a,ampl")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg,
             mea -> mea.acl.C_v / C_v_sat(mea.acl.T)),                                    lsub("Φ", "a,acl")),
    ]

    for (i, (y, lbl)) in enumerate(series)
        lines!(ax, t, y; color=palette[mod1(i, length(palette))], label=lbl)
    end
    # Saturation threshold: Φ = 1.
    lines!(ax, t, ones(length(t));
           color=:black, linestyle=:dash, linewidth=1.4, label=rich("Φ = 1 (saturation)"))
    _set_dense_ticks!(ax, t, vcat([s[1] for s in series], [ones(length(t))]))
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Anode relative humidity ", lsub("Φ", "a"), " (–)"),
                    legend=true,
                    legend_position=:rt)
    return nothing
end

"""Plot relative humidity through the cathode layers at the middle GC node as a function of time."""
function plot_Phi_c_1D_temporal(outputs::SimulationOutputs,
                                fc::AbstractFuelCell,
                                cd::AbstractCurrent,
                                cfg::SimulationConfig,
                                ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(fc)
    nb_mpl_mid = middle_mpl_index(fc)
    t = masked_time_history(outputs, cd, cfg)

    series = [
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg,
             mea -> mea.ccl.C_v / C_v_sat(mea.ccl.T)),                                    lsub("Φ", "c,ccl")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg,
             mea -> mea.cmpl[nb_mpl_mid].C_v / C_v_sat(mea.cmpl[nb_mpl_mid].T)),          lsub("Φ", "c,cmpl")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg,
             mea -> mea.cgdl[nb_gdl_mid].C_v / C_v_sat(mea.cgdl[nb_gdl_mid].T)),          lsub("Φ", "c,cgdl")),
        (extract_masked_mid_mea_series(outputs, fc, cd, cfg,
             mea -> mea.cgc.C_v / C_v_sat(mea.cgc.T)),                                    lsub("Φ", "c,cgc")),
    ]

    for (i, (y, lbl)) in enumerate(series)
        lines!(ax, t, y; color=palette[mod1(i, length(palette))], label=lbl)
    end
    # Saturation threshold: Φ = 1.
    lines!(ax, t, ones(length(t));
           color=:black, linestyle=:dash, linewidth=1.4, label=rich("Φ = 1 (saturation)"))
    _set_dense_ticks!(ax, t, vcat([s[1] for s in series], [ones(length(t))]))
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Cathode relative humidity ", lsub("Φ", "c"), " (–)"),
                    legend=true,
                    legend_position=:rt)
    return nothing
end

"""Plot anode and cathode gas velocities for each GC node as a function of time."""
function plot_v_1D_temporal(outputs::SimulationOutputs,
                            fc::AbstractFuelCell,
                            cd::AbstractCurrent,
                            cfg::SimulationConfig,
                            ax)
    palette = _publication_colors()
    nb_gc = fc.numerical_parameters.nb_gc
    t = masked_time_history(outputs, cd, cfg)

    v_series = Vector{Vector{Float64}}()
    for i in 1:nb_gc
        v_a_i = extract_masked_derived_gc_series(outputs, i, cd, cfg, x -> x.v_a)
        v_c_i = extract_masked_derived_gc_series(outputs, i, cd, cfg, x -> x.v_c)
        push!(v_series, v_a_i, v_c_i)
        lines!(ax, t, v_a_i; color=palette[mod1(i, length(palette))],
               linestyle=:solid, label=lsub("v", "a,$(i)"))
        lines!(ax, t, v_c_i; color=palette[mod1(i, length(palette))],
               linestyle=:dash,  label=lsub("v", "c,$(i)"))
    end
    _set_dense_ticks!(ax, t, v_series)
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Gas velocity ", lsub("v", ""), " (m·s⁻¹)"),
                    legend=true,
                    legend_position=:rt)
    return nothing
end

"""Plot Reynolds number at each GC node for both the AGC and CGC as a function of time.

The hydraulic diameter is computed as Dh = 2HW/(H+W) for a rectangular cross-section.
Gas-mixture density is derived from ideal-gas law; dynamic viscosity uses the
Wilke/linear mixture rule from `mu_mixture_gases`.
"""
function plot_Re_nb_1D_temporal(outputs::SimulationOutputs,
                                fc::AbstractFuelCell,
                                cd::AbstractCurrent,
                                cfg::SimulationConfig,
                                ax)
    palette = _publication_colors()
    nb_gc = fc.numerical_parameters.nb_gc
    pp = fc.physical_parameters
    t = masked_time_history(outputs, cd, cfg)
    n_t = length(t)

    # Hydraulic diameters of the rectangular gas channels.
    Dh_a = 2 * pp.Hagc * pp.Wagc / (pp.Hagc + pp.Wagc)
    Dh_c = 2 * pp.Hcgc * pp.Wcgc / (pp.Hcgc + pp.Wcgc)

    Re_series = Vector{Vector{Float64}}()
    for i in 1:nb_gc
        # Gas velocities.
        v_a_i = extract_masked_derived_gc_series(outputs, i, cd, cfg, x -> x.v_a)
        v_c_i = extract_masked_derived_gc_series(outputs, i, cd, cfg, x -> x.v_c)
        # Anode GC state.
        C_v_agc  = extract_masked_mea_series(outputs, i, cd, cfg, mea -> mea.agc.C_v)
        C_H2_agc = extract_masked_mea_series(outputs, i, cd, cfg, mea -> mea.agc.C_H2)
        T_agc    = extract_masked_mea_series(outputs, i, cd, cfg, mea -> mea.agc.T)
        # Cathode GC state.
        C_v_cgc  = extract_masked_mea_series(outputs, i, cd, cfg, mea -> mea.cgc.C_v)
        C_O2_cgc = extract_masked_mea_series(outputs, i, cd, cfg, mea -> mea.cgc.C_O2)
        C_N2_cgc = extract_masked_mea_series(outputs, i, cd, cfg, mea -> mea.cgc.C_N2)
        T_cgc    = extract_masked_mea_series(outputs, i, cd, cfg, mea -> mea.cgc.T)

        Re_a_i = Vector{Float64}(undef, n_t)
        Re_c_i = Vector{Float64}(undef, n_t)
        for j in 1:n_t
            # Anode: H₂O + H₂ mixture (same approximation as flow modules).
            x_H2O_a = C_v_agc[j] + C_H2_agc[j] > 0 ?
                      C_v_agc[j] / (C_v_agc[j] + C_H2_agc[j]) : 0.0
            rho_a   = C_v_agc[j] * M_H2O + C_H2_agc[j] * M_H2   # kg·m⁻³
            mu_a    = mu_mixture_gases(["H2O_v", "H2"], [x_H2O_a, 1 - x_H2O_a], T_agc[j])
            Re_a_i[j] = mu_a > 0 ? rho_a * abs(v_a_i[j]) * Dh_a / mu_a : 0.0

            # Cathode: H₂O + O₂ + N₂ mixture.
            C_dry_c  = max(C_O2_cgc[j] + C_N2_cgc[j], eps(Float64))
            y_O2_c   = C_O2_cgc[j] / C_dry_c
            x_H2O_c  = C_v_cgc[j] + C_O2_cgc[j] + C_N2_cgc[j] > 0 ?
                       C_v_cgc[j] / (C_v_cgc[j] + C_O2_cgc[j] + C_N2_cgc[j]) : 0.0
            x_O2_c   = y_O2_c * (1 - x_H2O_c)
            x_N2_c   = (1 - y_O2_c) * (1 - x_H2O_c)
            rho_c    = C_v_cgc[j] * M_H2O + C_O2_cgc[j] * M_O2 + C_N2_cgc[j] * M_N2
            mu_c     = mu_mixture_gases(["H2O_v", "O2", "N2"],
                                        [x_H2O_c, x_O2_c, x_N2_c], T_cgc[j])
            Re_c_i[j] = mu_c > 0 ? rho_c * abs(v_c_i[j]) * Dh_c / mu_c : 0.0
        end
        push!(Re_series, Re_a_i, Re_c_i)
        lines!(ax, t, Re_a_i; color=palette[mod1(i, length(palette))],
               linestyle=:solid, label=lsub("Re", "a,$(i)"))
        lines!(ax, t, Re_c_i; color=palette[mod1(i, length(palette))],
               linestyle=:dash,  label=lsub("Re", "c,$(i)"))
    end
    _set_dense_ticks!(ax, t, Re_series)
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Reynolds number ", lsub("Re", ""), " (–)"),
                    legend=true,
                    legend_position=:rt)
    return nothing
end

"""Plot pressure at middle AGC/CGC and inlets/outlets as a function of time."""
function plot_P_1D_temporal(outputs::SimulationOutputs,
                            fc::AbstractFuelCell,
                            cd::AbstractCurrent,
                            cfg::SimulationConfig,
                            ax)
    palette = _publication_colors()
    Pa_des = fc.operating_conditions.Pa_des
    Pc_des = fc.operating_conditions.Pc_des
    t = masked_time_history(outputs, cd, cfg)

    C_v_agc = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agc.C_v)
    C_H2_agc = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agc.C_H2)
    C_N2_agc = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agc.C_N2)
    T_agc = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.agc.T)
    C_v_cgc = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgc.C_v)
    C_O2_cgc = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgc.C_O2)
    C_N2_cgc = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgc.C_N2)
    T_cgc = extract_masked_mid_mea_series(outputs, fc, cd, cfg, mea -> mea.cgc.T)

    P_agc = (C_v_agc .+ C_H2_agc .+ C_N2_agc) .* R .* T_agc ./ 1e5
    P_cgc = (C_v_cgc .+ C_O2_cgc .+ C_N2_cgc) .* R .* T_cgc ./ 1e5
    Pa_in = extract_masked_derived_series(outputs, cd, cfg, x -> x.Pa_in) ./ 1e5
    Pc_in = extract_masked_derived_series(outputs, cd, cfg, x -> x.Pc_in) ./ 1e5
    Pa_out = fill(Pa_des / 1e5, length(t))
    Pc_out = fill(Pc_des / 1e5, length(t))

    lines!(ax, t, P_agc; color=palette[1], label=lsub("P", "agc"))
    lines!(ax, t, P_cgc; color=palette[2], label=lsub("P", "cgc"))
    lines!(ax, t, Pa_in; color=palette[3], label=lsub("P", "a,in"))
    lines!(ax, t, Pa_out; color=palette[4], label=lsub("P", "a,out"))
    lines!(ax, t, Pc_in; color=palette[5], label=lsub("P", "c,in"))
    lines!(ax, t, Pc_out; color=palette[6], label=lsub("P", "c,out"))

    _set_dense_ticks!(ax, t, [P_agc, P_cgc, Pa_in, Pa_out, Pc_in, Pc_out])
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Pressure ", lsub("P", ""), " (bar)"),
                    legend=true,
                    legend_position=:lb)
    return nothing
end


# ═══════════════════════════════════════════════════════════════════════════════
#  Polarization Curves
# ═══════════════════════════════════════════════════════════════════════════════

"""Plot polarization curve with unified CairoMakie style conventions."""
function plot_polarization_curve(outputs::SimulationOutputs,
                                 fc::AbstractFuelCell,
                                 cd::AbstractCurrent,
                                 cfg::SimulationConfig,
                                 ax)
    palette = _publication_colors()
    model_color = palette[1]
    exp_color = RGBf(0.10, 0.10, 0.10)
    model_label_base = _polarization_legend_base(cfg.type_fuel_cell;
                                                 simulation=true,
                                                 calibration=false)
    exp_label = _polarization_legend_base(cfg.type_fuel_cell;
                                          simulation=false,
                                          calibration=false)

    if cfg.display_timing == :postrun
        ifc_discretized, Ucell_discretized = _polarization_points(outputs, cd)

        y_series = [Ucell_discretized]
        sim_error = nothing
        if hasproperty(fc, :pola_exp_data)
            exp_data = getproperty(fc, :pola_exp_data)
            if hasproperty(exp_data, :i_exp) && hasproperty(exp_data, :U_exp)
                i_exp = getproperty(exp_data, :i_exp) ./ 1e4
                U_exp = getproperty(exp_data, :U_exp)
                scatter!(ax, i_exp, U_exp;
                         color=exp_color, marker=_experimental_marker(cfg.type_fuel_cell), markersize=12,
                         strokecolor=:white, strokewidth=0.5,
                         label=exp_label)
                push!(y_series, U_exp)
                _pola_rmse_enabled(cfg) && (sim_error = _polarization_rmse(ifc_discretized, Ucell_discretized, i_exp, U_exp))
            end
        end

        model_label = _label_with_rmse(model_label_base, sim_error)
        lines!(ax, ifc_discretized, Ucell_discretized;
               color=model_color, linewidth=3.0, label=model_label)

        _set_polarization_axis_limits!(ax)
        _set_polarization_fixed_ticks!(ax)
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Cell voltage ", lsub("U", "cell"), " (V)"),
                        title="Polarization curve",
                        legend=true,
                        legend_position=:rt)
    else
        # Dynamic mode: display the latest simulated operating point.
        t_hist = time_history(outputs)
        Ucell_t = derived_outputs(outputs).Ucell
        idx = lastindex(t_hist)
        ifc_last = current(cd, t_hist[idx]) / 1e4
        Ucell_last = Ucell_t[idx]
        scatter!(ax, [ifc_last], [Ucell_last];
                 color=:white, marker=:circle, markersize=8,
                 strokecolor=model_color, strokewidth=1.8,
                 label=model_label_base)

        _set_polarization_axis_limits!(ax)
        _set_polarization_fixed_ticks!(ax)
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Cell voltage ", lsub("U", "cell"), " (V)"),
                        title="Polarization point (dynamic)",
                        legend=true,
                        legend_position=:rt)
    end
    return nothing
end


"""Plot calibration polarization curve with unified CairoMakie style conventions."""
function plot_polarization_curve_for_cali(outputs::SimulationOutputs,
                                          fc::AbstractFuelCell,
                                          cd::AbstractCurrent,
                                          cfg::SimulationConfig,
                                          ax)
    palette = _publication_colors()
    model_color = palette[1]
    exp_color = RGBf(0.10, 0.10, 0.10)
    model_label_base = _polarization_legend_base(cfg.type_fuel_cell;
                                                 simulation=true,
                                                 calibration=true)
    exp_label = _polarization_legend_base(cfg.type_fuel_cell;
                                          simulation=false,
                                          calibration=true)

    if cfg.display_timing == :postrun
        ifc_discretized, Ucell_discretized = _polarization_points(outputs, cd)

        y_series = [Ucell_discretized]
        sim_error = nothing
        if hasproperty(fc, :pola_exp_data_cali)
            exp_data = getproperty(fc, :pola_exp_data_cali)
            if hasproperty(exp_data, :i_exp) && hasproperty(exp_data, :U_exp)
                i_exp = getproperty(exp_data, :i_exp) ./ 1e4
                U_exp = getproperty(exp_data, :U_exp)
                scatter!(ax, i_exp, U_exp;
                         color=exp_color, marker=_experimental_marker(cfg.type_fuel_cell), markersize=12,
                         strokecolor=:white, strokewidth=0.5,
                         label=exp_label)
                push!(y_series, U_exp)
                # Legacy calibration logic: RMSE is computed directly point-to-point.
                _pola_rmse_enabled(cfg) && (sim_error = calculate_simulation_error(Ucell_discretized, U_exp))
            end
        end

        model_label = _label_with_rmse(model_label_base, sim_error)
        lines!(ax, ifc_discretized, Ucell_discretized;
               color=model_color, linewidth=3.0, label=model_label)

        _set_polarization_axis_limits!(ax)
        _set_polarization_fixed_ticks!(ax)
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Cell voltage ", lsub("U", "cell"), " (V)"),
                        title="Polarization curve (calibration)",
                        legend=true,
                        legend_position=:rt)
    else
        t_hist = time_history(outputs)
        Ucell_t = derived_outputs(outputs).Ucell
        idx = lastindex(t_hist)
        ifc_last = current(cd, t_hist[idx]) / 1e4
        Ucell_last = Ucell_t[idx]
        scatter!(ax, [ifc_last], [Ucell_last];
                 color=:white, marker=:circle, markersize=8,
                 strokecolor=model_color, strokewidth=1.8,
                 label=model_label_base)
        _set_polarization_axis_limits!(ax)
        _set_polarization_fixed_ticks!(ax)
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Cell voltage ", lsub("U", "cell"), " (V)"),
                        title="Polarization point (calibration, dynamic)",
                        legend=true,
                        legend_position=:rt)
    end
    return nothing
end



# ═══════════════════════════════════════════════════════════════════════════════
#  Derived Polarization Curves
# ═══════════════════════════════════════════════════════════════════════════════

"""Plot power density curve with unified CairoMakie style conventions.

Power density is defined as P = U_cell × i_fc, expressed in W·cm⁻²."""
function plot_power_density_curve(outputs::SimulationOutputs,
                                  fc::AbstractFuelCell,
                                  cd::AbstractCurrent,
                                  cfg::SimulationConfig,
                                  ax)
    model_label = "Model ($(String(cfg.type_fuel_cell)))"

    if cfg.display_timing == :postrun
        ifc_discretized, Ucell_discretized = _polarization_points(outputs, cd)
        P_discretized = ifc_discretized .* Ucell_discretized   # W·cm⁻²
        lines!(ax, ifc_discretized, P_discretized;
               color=:black, linewidth=2.6, label=model_label)
        scatter!(ax, ifc_discretized, P_discretized;
                 color=:black, markersize=7)
        _set_dense_ticks!(ax, ifc_discretized, [P_discretized])
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Power density ", lsub("P", "cell"), " (W·cm⁻²)"),
                        title="Power density curve",
                        legend=true,
                        legend_position=:lt)
    else
        # Dynamic mode: display the latest simulated operating point.
        t_hist   = time_history(outputs)
        Ucell_t  = derived_outputs(outputs).Ucell
        idx      = lastindex(t_hist)
        ifc_last = current(cd, t_hist[idx]) / 1e4
        P_last   = ifc_last * Ucell_t[idx]
        scatter!(ax, [ifc_last], [P_last]; color=:black, markersize=8, label=model_label)
        _set_dense_ticks!(ax, [ifc_last], [[P_last]])
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Power density ", lsub("P", "cell"), " (W·cm⁻²)"),
                        title="Power density point (dynamic)",
                        legend=true,
                        legend_position=:rt)
    end
    return nothing
end


"""Plot cell voltage efficiency curve with unified CairoMakie style conventions.

Voltage efficiency is defined as η_v = U_cell / E_th, where E_th = 1.481 V is the
thermoneutral voltage based on the higher heating value (HHV) of hydrogen."""
function plot_cell_efficiency(outputs::SimulationOutputs,
                              fc::AbstractFuelCell,
                              cd::AbstractCurrent,
                              cfg::SimulationConfig,
                              ax)
    E_th = 1.481   # V — HHV-based thermoneutral voltage for H₂/O₂ PEM fuel cell.
    model_label = "Model ($(String(cfg.type_fuel_cell)))"

    if cfg.display_timing == :postrun
        ifc_discretized, Ucell_discretized = _polarization_points(outputs, cd)
        eta_discretized = Ucell_discretized ./ E_th
        lines!(ax, ifc_discretized, eta_discretized;
               color=:black, linewidth=2.6, label=model_label)
        scatter!(ax, ifc_discretized, eta_discretized;
                 color=:black, markersize=7)
        # Reference line at η = 1 (upper thermodynamic bound).
        lines!(ax, [ifc_discretized[1], ifc_discretized[end]], [1.0, 1.0];
               color=:gray, linestyle=:dot, linewidth=1.4, label="η = 1")
        _set_dense_ticks!(ax, ifc_discretized, [eta_discretized])
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Voltage efficiency ", lsub("η", "v"), " (–)"),
                        title="Cell efficiency curve",
                        legend=true,
                        legend_position=:lb)
    else
        t_hist   = time_history(outputs)
        Ucell_t  = derived_outputs(outputs).Ucell
        idx      = lastindex(t_hist)
        ifc_last = current(cd, t_hist[idx]) / 1e4
        eta_last = Ucell_t[idx] / E_th
        scatter!(ax, [ifc_last], [eta_last]; color=:black, markersize=8, label=model_label)
        _set_dense_ticks!(ax, [ifc_last], [[eta_last]])
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Voltage efficiency ", lsub("η", "v"), " (–)"),
                        title="Cell efficiency point (dynamic)",
                        legend=true,
                        legend_position=:lb)
    end
    return nothing
end


# ═══════════════════════════════════════════════════════════════════════════════
#  EIS Curves
# ═══════════════════════════════════════════════════════════════════════════════

"""Plot EIS Nyquist diagram with unified style conventions."""
function plot_EIS_curve_Nyquist(cd::AbstractCurrent,
                                Fourier_results::FourierOutputs,
                                ax)
    Z_real, minus_Z_imag, _f, _abs_Z, _phi = _eis_point(cd, Fourier_results)
    isfinite(Z_real) && isfinite(minus_Z_imag) || return nothing

    scatter!(ax, [Z_real], [minus_Z_imag];
             color=_publication_colors()[1], markersize=7,
             strokecolor=:white, strokewidth=0.4)

    # Keep Nyquist geometry visually correct for publication use.
    ax.aspect = DataAspect()
    _finalize_axis!(ax;
                    xlabel=rich("Real impedance ", lsub("Z", "real"), " (mΩ·cm²)"),
                    ylabel=rich("Imaginary impedance ", rich("-", lsub("Z", "imag")), " (mΩ·cm²)"),
                    title="EIS - Nyquist",
                    legend=false)
    return nothing
end


"""Plot EIS Bode amplitude diagram with unified style conventions."""
function plot_EIS_curve_Bode_amplitude(cd::AbstractCurrent,
                                       Fourier_results::FourierOutputs,
                                       ax)
    _Z_real, _minus_Z_imag, f, abs_Z, _phi = _eis_point(cd, Fourier_results)
    isfinite(f) && isfinite(abs_Z) || return nothing

    scatter!(ax, [f], [abs_Z];
             color=_publication_colors()[2], markersize=7,
             strokecolor=:white, strokewidth=0.4)

    ax.xscale = log10
    _finalize_axis!(ax;
                    xlabel=rich("Frequency ", lsub("f", ""), " (Hz)"),
                    ylabel=rich("Impedance amplitude ", lsub("|Z|", ""), " (mΩ·cm²)"),
                    title="EIS - Bode amplitude",
                    legend=false)
    return nothing
end


"""Plot EIS Bode phase diagram with unified style conventions."""
function plot_EIS_curve_Bode_angle(cd::AbstractCurrent,
                                   Fourier_results::FourierOutputs,
                                   ax)
    _Z_real, _minus_Z_imag, f, _abs_Z, phi_deg = _eis_point(cd, Fourier_results)
    isfinite(f) && isfinite(phi_deg) || return nothing

    scatter!(ax, [f], [phi_deg];
             color=_publication_colors()[3], markersize=7,
             strokecolor=:white, strokewidth=0.4)

    ax.xscale = log10
    _finalize_axis!(ax;
                    xlabel=rich("Frequency ", lsub("f", ""), " (Hz)"),
                    ylabel=rich("Phase ", lsub("ϕ", ""), " (°)"),
                    title="EIS - Bode phase",
                    legend=false)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Final 2D Thermal Map
# ═══════════════════════════════════════════════════════════════════════════════

"""Plot final pseudo-2D temperature map.

The target figure is passed explicitly so that the colorbar layout does not rely
on `ax.parent`, which is more fragile when the figure grid changes.
"""
function plot_T_pseudo_2D_final(outputs::SimulationOutputs,
                                fc::AbstractFuelCell,
                                fig,
                                ax)
    temp_matrix = final_temperature_matrix_celsius(outputs)
    n_rows, n_cols = size(temp_matrix)
    nb_gdl = fc.numerical_parameters.nb_gdl
    nb_mpl = fc.numerical_parameters.nb_mpl

    # Orange very light (cold) -> red very dark (hot).
    thermal_cmap = [
        RGBf(0.99, 0.93, 0.82),  # Very light orange (cold)
        RGBf(0.99, 0.86, 0.62),
        RGBf(0.99, 0.76, 0.40),
        RGBf(0.98, 0.62, 0.20),
        RGBf(0.92, 0.42, 0.10),
        RGBf(0.78, 0.24, 0.08),
        RGBf(0.58, 0.12, 0.06),
        RGBf(0.40, 0.05, 0.05),  # Dark red (hot)
    ]

    finite_vals = [float(v) for v in temp_matrix if isfinite(v)]
    isempty(finite_vals) && (finite_vals = [0.0])

    # Automatic temperature limits from real data.
    tmin = minimum(finite_vals)
    tmax = maximum(finite_vals)
    tmax <= tmin && (tmax = tmin + 1e-6)

    # Ensure every cell gets a valid color even if some values are NaN/Inf.
    # `map` preserves matrix shape (unlike a plain comprehension over the matrix).
    plot_matrix = map(v -> isfinite(v) ? float(v) : float(tmin), temp_matrix)

    # Force orientation: MEA through-plane along x, GC direction along y.
    hm = heatmap!(ax, 1:n_cols, 1:n_rows, permutedims(plot_matrix);
                  colormap=thermal_cmap,
                  colorrange=(float(tmin), float(tmax)))

    ax.title = "Final temperature distribution"
    ax.xlabel = rich("Through cell thickness"; font=:bold)
    ax.ylabel = rich("Through gas-channel"; font=:bold)
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xminorgridvisible = false
    ax.yminorgridvisible = false
    xlims!(ax, 0.5, n_cols + 0.5)
    ylims!(ax, 0.5, n_rows + 0.5)
    ax.yreversed = true
    # Labels without numeric prefix, with indices in subscript (agdl_1 -> agdl₁ style).
    x_labels_rich = Any[
        rich("agc"),
        [rich("agdl", subscript("$(j)")) for j in 1:nb_gdl]...,
        [rich("ampl", subscript("$(j)")) for j in 1:nb_mpl]...,
        rich("acl"), rich("mem"), rich("ccl"),
        [rich("cmpl", subscript("$(j)")) for j in 1:nb_mpl]...,
        [rich("cgdl", subscript("$(j)")) for j in 1:nb_gdl]...,
        rich("cgc"),
    ]
    if length(x_labels_rich) != n_cols
        # Safe fallback in case matrix geometry and expected stack differ.
        x_labels_rich = [rich("x", subscript("$(i)")) for i in 1:n_cols]
    end
    ax.xticks = (1:length(x_labels_rich), x_labels_rich)
    ax.yticks = (1:n_rows, string.(1:n_rows))
    ax.xticklabelrotation = π / 3
    ax.xticklabelsize = 10
    ax.yticklabelsize = 11

    # Add a temperature scale to compare values across the 2D map.
    cbar_ticks = _colorbar_ticks_auto(tmin, tmax)
    Colorbar(fig[1, 2], hm;
             label=rich("Temperature ", lsub("T", ""), " (°C)"),
             ticks=cbar_ticks)


    T_des = fc.operating_conditions.T_des - 273.15
    T_des_round = round(T_des; digits=1)
    # White boxed annotation for readability, tightly fitted around the text.
    t_des_plain = "T_des = $(T_des_round) °C"
    pad_x = 0.012f0
    box_h = 0.060f0
    box_w = clamp(0.0096f0 * length(t_des_plain) + 2f0 * pad_x, 0.17f0, 0.27f0)
    x_right = 0.985f0
    x_left = x_right - box_w
    y_bottom = 0.008f0
    poly!(ax, Rect2f(x_left, y_bottom, box_w, box_h);
          color=(:white, 0.92), strokecolor=(:black, 0.75), strokewidth=0.8, space=:relative)
    text!(ax, x_left + pad_x, y_bottom + box_h / 2;
          text=rich("T", subscript("des"), " = $(T_des_round) °C"),
          align=(:left, :center), space=:relative, color=:black)
    return nothing
end

