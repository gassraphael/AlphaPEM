# plot.jl
#
# Native CairoMakie plotting recipes for AlphaPEM.


using .PlotHelpers: _publication_colors,
                    _macroscopic_color,
                    _cell_current_color,
                    _cell_voltage_color,
                    _finalize_axis!,
                    _label_with_rmse,
                    _polarization_legend_base,
                    _experimental_marker,
                    _set_polarization_axis_limits!,
                    _set_polarization_fixed_ticks!,
                    _set_initial_temporal_limits!,
                    _format_fixed,
                    _nice_tick_step,
                    _colorbar_ticks_auto,
                    _set_dense_ticks!,
                    _plot_final_profile_along_gc!,
                    gc_direction_labels,
                    lsub

# ═══════════════════════════════════════════════════════════════════════════════
#  Internal States - Temporal Plots
# ═══════════════════════════════════════════════════════════════════════════════

"""Plot current density through the gas-channel direction as a function of time."""
function plot_ifc_1D_temporal(outputs::SimulationOutputs,
                              cd::AbstractCurrent,
                              cfg::SimulationConfig,
                              ax)
    palette = _publication_colors()
    t = time_history(outputs)
    nb_gc = cfg.numerical_parameters.nb_gc

    i_fc_cell_t = current(cd, t) ./ 1e4
    lines!(ax, t, i_fc_cell_t; color=_cell_current_color(nb_gc), linewidth=2.8, label=lsub("i", "fc,cell"))

    all_series = [i_fc_cell_t]
    if nb_gc > 1
        i_fc_series = [extract_derived_gc_series(outputs, i, x -> x.i_fc) ./ 1e4 for i in 1:nb_gc]
        for i in 1:nb_gc
            lines!(ax, t, i_fc_series[i]; color=palette[mod1(i, length(palette))], label=lsub("i", "fc,$(i)"))
        end
        append!(all_series, i_fc_series)
        _set_dense_ticks!(ax, t, all_series)
    else
        _set_dense_ticks!(ax, t, all_series)
    end

    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Current density ", rich(lsub("i", "fc"), "\n"), " (A·cm⁻²)"),
                    legend=true,
                    legend_position=:lb)
    _set_initial_temporal_limits!(ax, t, all_series, initial_time_range(outputs, cd, cfg))

    if nb_gc > 1
        # Add flow configuration annotation in bottom-right corner
        flow_note = cfg.type_flow == :counter_flow ?
            rich(lsub("i", "fc,1"), " = cathode inlet\n(counter-flow)") :
            rich(lsub("i", "fc,1"), " = inlet\n(co-flow)")
        pad_x = 0.012f0
        box_h = 0.130f0
        box_w = clamp(0.29f0, 0.22f0, 0.36f0)
        x_right = 0.975f0
        x_left = x_right - box_w
        y_bottom = 0.002f0
        poly!(ax, Rect2f(x_left, y_bottom, box_w, box_h);
              color=(:white, 0.92), strokecolor=(:black, 0.75), strokewidth=0.5, space=:relative)
        text!(ax, x_left + pad_x, y_bottom + box_h / 2;
              text=flow_note,
              align=(:left, :center), space=:relative, color=:black, fontsize=11)
    end

    return nothing
end

"""Plot vapour concentration through the MEA middle channel as a function of time."""
function plot_C_v_1D_temporal(outputs::SimulationOutputs,
                              cd::AbstractCurrent,
                              cfg::SimulationConfig,
                              ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(cfg)
    nb_mpl_mid = middle_mpl_index(cfg)
    t = time_history(outputs)
    T_ccl = extract_mid_mea_series(outputs, cfg, mea -> mea.ccl.T)
    C_v_sat_ccl = [C_v_sat(T) for T in T_ccl]

    series = [
        (extract_mid_mea_series(outputs, cfg, mea -> mea.agc.C_v), lsub("C", "v,agc")),
        (extract_mid_mea_series(outputs, cfg, mea -> mea.agdl[nb_gdl_mid].C_v), lsub("C", "v,agdl")),
        (extract_mid_mea_series(outputs, cfg, mea -> mea.ampl[nb_mpl_mid].C_v), lsub("C", "v,ampl")),
        (extract_mid_mea_series(outputs, cfg, mea -> mea.acl.C_v), lsub("C", "v,acl")),
        (extract_mid_mea_series(outputs, cfg, mea -> mea.ccl.C_v), lsub("C", "v,ccl")),
        (extract_mid_mea_series(outputs, cfg, mea -> mea.cmpl[nb_mpl_mid].C_v), lsub("C", "v,cmpl")),
        (extract_mid_mea_series(outputs, cfg, mea -> mea.cgdl[nb_gdl_mid].C_v), lsub("C", "v,cgdl")),
        (extract_mid_mea_series(outputs, cfg, mea -> mea.cgc.C_v), lsub("C", "v,cgc")),
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
    _set_initial_temporal_limits!(ax, t, vcat([s[1] for s in series], [C_v_sat_ccl]), initial_time_range(outputs, cd, cfg))
    return nothing
end

"""Plot membrane/CL water content as a function of time."""
function plot_lambda_1D_temporal(outputs::SimulationOutputs,
                                 cd::AbstractCurrent,
                                 cfg::SimulationConfig,
                                 ax)
    palette = _publication_colors()
    nb_gc = cfg.numerical_parameters.nb_gc
    mid_gc = middle_gas_channel_index(cfg)
    if cfg.display_timing == :postrun && cfg.type_current isa Union{PolarizationParams, PolarizationCalibrationParams}
        t_hist = time_history(outputs)
        sample_indices = polarisation_sampling_indices(outputs, cd)
        ifc_discretized = [current(cd, t_hist[idx]) / 1e4 for idx in sample_indices]

        lambda_acl_t = extract_mid_mea_series(outputs, cfg, mea -> mea.acl.lambda)[sample_indices]
        lambda_mem_t = extract_mid_mea_series(outputs, cfg, mea -> mea.mem.lambda)[sample_indices]
        lambda_ccl_t = extract_mid_mea_series(outputs, cfg, mea -> mea.ccl.lambda)[sample_indices]

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
        t = time_history(outputs)
        lines!(ax, t, extract_mid_mea_series(outputs, cfg, mea -> mea.acl.lambda);
               color=palette[3], label=lsub("λ", "acl"))
        lines!(ax, t, extract_mid_mea_series(outputs, cfg, mea -> mea.mem.lambda);
               color=palette[4], label=lsub("λ", "mem"))
        lambda_ccl_in_t = extract_mea_series(outputs, 1, mea -> mea.ccl.lambda)
        lambda_ccl_mid_t = extract_mea_series(outputs, mid_gc, mea -> mea.ccl.lambda)
        lambda_ccl_out_t = extract_mea_series(outputs, nb_gc, mea -> mea.ccl.lambda)
        lines!(ax, t, lambda_ccl_in_t;
               color=palette[5], linestyle=:dash, label=lsub("λ", "ccl,in"))
        mid_gc != 1 && mid_gc != nb_gc &&
            lines!(ax, t, lambda_ccl_mid_t;
                   color=palette[5], linestyle=:solid, label=lsub("λ", "ccl,mid"))
        nb_gc != 1 &&
            lines!(ax, t, lambda_ccl_out_t;
                   color=palette[5], linestyle=:dot, label=lsub("λ", "ccl,out"))
        lambda_acl_t = extract_mid_mea_series(outputs, cfg, mea -> mea.acl.lambda)
        lambda_mem_t = extract_mid_mea_series(outputs, cfg, mea -> mea.mem.lambda)
        y_series = [lambda_acl_t, lambda_mem_t, lambda_ccl_in_t]
        mid_gc != 1 && mid_gc != nb_gc && push!(y_series, lambda_ccl_mid_t)
        nb_gc != 1 && push!(y_series, lambda_ccl_out_t)
        _set_dense_ticks!(ax, t, y_series)
        _finalize_axis!(ax;
                        xlabel=rich("Time ", lsub("t", ""), " (s)"),
                        ylabel=rich("Water content ", lsub("λ", ""), " (-)"),
                        legend=true,
                        legend_position=:lb)
        _set_initial_temporal_limits!(ax, t, y_series, initial_time_range(outputs, cd, cfg))
    end
    return nothing
end

"""Plot liquid water saturation as a function of time."""
function plot_s_1D_temporal(outputs::SimulationOutputs,
                            cd::AbstractCurrent,
                            cfg::SimulationConfig,
                            ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(cfg)
    nb_mpl_mid = middle_mpl_index(cfg)
    if cfg.display_timing == :postrun && cfg.type_current isa Union{PolarizationParams, PolarizationCalibrationParams}
        t_hist = time_history(outputs)
        sample_indices = polarisation_sampling_indices(outputs, cd)
        ifc_discretized = [current(cd, t_hist[idx]) / 1e4 for idx in sample_indices]

        series = [
            (extract_mid_mea_series(outputs, cfg, mea -> mea.agc.s)[sample_indices], lsub("s", "agc")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.agdl[nb_gdl_mid].s)[sample_indices], lsub("s", "agdl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.ampl[nb_mpl_mid].s)[sample_indices], lsub("s", "ampl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.acl.s)[sample_indices], lsub("s", "acl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.ccl.s)[sample_indices], lsub("s", "ccl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cmpl[nb_mpl_mid].s)[sample_indices], lsub("s", "cmpl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cgdl[nb_gdl_mid].s)[sample_indices], lsub("s", "cgdl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cgc.s)[sample_indices], lsub("s", "cgc")),
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
        t = time_history(outputs)

        series = [
            (extract_mid_mea_series(outputs, cfg, mea -> mea.agc.s), lsub("s", "agc")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.agdl[nb_gdl_mid].s), lsub("s", "agdl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.ampl[nb_mpl_mid].s), lsub("s", "ampl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.acl.s), lsub("s", "acl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.ccl.s), lsub("s", "ccl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cmpl[nb_mpl_mid].s), lsub("s", "cmpl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cgdl[nb_gdl_mid].s), lsub("s", "cgdl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cgc.s), lsub("s", "cgc")),
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
        _set_initial_temporal_limits!(ax, t, [s[1] for s in series], initial_time_range(outputs, cd, cfg))
    end
    return nothing
end

"""Plot hydrogen concentration as a function of time."""
function plot_C_H2_1D_temporal(outputs::SimulationOutputs,
                               cd::AbstractCurrent,
                               cfg::SimulationConfig,
                               ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(cfg)
    nb_mpl_mid = middle_mpl_index(cfg)
    t = time_history(outputs)

    lines!(ax, t, extract_mid_mea_series(outputs, cfg, mea -> mea.agc.C_H2);
           color=palette[1], label=lsub("C", "H₂,agc"))
    lines!(ax, t, extract_mid_mea_series(outputs, cfg, mea -> mea.agdl[nb_gdl_mid].C_H2);
           color=palette[2], label=lsub("C", "H₂,agdl"))
    lines!(ax, t, extract_mid_mea_series(outputs, cfg, mea -> mea.ampl[nb_mpl_mid].C_H2);
           color=palette[3], label=lsub("C", "H₂,ampl"))
    lines!(ax, t, extract_mid_mea_series(outputs, cfg, mea -> mea.acl.C_H2);
           color=palette[4], label=lsub("C", "H₂,acl"))

    C_H2_agc_t = extract_mid_mea_series(outputs, cfg, mea -> mea.agc.C_H2)
    C_H2_agdl_t = extract_mid_mea_series(outputs, cfg, mea -> mea.agdl[nb_gdl_mid].C_H2)
    C_H2_ampl_t = extract_mid_mea_series(outputs, cfg, mea -> mea.ampl[nb_mpl_mid].C_H2)
    C_H2_acl_t = extract_mid_mea_series(outputs, cfg, mea -> mea.acl.C_H2)
    _set_dense_ticks!(ax, t, [C_H2_agc_t, C_H2_agdl_t, C_H2_ampl_t, C_H2_acl_t])
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Hydrogen concentration ", lsub("C", "H₂"), " (mol·m⁻³)"),
                    legend=true,
                    legend_position=:rt)
    _set_initial_temporal_limits!(ax, t, [C_H2_agc_t, C_H2_agdl_t, C_H2_ampl_t, C_H2_acl_t], initial_time_range(outputs, cd, cfg))
    return nothing
end

"""Plot oxygen concentration as a function of time."""
function plot_C_O2_1D_temporal(outputs::SimulationOutputs,
                               cd::AbstractCurrent,
                               cfg::SimulationConfig,
                               ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(cfg)
    nb_mpl_mid = middle_mpl_index(cfg)
    nb_gc = cfg.numerical_parameters.nb_gc
    mid_gc = middle_gas_channel_index(cfg)
    t = time_history(outputs)

    C_O2_Pt_t = extract_mid_derived_gc_series(outputs, cfg, x -> x.C_O2_Pt)
    lines!(ax, t, C_O2_Pt_t;
           color=palette[10], label=lsub("C", "O₂,Pt"))
    lines!(ax, t, extract_mid_mea_series(outputs, cfg, mea -> mea.ccl.C_O2);
           color=palette[5], label=lsub("C", "O₂,ccl"))
    lines!(ax, t, extract_mid_mea_series(outputs, cfg, mea -> mea.cmpl[nb_mpl_mid].C_O2);
           color=palette[6], label=lsub("C", "O₂,cmpl"))
    lines!(ax, t, extract_mid_mea_series(outputs, cfg, mea -> mea.cgdl[nb_gdl_mid].C_O2);
           color=palette[7], label=lsub("C", "O₂,cgdl"))
    C_O2_cgc_in_t = extract_mea_series(outputs, 1, mea -> mea.cgc.C_O2)
    C_O2_cgc_mid_t = extract_mea_series(outputs, mid_gc, mea -> mea.cgc.C_O2)
    C_O2_cgc_out_t = extract_mea_series(outputs, nb_gc, mea -> mea.cgc.C_O2)
    lines!(ax, t, C_O2_cgc_in_t;
           color=palette[8], linestyle=:dash, label=lsub("C", "O₂,cgc,in"))
    mid_gc != 1 && mid_gc != nb_gc &&
        lines!(ax, t, C_O2_cgc_mid_t;
               color=palette[8], linestyle=:solid, label=lsub("C", "O₂,cgc,mid"))
    nb_gc != 1 &&
        lines!(ax, t, C_O2_cgc_out_t;
               color=palette[8], linestyle=:dot, label=lsub("C", "O₂,cgc,out"))

    C_O2_ccl_t = extract_mid_mea_series(outputs, cfg, mea -> mea.ccl.C_O2)
    C_O2_cmpl_t = extract_mid_mea_series(outputs, cfg, mea -> mea.cmpl[nb_mpl_mid].C_O2)
    C_O2_cgdl_t = extract_mid_mea_series(outputs, cfg, mea -> mea.cgdl[nb_gdl_mid].C_O2)
    y_series = [C_O2_Pt_t, C_O2_ccl_t, C_O2_cmpl_t, C_O2_cgdl_t, C_O2_cgc_in_t]
    mid_gc != 1 && mid_gc != nb_gc && push!(y_series, C_O2_cgc_mid_t)
    nb_gc != 1 && push!(y_series, C_O2_cgc_out_t)
    _set_dense_ticks!(ax, t, y_series)
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Oxygen concentration ", lsub("C", "O₂"), " (mol·m⁻³)"),
                    legend=true,
                    legend_position=:lt)
    _set_initial_temporal_limits!(ax, t, y_series, initial_time_range(outputs, cd, cfg))
    return nothing
end

"""Plot temperature through the MEA as a function of time."""
function plot_T_1D_temporal(outputs::SimulationOutputs,
                            fc::AbstractFuelCell,
                            cd::AbstractCurrent,
                            cfg::SimulationConfig,
                            ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(cfg)
    nb_mpl_mid = middle_mpl_index(cfg)
    T_des = fc.operating_conditions.T_des - 273.15
    if cfg.display_timing == :postrun && cfg.type_current isa Union{PolarizationParams, PolarizationCalibrationParams}
        t_hist = time_history(outputs)
        sample_indices = polarisation_sampling_indices(outputs, cd)
        ifc_discretized = [current(cd, t_hist[idx]) / 1e4 for idx in sample_indices]

        series = [
            (extract_mid_mea_series(outputs, cfg, mea -> mea.agc.T)[sample_indices] .- 273.15, lsub("T", "agc")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.agdl[nb_gdl_mid].T)[sample_indices] .- 273.15, lsub("T", "agdl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.ampl[nb_mpl_mid].T)[sample_indices] .- 273.15, lsub("T", "ampl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.acl.T)[sample_indices] .- 273.15, lsub("T", "acl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.mem.T)[sample_indices] .- 273.15, lsub("T", "mem")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.ccl.T)[sample_indices] .- 273.15, lsub("T", "ccl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cmpl[nb_mpl_mid].T)[sample_indices] .- 273.15, lsub("T", "cmpl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cgdl[nb_gdl_mid].T)[sample_indices] .- 273.15, lsub("T", "cgdl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cgc.T)[sample_indices] .- 273.15, lsub("T", "cgc")),
        ]

        for (i, (y, lbl)) in enumerate(series)
            lines!(ax, ifc_discretized, y; color=palette[mod1(i, length(palette))], linewidth=2.4, label=lbl)
        end

        _set_dense_ticks!(ax, ifc_discretized, [s[1] for s in series])
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Temperature ", lsub("T", ""), " (°C)"),
                        legend=true,
                        legend_position=:lt)
    else
        t = time_history(outputs)

        series = [
            (extract_mid_mea_series(outputs, cfg, mea -> mea.agc.T) .- 273.15, lsub("T", "agc")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.agdl[nb_gdl_mid].T) .- 273.15, lsub("T", "agdl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.ampl[nb_mpl_mid].T) .- 273.15, lsub("T", "ampl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.acl.T) .- 273.15, lsub("T", "acl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.mem.T) .- 273.15, lsub("T", "mem")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.ccl.T) .- 273.15, lsub("T", "ccl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cmpl[nb_mpl_mid].T) .- 273.15, lsub("T", "cmpl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cgdl[nb_gdl_mid].T) .- 273.15, lsub("T", "cgdl")),
            (extract_mid_mea_series(outputs, cfg, mea -> mea.cgc.T) .- 273.15, lsub("T", "cgc")),
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
        _set_initial_temporal_limits!(ax, t, [s[1] for s in series], initial_time_range(outputs, cd, cfg))
    end
    return nothing
end

"""Plot cell voltage as a function of time."""
function plot_Ucell(outputs::SimulationOutputs,
                    cd::AbstractCurrent,
                    cfg::SimulationConfig,
                    ax)
    t = time_history(outputs)
    Ucell_t = extract_derived_series(outputs, x -> x.Ucell)
    lines!(ax, t, Ucell_t; color=_cell_voltage_color(), label=lsub("U", "cell"))
    _set_dense_ticks!(ax, t, [Ucell_t])
    _finalize_axis!(ax;
                    xlabel=rich("Time ", lsub("t", ""), " (s)"),
                    ylabel=rich("Cell voltage ", lsub("U", "cell"), " (V)"),
                    legend=true,
                    legend_position=:lt)
    _set_initial_temporal_limits!(ax, t, [Ucell_t], initial_time_range(outputs, cd, cfg))
    return nothing
end

"""Plot nitrogen concentration at AGC and CGC for each GC node as a function of time."""
function plot_C_N2_1D_temporal(outputs::SimulationOutputs,
                               cd::AbstractCurrent,
                               cfg::SimulationConfig,
                               ax)
    palette = _publication_colors()
    nb_gc = cfg.numerical_parameters.nb_gc
    t = time_history(outputs)

    C_N2_series = Vector{Vector{Float64}}()
    for i in 1:nb_gc
        C_N2_agc_i = extract_mea_series(outputs, i, mea -> mea.agc.C_N2)
        C_N2_cgc_i = extract_mea_series(outputs, i, mea -> mea.cgc.C_N2)
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
    _set_initial_temporal_limits!(ax, t, C_N2_series, initial_time_range(outputs, cd, cfg))
    return nothing
end

"""Plot relative humidity through the anode layers at the middle GC node as a function of time."""
function plot_Phi_a_1D_temporal(outputs::SimulationOutputs,
                                cd::AbstractCurrent,
                                cfg::SimulationConfig,
                                ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(cfg)
    nb_mpl_mid = middle_mpl_index(cfg)
    t = time_history(outputs)

    series = [
        (extract_mid_mea_series(outputs, cfg,
             mea -> mea.agc.C_v / C_v_sat(mea.agc.T)),                                    lsub("Φ", "a,agc")),
        (extract_mid_mea_series(outputs, cfg,
             mea -> mea.agdl[nb_gdl_mid].C_v / C_v_sat(mea.agdl[nb_gdl_mid].T)),          lsub("Φ", "a,agdl")),
        (extract_mid_mea_series(outputs, cfg,
             mea -> mea.ampl[nb_mpl_mid].C_v / C_v_sat(mea.ampl[nb_mpl_mid].T)),          lsub("Φ", "a,ampl")),
        (extract_mid_mea_series(outputs, cfg,
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
    _set_initial_temporal_limits!(ax, t, vcat([s[1] for s in series], [ones(length(t))]), initial_time_range(outputs, cd, cfg))
    return nothing
end

"""Plot relative humidity through the cathode layers at the middle GC node as a function of time."""
function plot_Phi_c_1D_temporal(outputs::SimulationOutputs,
                                cd::AbstractCurrent,
                                cfg::SimulationConfig,
                                ax)
    palette = _publication_colors()
    nb_gdl_mid = middle_gdl_index(cfg)
    nb_mpl_mid = middle_mpl_index(cfg)
    t = time_history(outputs)

    series = [
        (extract_mid_mea_series(outputs, cfg,
             mea -> mea.ccl.C_v / C_v_sat(mea.ccl.T)),                                    lsub("Φ", "c,ccl")),
        (extract_mid_mea_series(outputs, cfg,
             mea -> mea.cmpl[nb_mpl_mid].C_v / C_v_sat(mea.cmpl[nb_mpl_mid].T)),          lsub("Φ", "c,cmpl")),
        (extract_mid_mea_series(outputs, cfg,
             mea -> mea.cgdl[nb_gdl_mid].C_v / C_v_sat(mea.cgdl[nb_gdl_mid].T)),          lsub("Φ", "c,cgdl")),
        (extract_mid_mea_series(outputs, cfg,
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
    _set_initial_temporal_limits!(ax, t, vcat([s[1] for s in series], [ones(length(t))]), initial_time_range(outputs, cd, cfg))
    return nothing
end

"""Plot anode and cathode gas velocities for each GC node as a function of time."""
function plot_v_1D_temporal(outputs::SimulationOutputs,
                            cd::AbstractCurrent,
                            cfg::SimulationConfig,
                            ax)
    palette = _publication_colors()
    nb_gc = cfg.numerical_parameters.nb_gc
    t = time_history(outputs)

    v_series = Vector{Vector{Float64}}()
    for i in 1:nb_gc
        v_a_i = extract_derived_gc_series(outputs, i, x -> x.v_a)
        v_c_i = extract_derived_gc_series(outputs, i, x -> x.v_c)
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
    _set_initial_temporal_limits!(ax, t, v_series, initial_time_range(outputs, cd, cfg))
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
    nb_gc = cfg.numerical_parameters.nb_gc
    t = time_history(outputs)

    # Use the centralized calculation function
    Re_a_full, Re_c_full = calculate_reynolds_numbers(outputs, fc)

    Re_series = Vector{Vector{Float64}}()
    for i in 1:nb_gc
        Re_a_i = Re_a_full[i]
        Re_c_i = Re_c_full[i]

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
    _set_initial_temporal_limits!(ax, t, Re_series, initial_time_range(outputs, cd, cfg))
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
    t = time_history(outputs)

    C_v_agc = extract_mid_mea_series(outputs, cfg, mea -> mea.agc.C_v)
    C_H2_agc = extract_mid_mea_series(outputs, cfg, mea -> mea.agc.C_H2)
    C_N2_agc = extract_mid_mea_series(outputs, cfg, mea -> mea.agc.C_N2)
    T_agc = extract_mid_mea_series(outputs, cfg, mea -> mea.agc.T)
    C_v_cgc = extract_mid_mea_series(outputs, cfg, mea -> mea.cgc.C_v)
    C_O2_cgc = extract_mid_mea_series(outputs, cfg, mea -> mea.cgc.C_O2)
    C_N2_cgc = extract_mid_mea_series(outputs, cfg, mea -> mea.cgc.C_N2)
    T_cgc = extract_mid_mea_series(outputs, cfg, mea -> mea.cgc.T)

    P_agc = (C_v_agc .+ C_H2_agc .+ C_N2_agc) .* R .* T_agc ./ 1e5
    P_cgc = (C_v_cgc .+ C_O2_cgc .+ C_N2_cgc) .* R .* T_cgc ./ 1e5
    Pa_in = extract_derived_series(outputs, x -> x.Pa_in) ./ 1e5
    Pc_in = extract_derived_series(outputs, x -> x.Pc_in) ./ 1e5
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
    _set_initial_temporal_limits!(ax, t, [P_agc, P_cgc, Pa_in, Pa_out, Pc_in, Pc_out], initial_time_range(outputs, cd, cfg))
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
    model_color = _cell_voltage_color()
    exp_color = RGBf(0.10, 0.10, 0.10)
    model_label_base = _polarization_legend_base(cfg.type_fuel_cell;
                                                 simulation=true,
                                                 calibration=false)
    exp_label = _polarization_legend_base(cfg.type_fuel_cell;
                                          simulation=false,
                                          calibration=false)

    if cfg.display_timing == :postrun
        ifc_discretized, Ucell_discretized = _polarization_points(outputs, cd; average=true)

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
        sc = scatter!(ax, [ifc_last], [Ucell_last];
                 color=:white, marker=:circle, markersize=8,
                 strokecolor=model_color, strokewidth=1.8)

        _set_polarization_axis_limits!(ax)
        _set_polarization_fixed_ticks!(ax)
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Cell voltage ", lsub("U", "cell"), " (V)"),
                        title="Polarization curve",
                        legend=false)   # skip auto-legend (would create N entries)
        # Single legend entry representing the model curve style.
        axislegend(ax, [sc], [model_label_base];
                   position=:rt, framevisible=true,
                   framewidth=0.8, framecolor=(:black, 0.35),
                   backgroundcolor=(:white, 0.88), padding=(6, 6, 6, 6))
    end
    return nothing
end


"""Plot polarization hysteresis curve (no averaging)."""
function plot_polarization_hysteresis(outputs::SimulationOutputs,
                                      cd::AbstractCurrent,
                                      cfg::SimulationConfig,
                                      ax)
    model_color = _cell_voltage_color()
    model_label_base = _polarization_legend_base(cfg.type_fuel_cell;
                                                 simulation=true,
                                                 calibration=false)

    if cfg.display_timing == :postrun
        ifc_full, Ucell_full = _polarization_points(outputs, cd; average=false)

        lines!(ax, ifc_full, Ucell_full;
               color=model_color, linewidth=2.5, label=model_label_base)

        _set_polarization_axis_limits!(ax)
        _set_polarization_fixed_ticks!(ax)
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Cell voltage ", lsub("U", "cell"), " (V)"),
                        title="Polarization curve with hysteresis",
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
    model_color = _cell_voltage_color()
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
        sc = scatter!(ax, [ifc_last], [Ucell_last];
                 color=:white, marker=:circle, markersize=8,
                 strokecolor=model_color, strokewidth=1.8)
        _set_polarization_axis_limits!(ax)
        _set_polarization_fixed_ticks!(ax)
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Cell voltage ", lsub("U", "cell"), " (V)"),
                        title="Polarization point (calibration, dynamic)",
                        legend=false)   # skip auto-legend
        axislegend(ax, [sc], [model_label_base];
                   position=:rt, framevisible=true,
                   framewidth=0.8, framecolor=(:black, 0.35),
                   backgroundcolor=(:white, 0.88), padding=(6, 6, 6, 6))
    end
    return nothing
end



# ═══════════════════════════════════════════════════════════════════════════════
#  Derived Polarization Curves
# ═══════════════════════════════════════════════════════════════════════════════

"""Plot power density curve with unified CairoMakie style conventions.

Power density is defined as P = U_cell × i_fc, expressed in W·cm⁻²."""
function plot_power_density_curve(outputs::SimulationOutputs,
                                  cd::AbstractCurrent,
                                  cfg::SimulationConfig,
                                  ax)
    model_color = _macroscopic_color()
    model_label = "Model ($(String(cfg.type_fuel_cell)))"

    if cfg.display_timing == :postrun
        ifc_discretized, Ucell_discretized = _polarization_points(outputs, cd)
        P_discretized = ifc_discretized .* Ucell_discretized   # W·cm⁻²
        lines!(ax, ifc_discretized, P_discretized;
               color=model_color, linewidth=3.0, label=model_label)
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
        scatter!(ax, [ifc_last], [P_last]; color=model_color, markersize=8, label=model_label)
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
                              cd::AbstractCurrent,
                              cfg::SimulationConfig,
                              ax)
    E_th = 1.481   # V — HHV-based thermoneutral voltage for H₂/O₂ PEM fuel cell.
    model_color = _macroscopic_color()
    model_label = "Model ($(String(cfg.type_fuel_cell)))"

    if cfg.display_timing == :postrun
        ifc_discretized, Ucell_discretized = _polarization_points(outputs, cd)
        eta_discretized = Ucell_discretized ./ E_th
        lines!(ax, ifc_discretized, eta_discretized;
               color=model_color, linewidth=3.0, label=model_label)
        _set_dense_ticks!(ax, ifc_discretized, [eta_discretized])
        # Ensure efficiency curve is properly framed around its values.
        ylims!(ax, minimum(eta_discretized) * 0.95, min(1.0, maximum(eta_discretized) * 1.05))
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Voltage efficiency ", lsub("η", "v"), " (–)"),
                        title="Cell efficiency curve",
                        legend=true,
                        legend_position=:rt)
    else
        t_hist   = time_history(outputs)
        Ucell_t  = derived_outputs(outputs).Ucell
        idx      = lastindex(t_hist)
        ifc_last = current(cd, t_hist[idx]) / 1e4
        eta_last = Ucell_t[idx] / E_th
        scatter!(ax, [ifc_last], [eta_last]; color=model_color, markersize=8, label=model_label)
        _set_dense_ticks!(ax, [ifc_last], [[eta_last]])
        _finalize_axis!(ax;
                        xlabel=rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                        ylabel=rich("Voltage efficiency ", lsub("η", "v"), " (–)"),
                        title="Cell efficiency point (dynamic)",
                        legend=true,
                        legend_position=:rt)
    end
    return nothing
end


# ═══════════════════════════════════════════════════════════════════════════════
#  EIS Curves
# ═══════════════════════════════════════════════════════════════════════════════

"""Plot EIS Nyquist diagram with unified style conventions."""
function plot_EIS_curve_Nyquist(cd::AbstractCurrent,
                                Fourier_results::AbstractVector{FourierOutputs},
                                ax)
    
    Z_reals = Float64[]
    minus_Z_imags = Float64[]
    
    for res in Fourier_results
        Z_real, minus_Z_imag, _, _, _ = _eis_point(cd, res)
        if isfinite(Z_real) && isfinite(minus_Z_imag)
            push!(Z_reals, Z_real)
            push!(minus_Z_imags, minus_Z_imag)
        end
    end
    
    isempty(Z_reals) && return nothing

    scatter!(ax, Z_reals, minus_Z_imags;
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
                                       Fourier_results::AbstractVector{FourierOutputs},
                                       ax)
    
    fs = Float64[]
    abs_Zs = Float64[]
    
    for res in Fourier_results
        _, _, f, abs_Z, _ = _eis_point(cd, res)
        if isfinite(f) && isfinite(abs_Z)
            push!(fs, f)
            push!(abs_Zs, abs_Z)
        end
    end
    
    isempty(fs) && return nothing

    scatter!(ax, fs, abs_Zs;
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
                                   Fourier_results::AbstractVector{FourierOutputs},
                                   ax)
    
    fs = Float64[]
    phi_degs = Float64[]
    
    for res in Fourier_results
        _, _, f, _, phi_deg = _eis_point(cd, res)
        if isfinite(f) && isfinite(phi_deg)
            push!(fs, f)
            push!(phi_degs, phi_deg)
        end
    end
    
    isempty(fs) && return nothing

    scatter!(ax, fs, phi_degs;
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
#  Final GC Maps and Profiles
# ═══════════════════════════════════════════════════════════════════════════════



"""Plot the final current-density distribution along the gas channel."""
function plot_ifc_GC_final(outputs::SimulationOutputs,
                           cd::AbstractCurrent,
                           cfg::SimulationConfig,
                           ax)
    palette = _publication_colors()
    nb_gc = cfg.numerical_parameters.nb_gc
    x_gc = collect(1:nb_gc)

    i_fc_hist = derived_outputs(outputs).i_fc
    i_fc_final = [i_fc_hist[k][end] for k in 1:nb_gc] ./ 1e4 # A·cm⁻²
    i_fc_cell_final = current(cd, time_history(outputs)[end]) / 1e4

    _plot_final_profile_along_gc!(ax, x_gc, i_fc_final, cfg;
                                  color=palette[1],
                                  series_label=rich(lsub("i", "fc,1D,final")),
                                  ylabel=rich("Final current density ", lsub("i", "fc"), " (A·cm⁻²)"),
                                  reference_value=i_fc_cell_final,
                                  reference_label=rich(lsub("i", "fc,cell,final")),
                                  legend_position=:rt)
    return nothing
end

"""Plot the final dissolved-oxygen concentration at Pt sites along the gas channel."""
function plot_C_O2_Pt_GC_final(outputs::SimulationOutputs,
                               cfg::SimulationConfig,
                               ax)
    palette = _publication_colors()
    nb_gc = cfg.numerical_parameters.nb_gc
    x_gc = collect(1:nb_gc)

    C_O2_Pt_hist = derived_outputs(outputs).C_O2_Pt
    C_O2_Pt_final = [C_O2_Pt_hist[k][end] for k in 1:nb_gc]

    _plot_final_profile_along_gc!(ax, x_gc, C_O2_Pt_final, cfg;
                                  color=palette[10],
                                  series_label=rich(lsub("C", "O₂,Pt,final")),
                                  ylabel=rich("Final oxygen concentration ", lsub("C", "O₂,Pt"), " (mol·m⁻³)"),
                                  legend_position=:rt)
    return nothing
end

"""Plot the final membrane water-content distribution along the gas channel."""
function plot_lambda_mem_GC_final(outputs::SimulationOutputs,
                                  cfg::SimulationConfig,
                                  ax)
    palette = _publication_colors()
    nb_gc = cfg.numerical_parameters.nb_gc
    x_gc = collect(1:nb_gc)

    last_state = solver_state_history(outputs)[end]
    lambda_mem_final = [mea_state_at(last_state, k).mem.lambda for k in 1:nb_gc]

    _plot_final_profile_along_gc!(ax, x_gc, lambda_mem_final, cfg;
                                  color=palette[4],
                                  series_label=rich(lsub("λ", "mem,final")),
                                  ylabel=rich("Final water content ", lsub("λ", "mem"), " (-)"),
                                  legend_position=:rt)
    return nothing
end

"""Plot final pseudo-2D temperature map.

The target figure is passed explicitly so that the colorbar layout does not rely
on `ax.parent`, which is more fragile when the figure grid changes.
"""
function plot_T_pseudo_2D_final(outputs::SimulationOutputs,
                                fc::AbstractFuelCell,
                                fig,
                                ax,
                                cfg::SimulationConfig)
    temp_matrix = final_temperature_matrix_celsius(outputs)
    n_rows, n_cols = size(temp_matrix)
    nb_gdl = cfg.numerical_parameters.nb_gdl
    nb_mpl = cfg.numerical_parameters.nb_mpl

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

    # Build Y-axis labels: annotate the two extremities while keeping numeric labels inside.
    y_labels_rich = gc_direction_labels(cfg, n_rows)

    ax.xticks = (1:length(x_labels_rich), x_labels_rich)
    ax.yticks = (1:n_rows, y_labels_rich)
    ax.xticklabelrotation = π / 3
    ax.xticklabelsize = 10
    ax.yticklabelsize = 10

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
     box_w = clamp(0.0096f0 * length(t_des_plain) + 2f0 * pad_x, 0.20f0, 0.32f0)
     x_right = 0.985f0
     x_left = x_right - box_w
     y_bottom_text = 0.008f0
     poly!(ax, Rect2f(x_left, y_bottom_text, box_w, box_h);
           color=(:white, 0.92), strokecolor=(:black, 0.75), strokewidth=0.8, space=:relative)
     text!(ax, x_left + pad_x, y_bottom_text + box_h / 2;
           text=rich("T", subscript("des"), " = $(T_des_round) °C"),
           align=(:left, :center), space=:relative, color=:black)
    return nothing
end

