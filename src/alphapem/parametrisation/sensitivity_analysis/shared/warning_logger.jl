# -*- coding: utf-8 -*-

"""
    _SimulationWarningLogger

Custom logger that counts known AlphaPEM safety-stop and timeout messages during
batch simulations instead of printing them to the console.

Matched messages are silently counted; all other messages are forwarded to the
base logger. The counters are thread-safe so they can be incremented from
parallel simulation tasks.
"""
struct _SimulationWarningLogger <: AbstractLogger
    counters::Dict{Symbol, Threads.Atomic{Int}}
    base_logger::AbstractLogger
end

Logging.handle_message(logger::_SimulationWarningLogger, level, message, _module, group, id, file, line; kwargs...) = begin
    msg = string(message)
    if level == Logging.Warn && occursin("Safety stop: oxygen concentration", msg)
        Threads.atomic_add!(logger.counters[:oxygen_starvation], 1)
        return nothing
    elseif level == Logging.Warn && occursin("Safety stop: cell voltage", msg)
        Threads.atomic_add!(logger.counters[:voltage_limit], 1)
        return nothing
    elseif level == Logging.Error && occursin("Simulation exceeded maximum runtime", msg)
        Threads.atomic_add!(logger.counters[:max_runtime], 1)
        return nothing
    end
    Logging.handle_message(logger.base_logger, level, message, _module, group, id, file, line; kwargs...)
end

Logging.shouldlog(logger::_SimulationWarningLogger, level, _module, group, id) = true
Logging.min_enabled_level(logger::_SimulationWarningLogger) = Logging.min_enabled_level(logger.base_logger)
Logging.catch_exceptions(logger::_SimulationWarningLogger) = Logging.catch_exceptions(logger.base_logger)


"""
    _new_warning_counters()

Return a fresh counter dictionary with the keys expected by `_SimulationWarningLogger`.
"""
_new_warning_counters()::Dict{Symbol, Threads.Atomic{Int}} = Dict{Symbol, Threads.Atomic{Int}}(
    :oxygen_starvation => Threads.Atomic{Int}(0),
    :voltage_limit => Threads.Atomic{Int}(0),
    :max_runtime => Threads.Atomic{Int}(0),
)


"""
    _print_warning_summary(counters)

Print a single summary line for each warning type that occurred during the
simulation loop.
"""
function _print_warning_summary(counters::Dict{Symbol, Threads.Atomic{Int}})
    total = sum(c -> c[], values(counters))
    total == 0 && return nothing
    @info "Simulation warnings summary:" [
        "oxygen_starvation" => counters[:oxygen_starvation][],
        "voltage_limit" => counters[:voltage_limit][],
        "max_runtime" => counters[:max_runtime][],
    ]
    return nothing
end
