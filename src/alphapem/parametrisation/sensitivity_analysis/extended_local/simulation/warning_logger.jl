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

"""
    Logging.handle_message(logger::_SimulationWarningLogger, level, message, _module, group, id, file, line; kwargs...)

Count selected simulation warnings and forward all other log messages.

Parameters
----------
logger : _SimulationWarningLogger
    Logger holding thread-safe warning counters and the wrapped logger.
level, message, _module, group, id, file, line
    Log record fields supplied by Julia's logging system.
kwargs
    Additional log metadata forwarded to the wrapped logger.

Returns
-------
Nothing
    Returns `nothing` for counted messages; otherwise returns the wrapped logger's result.
"""
Logging.handle_message(logger::_SimulationWarningLogger, level, message, _module, group, id, file, line; kwargs...) = begin
    # Match and count the known simulation safety-stop and timeout messages
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

    # Preserve normal logging for all messages not counted by this logger
    Logging.handle_message(logger.base_logger, level, message, _module, group, id, file, line; kwargs...)
end

"""
    Logging.shouldlog(logger::_SimulationWarningLogger, level, _module, group, id)

Allow the logger to receive every log level so it can count selected warnings.

Parameters
----------
logger : _SimulationWarningLogger
    Logger queried by the logging system.
level, _module, group, id
    Log record metadata supplied by Julia's logging system.

Returns
-------
Bool
    Always returns `true`.
"""
Logging.shouldlog(logger::_SimulationWarningLogger, level, _module, group, id) = true

"""
    Logging.min_enabled_level(logger::_SimulationWarningLogger)

Return the minimum enabled level of the wrapped logger.

Parameters
----------
logger : _SimulationWarningLogger
    Logger whose wrapped logger determines the minimum level.

Returns
-------
LogLevel
    Minimum level enabled by the wrapped logger.
"""
Logging.min_enabled_level(logger::_SimulationWarningLogger) = Logging.min_enabled_level(logger.base_logger)

"""
    Logging.catch_exceptions(logger::_SimulationWarningLogger)

Return whether the wrapped logger catches exceptions during log handling.

Parameters
----------
logger : _SimulationWarningLogger
    Logger whose wrapped logger determines exception handling.

Returns
-------
Bool
    Exception-handling policy of the wrapped logger.
"""
Logging.catch_exceptions(logger::_SimulationWarningLogger) = Logging.catch_exceptions(logger.base_logger)


"""
    Create fresh thread-safe counters for known simulation warnings.

    Parameters
    ----------
    None

    Returns
    -------
    Dict{Symbol, Threads.Atomic{Int}}
        Zeroed counters keyed by warning type.
"""
_new_warning_counters()::Dict{Symbol, Threads.Atomic{Int}} = Dict{Symbol, Threads.Atomic{Int}}(
    :oxygen_starvation => Threads.Atomic{Int}(0),
    :voltage_limit => Threads.Atomic{Int}(0),
    :max_runtime => Threads.Atomic{Int}(0),
)


"""
    _print_warning_summary(counters)

Print a summary of known warnings recorded during a simulation loop.

Parameters
----------
counters : Dict{Symbol, Threads.Atomic{Int}}
    Warning counters created by [`_new_warning_counters`](@ref).

Returns
-------
Nothing
    Returns `nothing` and emits no summary when all counters are zero.
"""
function _print_warning_summary(counters::Dict{Symbol, Threads.Atomic{Int}})
    # Skip output when no known warning occurred
    total = sum(c -> c[], values(counters))
    total == 0 && return nothing

    # Report the accumulated warning counts once after the simulation batch
    @info "Simulation warnings summary:" [
        "oxygen_starvation" => counters[:oxygen_starvation][],
        "voltage_limit" => counters[:voltage_limit][],
        "max_runtime" => counters[:max_runtime][],
    ]
    return nothing
end
