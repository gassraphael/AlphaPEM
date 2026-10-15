# -*- coding: utf-8 -*-

"""
    _save_jld2_atomic(path, label; data...)

Atomically persist `data` to the JLD2 file `path`.

The data is first written to a temporary file and then renamed so that an
interrupted write never leaves a corrupt checkpoint behind. Failures only emit a
warning (`label` names the analysis) so that a checkpoint problem never aborts a run.
"""
function _save_jld2_atomic(path::String, label::String; data...)
    mkpath(dirname(path))
    tmp_path = path * ".tmp"
    try
        JLD2.jldsave(tmp_path; data...)
        mv(tmp_path, path; force = true)
    catch e
        @warn "Failed to save $label checkpoint: $e"
    end
    return nothing
end


"""
    _load_jld2_checked(path, expected_hash, label)

Load a JLD2 checkpoint and return its content as a `Dict`.

Returns `nothing` if the file does not exist or if its stored `config_hash`
differs from `expected_hash` (a warning is emitted in the latter case).
"""
function _load_jld2_checked(path::String, expected_hash::String, label::String)
    isfile(path) || return nothing
    data = JLD2.load(path)
    if data["config_hash"] != expected_hash
        @warn "$label checkpoint config hash mismatch. Ignoring checkpoint and starting fresh."
        return nothing
    end
    return data
end


"""
    _find_missing_sample_ids(df_curves, n_total)

Return the list of `sample_id` values in `1:n_total` that have not yet been completed.

A sample is considered complete as soon as it has a row in `df_curves`
(i.e. the simulation was attempted). Missing rows are those that must still be run.
"""
function _find_missing_sample_ids(df_curves::DataFrame, n_total::Int)::Vector{Int}
    completed = Set{Int}(df_curves.sample_id)
    return Int[i for i in 1:n_total if !(i in completed)]
end

_find_missing_sample_ids(::Nothing, n_total::Int)::Vector{Int} = collect(1:n_total)


"""
    _remove_checkpoint(path)

Delete the checkpoint file at `path` if it exists. Called after a successful run.
"""
function _remove_checkpoint(path::String)
    isfile(path) && rm(path)
    return nothing
end
