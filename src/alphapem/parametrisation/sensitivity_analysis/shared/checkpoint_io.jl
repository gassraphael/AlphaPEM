# -*- coding: utf-8 -*-

"""
    _save_jld2_atomic(path, label; data...)

Atomically persist named data to a JLD2 checkpoint.

The data is first written to a temporary file and then renamed so that an
interrupted write never leaves a corrupt checkpoint behind. Failures only emit a
warning (`label` names the analysis) so that a checkpoint problem never aborts a run.

Parameters
----------
path : String
    Destination checkpoint path.
label : String
    Analysis label included in warning messages.
data : NamedTuple
    Named values to write to the checkpoint.

Returns
-------
Nothing
    Always returns `nothing`; persistence failures are reported with a warning.
"""
function _save_jld2_atomic(path::String, label::String; data...)
    # Prepare a temporary file beside the destination for an atomic replacement
    mkpath(dirname(path))
    tmp_path = path * ".tmp"

    # Write all checkpoint data before replacing the destination
    try
        JLD2.jldsave(tmp_path; data...)
        mv(tmp_path, path; force = true)
    catch e
        # A failed checkpoint must not interrupt the analysis run.
        @warn "Failed to save $label checkpoint: $e"
    end
    return nothing
end


"""
    _load_jld2_checked(path, expected_hash, label)

Load and validate a JLD2 checkpoint.

Parameters
----------
path : String
    Checkpoint path to read.
expected_hash : String
    Configuration hash expected in the checkpoint.
label : String
    Analysis label included in a mismatch warning.

Returns
-------
Dict or nothing
    Checkpoint contents, or `nothing` if the file is absent or its configuration hash differs.
"""
function _load_jld2_checked(path::String, expected_hash::String, label::String)
    # Skip loading when there is no checkpoint to resume
    isfile(path) || return nothing
    data = JLD2.load(path)

    # Reject results produced with a different configuration
    if data["config_hash"] != expected_hash
        @warn "$label checkpoint config hash mismatch. Ignoring checkpoint and starting fresh."
        return nothing
    end

    # Return only a checkpoint matching the current analysis configuration
    return data
end


"""
    _find_missing_sample_ids(df_curves, n_total)

Find sample IDs that have not yet been attempted.

A sample is considered complete as soon as it has a row in `df_curves`.

Parameters
----------
df_curves : DataFrame or nothing
    Existing curve results; `nothing` means no samples have been attempted.
n_total : Int
    Total number of sample IDs.

Returns
-------
Vector{Int}
    IDs in `1:n_total` that are absent from `df_curves`.
"""
function _find_missing_sample_ids(df_curves::DataFrame, n_total::Int)::Vector{Int}
    completed = Set{Int}(df_curves.sample_id)
    return Int[i for i in 1:n_total if !(i in completed)]
end

_find_missing_sample_ids(::Nothing, n_total::Int)::Vector{Int} = collect(1:n_total)


"""
    _remove_checkpoint(path)

Remove a checkpoint file if it exists.

Parameters
----------
path : String
    Checkpoint path to remove.

Returns
-------
Nothing
    Always returns `nothing`.
"""
function _remove_checkpoint(path::String)
    isfile(path) && rm(path)
    return nothing
end
