# -*- coding: utf-8 -*-

"""
    CVExtraction I/O

Read cyclic voltammetry files in ZAHNER `.txt` and BioLogic `.mpt` formats.
"""

using CSV
using DataFrames

"""
    _count_header_lines(path::String) -> Int

Return the number of lines to skip before the column-header line in a text
data file. The first fully numeric line marks the beginning of the data.
"""
function _count_header_lines(path::String)::Int
    open(path, "r") do io
        n = 0
        for line in eachline(io)
            tokens = split(line)
            if isempty(tokens)
                n += 1
                continue
            end
            # A data line is assumed to contain only numeric tokens
            numeric = true
            for tok in tokens
                if tryparse(Float64, tok) === nothing
                    numeric = false
                    break
                end
            end
            if numeric
                return n - 1
            end
            n += 1
        end
        return n - 1
    end
end

"""
    read_cv_file(path::String) -> CVData

Read a CV file and return a `CVData` object.

Supported formats:
- ZAHNER `.txt` with columns `Time/s`, `Potential/V`, `Current/A`
- BioLogic `.mpt` with columns `time/s`, `Ewe/V`, `<I>/mA`
"""
function read_cv_file(path::String)::CVData
    ext = lowercase(splitext(path)[2])
    n_header = _count_header_lines(path)

    if ext == ".mpt"
        return _read_biologic_mpt(path, n_header)
    elseif ext == ".txt"
        return _read_zahner_txt(path, n_header)
    else
        error("Unsupported CV file extension: $ext (expected .txt or .mpt)")
    end
end

function _read_zahner_txt(path::String, n_header::Int)::CVData
    df = CSV.read(
        path,
        DataFrame;
        header = n_header + 1,
        delim = ' ',
        ignorerepeated = true,
        stripwhitespace = true,
        types = Float64,
    )

    colnames = names(df)
    time_col = _find_column(colnames, ["Time/s", "time/s", "t/s", "Time"])
    potential_col = _find_column(colnames, ["Potential/V", "potential/V", "U/V", "UinV"])
    current_col = _find_column(colnames, ["Current/A", "current/A", "I/A", "IinA"])

    t = df[:, time_col]
    U = df[:, potential_col]
    I = df[:, current_col]

    return CVData(t, U, I)
end

function _read_biologic_mpt(path::String, n_header::Int)::CVData
    df = CSV.read(
        path,
        DataFrame;
        header = n_header + 1,
        delim = '\t',
        stripwhitespace = true,
        types = Float64,
    )

    colnames = names(df)
    time_col = _find_column(colnames, ["time/s", "Time/s", "t/s"])
    potential_col = _find_column(colnames, ["Ewe/V", "Potential/V", "U/V"])
    current_col = _find_column(colnames, ["<I>/mA", "I/mA", "Current/mA"])

    t = df[:, time_col]
    U = df[:, potential_col]
    I = df[:, current_col] ./ 1000.0  # mA -> A

    return CVData(t, U, I)
end

function _find_column(colnames::Vector{String}, candidates::Vector{String})::String
    for cand in candidates
        idx = findfirst(==(cand), colnames)
        if idx !== nothing
            return colnames[idx]
        end
    end
    error("Could not find any of the expected columns $(candidates) in $colnames")
end
