# Native integration-Newton verification; no physical rates or states are modified.
#
# IDA permits a user convergence test on an attached SUNNonlinSol_Newton:
# https://sundials.readthedocs.io/en/latest/sunnonlinsol/SUNNonlinSol_links.html
# Retaining IDA's original test additionally requires its generated Newton-content
# binding. That ABI is explicitly restricted below; unsupported versions fail
# before stepping, rather than silently omitting verification.
#
# This is a discrete nonlinear-solve accuracy check, not a temporal error estimate,
# a guarantee of global conservation, or a proof of nonnegative gas inventories.

import Sundials

Base.@kwdef mutable struct DAENewtonVerification
    owner::Any
    verify::Any
    n::Int
    nls::Any = nothing
    callback::Any = nothing
    original_test::Ptr{Cvoid} = C_NULL
    original_data::Ptr{Cvoid} = C_NULL
    installed::Bool = false
    released::Bool = false
    failure::Union{Nothing,String} = nothing
    calls::Int = 0
    checks::Int = 0
    vetoes::Int = 0
    verified_successes::Int = 0
    max_verified_wrms::Float64 = 0.0
    max_verified_component::Float64 = 0.0
    accepted_steps::Int = 0
    last_step::Int = 0
    max_endpoint_difference_scaled::Float64 = 0.0
    last_verified::Any = nothing
end

"""Estimate the remaining correction with private residual/Jacobian workspaces.

The coupled DAE Jacobian includes cj*dF/d(du). Both the WRMS nonlinear budget
and a unit component budget must hold. No predicted phase-sign veto is used.
The estimate is never applied to the integrated state.
"""
function _build_newton_verifier(residual!, packed, prototype, differential_vars)
    J = copy(prototype)
    cache = _build_jacobian_coloring_cache(prototype, differential_vars)
    res = Vector{Float64}(undef, size(J, 1))
    return function (y, dy, cj, t, weights, tolerance)
        all(isfinite, y) && all(isfinite, dy) &&
            all(x -> isfinite(x) && x > 0, weights) ||
            error("Nonfinite DAE Newton verification input.")
        isfinite(tolerance) && tolerance > 0 ||
            error("Invalid native nonlinear convergence tolerance.")
        residual!(res, dy, y, packed, t)
        all(isfinite, res) || error("Nonfinite DAE residual during Newton verification.")
        _dae_jacobian_fd!(J, dy, y, packed, cj, t, residual!, differential_vars, cache)
        correction = -(J \ res)
        all(isfinite, correction) || error("Nonfinite fresh DAE Newton correction.")
        q = abs.(correction) .* weights
        wrms = sqrt(sum(abs2, q) / length(q))
        component = maximum(q)
        return (; passed=wrms <= tolerance && component <= 1.0, wrms, component)
    end
end

function _newton_native_vector(ptr, n)
    ptr == C_NULL && error("Null native vector during Newton verification.")
    vector = Sundials.N_Vector(ptr)
    Sundials.N_VGetLength(vector) == n || error("Newton verification vector length mismatch.")
    data = Sundials.N_VGetArrayPointer(vector)
    data == C_NULL && error("Null native Newton vector data.")
    return unsafe_wrap(Vector{Float64}, data, n; own=false)
end

function _newton_system_data(g)
    t = Ref{Cdouble}(0)
    cj = Ref{Cdouble}(0)
    yp = Ref{Ptr{Cvoid}}(C_NULL); dp = Ref{Ptr{Cvoid}}(C_NULL)
    yn = Ref{Ptr{Cvoid}}(C_NULL); dn = Ref{Ptr{Cvoid}}(C_NULL)
    rr = Ref{Ptr{Cvoid}}(C_NULL); user = Ref{Ptr{Cvoid}}(C_NULL)
    flag = ccall((:IDAGetNonlinearSystemData, Sundials.libsundials_idas), Cint,
        (Sundials.IDAMemPtr, Ref{Cdouble}, Ref{Ptr{Cvoid}}, Ref{Ptr{Cvoid}},
         Ref{Ptr{Cvoid}}, Ref{Ptr{Cvoid}}, Ref{Ptr{Cvoid}}, Ref{Cdouble},
         Ref{Ptr{Cvoid}}),
        g.owner.mem, t, yp, dp, yn, dn, rr, cj, user)
    flag == 0 || error("IDAGetNonlinearSystemData failed: $flag")
    return (; t=Float64(t[]), cj=Float64(cj[]),
        predictor=_newton_native_vector(yp[], g.n),
        derivative_predictor=_newton_native_vector(dp[], g.n))
end

function _dae_newton_convergence_test(g, nls, ycor, delta, tol, ewt, data)
    # No Julia exception may unwind through C. The integration wrapper surfaces
    # a stored diagnostic after the native solve stops.
    try
        g.calls += 1
        original = ccall(g.original_test, Cint,
            (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Cdouble, Ptr{Cvoid}, Ptr{Cvoid}),
            nls, ycor, delta, tol, ewt, g.original_data)
        original == 0 || return Cint(original)
        g.checks += 1
        info = _newton_system_data(g)
        weights = copy(_newton_native_vector(ewt, g.n))
        correction = copy(_newton_native_vector(ycor, g.n))
        update = copy(_newton_native_vector(delta, g.n))
        y = copy(info.predictor) .+ correction
        dy = copy(info.derivative_predictor) .+ info.cj .* correction
        assessment = g.verify(y, dy, info.cj, info.t, weights, tol)
        isequal(_newton_native_vector(ycor, g.n), correction) &&
            isequal(_newton_native_vector(delta, g.n), update) &&
            isequal(_newton_native_vector(ewt, g.n), weights) ||
            error("Private verification changed native Newton vectors.")
        if !assessment.passed
            g.vetoes += 1
            return Cint(Sundials.SUN_NLS_CONTINUE)
        end
        g.verified_successes += 1
        g.max_verified_wrms = max(g.max_verified_wrms, assessment.wrms)
        g.max_verified_component = max(g.max_verified_component, assessment.component)
        g.last_verified = (; t=info.t, tolerance=Float64(tol),
            wrms=assessment.wrms, component=assessment.component, y, dy)
        return Cint(original)
    catch err
        g.failure = sprint(showerror, err)
        return Cint(-1)
    end
end

"""Check the explicitly reviewed native ABI before attaching a custom solver.

Project compatibility alone does not establish the layout of native content.
Extending this list requires testing the corresponding Julia and native binding.
No unchecked fallback is permitted.
"""
function _check_newton_verification_backend()
    Base.pkgversion(Sundials) == v"5.3.0" &&
        Base.pkgversion(Sundials.Sundials_jll) == v"7.5.0+0" ||
        error("DAE Newton verification currently supports Sundials.jl 5.3.0 / " *
              "SUNDIALS 7.5.0+0; review the native convergence-test binding for this installation.")
    required = (:SUNNonlinSol_Newton, :IDASetNonlinearSolver,
        :SUNNonlinSolSetConvTestFn, :SUNNonlinSolFree,
        :_SUNNonlinearSolverContent_Newton, :IDASetMaxNonlinIters)
    all(name -> isdefined(Sundials, name), required) ||
        error("Required native Newton verification interface is unavailable.")
    fieldnames(Sundials._SUNNonlinearSolverContent_Newton) ==
        (:Sys, :LSetup, :LSolve, :CTest, :delta, :jcur, :curiter, :maxiters,
         :niters, :nconvfails, :ctest_data) ||
        error("Unexpected native Newton convergence-test layout.")
    return nothing
end

"""Attach verification after IDA initialization and before its first step."""
function _install_newton_verification!(slot, integrator, verifier)
    _check_newton_verification_backend()
    slot[] === nothing || error("Newton verification was already installed.")
    g = DAENewtonVerification(owner=integrator, verify=verifier, n=length(integrator.u))
    slot[] = g # Also allows exception-safe cleanup after a partial installation.
    g.nls = Sundials.SUNNonlinSol_Newton(integrator.u_nvec, integrator.u_nvec.ctx)
    g.nls == C_NULL && error("SUNNonlinSol_Newton allocation failed.")
    flag = Sundials.IDASetNonlinearSolver(integrator.mem, g.nls)
    flag == 0 || error("IDASetNonlinearSolver failed: $flag")
    content = unsafe_load(Ptr{Sundials._SUNNonlinearSolverContent_Newton}(
        unsafe_load(g.nls).content))
    g.original_test = content.CTest
    g.original_data = content.ctest_data
    g.original_test != C_NULL && g.original_data != C_NULL ||
        error("Original IDA convergence callback is unavailable.")
    flag = Sundials.IDASetMaxNonlinIters(integrator.mem,
        Cint(integrator.alg.max_nonlinear_iters))
    flag == 0 || error("Restoring IDA nonlinear iteration limit failed: $flag")
    cb = (a,b,c,d,e,f) -> _dae_newton_convergence_test(g,a,b,c,d,e,f)
    g.callback = @cfunction($cb, Cint,
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Cdouble, Ptr{Cvoid}, Ptr{Cvoid}))
    flag = Sundials.SUNNonlinSolSetConvTestFn(g.nls,
        Base.unsafe_convert(Ptr{Cvoid}, g.callback), g.original_data)
    flag == 0 || error("Installing DAE Newton verification failed: $flag")
    g.installed = true
    return g
end

"""Verify linkage to an accepted native endpoint, after native constraint processing."""
function _check_verified_newton_endpoint!(slot, integrator, native_y, native_dy)
    g = slot[]
    g === nothing && error("Accepted endpoint has no installed Newton verification.")
    g.failure === nothing || error("DAE Newton verification failed: " * g.failure)
    v = g.last_verified
    v === nothing && error("Accepted endpoint has no verified Newton success.")
    t = Float64(integrator.t)
    abs(v.t - t) <= 64 * eps(Float64) * max(1.0, abs(t), abs(v.t)) ||
        error("Accepted endpoint has no same-time verified Newton success.")
    v.wrms <= v.tolerance && v.component <= 1.0 ||
        error("Accepted endpoint failed the Newton correction budgets.")
    step = Ref{Clong}(0)
    flag = Sundials.IDAGetNumSteps(integrator.mem, step)
    flag == 0 || error("IDAGetNumSteps failed: $flag")
    current = Int(step[])
    if current != g.last_step
        current == g.last_step + 1 || error("Newton verification missed an accepted step.")
        g.accepted_steps += 1
        g.last_step = current
    end
    # Constraint correction is native, not applied by this verifier. Its effect
    # must remain visible; pre-constraint success is not an endpoint residual test.
    g.max_endpoint_difference_scaled = max(g.max_endpoint_difference_scaled,
        maximum(abs.(native_y .- v.y)))
    return nothing
end

function _release_newton_verification!(slot)
    g = slot[]
    g === nothing && return nothing
    g.released && return nothing
    try
        if g.nls !== nothing && g.nls != C_NULL
            owner = g.owner
            flag = GC.@preserve owner g begin
                Sundials.SUNNonlinSolFree(g.nls)
            end
            flag == 0 || error("SUNNonlinSolFree failed: $flag")
            g.nls = nothing
        end
        g.released = true
    finally
        # IDA borrows the custom NLS. Free it while its SUNContext owner is alive.
        # No further IDASolve call is allowed on this integrator after release.
        g.callback = nothing
        g.owner = nothing
        g.verify = nothing
        g.last_verified = nothing
    end
    return nothing
end

function _newton_verification_stats(slot)
    g = slot[]
    g === nothing && return nothing
    return (; installed=g.installed, released=g.released, failure=g.failure,
        callbacks=g.calls, checks=g.checks, vetoes=g.vetoes,
        verified_successes=g.verified_successes, accepted_steps=g.accepted_steps,
        maximum_verified_wrms=g.max_verified_wrms,
        maximum_verified_component=g.max_verified_component,
        maximum_candidate_endpoint_difference_scaled=g.max_endpoint_difference_scaled,
        implementation="productive DAE Newton verification")
end
