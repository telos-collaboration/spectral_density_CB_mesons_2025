function madras_sokal_estimator_fixedt(x, t)
    m = mean(x)
    Γ = zero(eltype(x))
    N = length(x)
    for i in 1:N-t
        Γ += (x[i]-m)*(x[i+t]-m)/(N-t)
    end
    return Γ 
end
# With the default maximal window size tmax=length(x)÷10,  
# the Madras-Sokal variance estimate is such that at tmax
# dτ ÷ τ = 1/sqrt(2.5) ≈ 0.63
function madras_sokal_estimator_windows(x;max_window=length(x)÷10)
    Γ = zeros(eltype(x),max_window)
    for t in 1:max_window
        Γ[t] = madras_sokal_estimator_fixedt(x, t)
    end
    return Γ 
end
function madras_sokal_windows(x;kws...)
    Γ  = madras_sokal_estimator_windows(x;kws...)
    τ  = 1/2 .+ cumsum(Γ/Γ[1])
    dτ = similar(τ)
    N  = length(x)
    for i in eachindex(τ)
        dτ[i] = sqrt(τ[i]^2 * (4i+2)/N)
    end
    return τ, dτ
end
function madras_sokal_time(x,therms;stop=length(x).-therms,kws...)
    τ  = zeros(Float64,size(therms))
    dτ = zeros(Float64,size(therms))
    for j in eachindex(therms)
        # default value of stop is chosen, such that it is the end of x
        xc = x[therms[j]:therms[j]+stop[j]]
        τ_windows, dτ_windows = madras_sokal_windows(xc;kws...)
        τ[j], W = findmax(τ_windows)
        dτ[j]   = dτ_windows[W]
    end
    return τ, dτ
end
function madras_sokal_time(x;kws...)
    τ_windows, dτ_windows = madras_sokal_windows(x;kws...)
    τ, W = findmax(τ_windows)
    dτ   = dτ_windows[W]
    return τ, dτ
end
