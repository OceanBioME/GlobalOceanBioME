@inline function time_indices_weights(time, times)
    day_number = mod(time / day, 365)

    if day_number <= times[1]
        t₁, t₂ = times[12], times[1] 
        Δt = t₂ + 365 - t₁

        return 12, 1, (day_number + 365 - t₁) / Δt
    elseif day_number >= times[end]
        t₁, t₂ = times[12], times[1] 
        Δt = t₂ + 365 - t₁
        return 12, 1, (day_number - t₁) / Δt
    else
        n₁ = 1
        while day_number > times[n₁ + 1]
            n₁ += 1
        end

        t₁, t₂ = times[n₁], times[n₁ + 1] 

        Δt = t₂ - t₁

        return n₁, n₁ + 1, (day_number - t₁) / Δt
    end
end

@inline function cyclic_interpolate(data, i, j, k, time, times)
    n₁, n₂, weight = time_indices_weights(time, times)

    @inbounds begin
        d₁ = data[i, j, k, n₁]
        d₂ = data[i, j, k, n₂]
    end

    return d₁ + (d₂ - d₁) * weight
end

@inline function cyclic_interpolate(data, i, j, time, times)
    n₁, n₂, weight = time_indices_weights(time, times)

    @inbounds begin
        d₁ = data[i, j, n₁]
        d₂ = data[i, j, n₂]
    end

    return d₁ + (d₂ - d₁) * weight
end