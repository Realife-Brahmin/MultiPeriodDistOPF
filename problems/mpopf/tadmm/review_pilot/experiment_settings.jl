# July release 50b8901c324ba9eefd0d672c6e484e0ededc7a57: config.jl.
const JULY_CB_SOURCE = "50b8901c324ba9eefd0d672c6e484e0ededc7a57"

function july_battery_cb(T::Int)
    T>=1 || error("T must be positive")
    prices=T==1 ? [0.14] : 0.08 .+ 0.12 .* (sin.(range(0,2pi;length=T)) .+ 1) ./ 2
    return 1e-6*minimum(prices)
end

function apply_review_cb!(d,system;profile=get(ENV,"TADMM_CB_PROFILE","july"))
    profile in ("july","ddp") || error("TADMM_CB_PROFILE must be july or ddp")
    d[:C_B]=profile=="july" ? july_battery_cb(d[:T]) : battery_cb(system,"system")
    return profile
end
