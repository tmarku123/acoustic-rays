struct Property
    velocity::Union{Float64,Int64}
    density::Union{Float64,Int64}
end

struct PolarDomain
    r_max::Union{Float64,Int64}
    theta_max::Union{Float64,Int64}
    no_layers::Int64
    properties::Vector{Property}
    minimum_amplitude::Float64
end

