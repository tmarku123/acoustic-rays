include("cartesian_structs.jl")
include("polar_structs.jl")

struct InitialConditions
    position::Vector{Union{Float64,Int}}
    angle::Union{Float64,Int}
    amplitude::Union{Float64,Int}
end

struct Ray
    position::Vector{Union{Float64, Int}}
    direction::Vector{Union{Float64, Int}}
    amplitude::Union{Float64, Int}
    length::Union{Float64, Int}
    initial_conditions::InitialConditions
    status::String
end

struct Parameters
    time::Vector{Float64}
end

struct Simulation
    parameters::Parameters
    domain::Union{Domain,PolarDomain}
    rays::Vector{Ray}
    data::Array{Union{Float64,Int}}
end

struct Wavefront
    snapshot::Matrix{Float64}
    time::Float64
end
