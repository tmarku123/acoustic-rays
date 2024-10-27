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


struct MaterialProperties
    wave_velocity::Float64
    density::Float64
end

struct Zone
    id::Int64
    properties::MaterialProperties
end

struct Boundaries
    indices::Matrix{Union{Float64,Int}}
    normals::Matrix{Vector{Union{Float64,Int}}}
end

struct Domain
    zones::Vector{Zone}
    boundaries::Boundaries
    domain::Matrix{Union{Float64,Int}}
    scaling_factor::Int64
end

struct Parameters
    time::Vector{Float64}
end

struct Simulation
    parameters::Parameters
    domain::Domain
    rays::Vector{Ray}
    data::Array{Union{Float64,Int}}
end

struct Wavefront
    snapshot::Matrix{Float64}
    time::Float64
end



    



