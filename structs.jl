struct InitialConditions
    position::Vector{Float64}
    angle::Float64
    amplitude::Float64
end

struct Ray
    position::Vector{Float64}
    direction::Vector{Float64}
    amplitude::Float64
    length::Float64
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
    indices::Matrix{Float64}
    normals::Matrix{Vector{Float64}}
end

struct Domain
    zones::Vector{Zone}
    boundaries::Boundaries
    domain::Matrix{Float64}
    scaling_factor::Int64
end

struct Parameters
    time::Vector{Float64}
end

struct Simulation
    parameters::Parameters
    domain::Domain
    rays::Vector{Ray}
    data::Array{Float64}
end

struct Wavefront
    snapshot::Matrix{Float64}
    time::Float64
end



    



