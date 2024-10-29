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







