include("../structs/general_structs.jl")

using LinearAlgebra

function get_zone_id(ray::Ray,domain::Domain)
    (z,x) = ray.position
    (nz,nx) = size(domain.domain)
    scaling_factor = domain.scaling_factor
    x = clamp(Int(round(x*scaling_factor)),1,nx);
    z = clamp(Int(round(z*scaling_factor)),1,nz);
    id = domain.domain[z,x]
    return Int(id)
end;

function get_velocity(ray,domain::Domain)
    zone_id = get_zone_id(ray,domain)
    velocity = domain.zones[zone_id].properties.wave_velocity
    return velocity
end;

function compute_boundaries(domain::Matrix{Float64})
    nz,nx = size(domain)
    boundaries = zeros(nz,nx)
    for i = 2:nz-1, j = 2:nx-1
        submatrix = domain[(i-1:i+1),(j-1:j+1)]
        if sum(submatrix) != domain[i,j]*9
            boundaries[i,j] = 1
        end
    end
    return boundaries
end;

function get_active_rays(rays::Vector{Ray})
    status = [ray.status for ray in rays] .== "active"
    return status
end;

function get_empty_ray(rays::Vector{Ray})
    status = [ray.status for ray in rays] .== "inactive"
    indc = findfirst(status .== 1)
    return indc 
end;

function complete_ray(ray::Ray)
    pos = ray.position
    dir = ray.direction
    amp = ray.amplitude
    length = ray.length
    IC = ray.initial_conditions
    completed_ray = Ray(pos,dir,amp,length,IC,"completed")
    return completed_ray
end;

function save_ray(ray::Ray,data::Array{},ray_index::Int64,time_index::Int64)
    (z,x) = ray.position 
    amp = ray.amplitude
    data[time_index,:,ray_index] .= [z,x,amp]
    return data
end;

function get_scaled_points(x,r_max::Float64,nx::Int64)
    scaling = (nx/(2*r_max))
    x = (x + r_max)*scaling + 1
    return x
end

function blur_amps(z::Float64,x::Float64,amp::Float64,snapshot::Array{Float64},r_max::Float64)
    nx=size(snapshot,1)
    x_floor, z_floor = floor(x), floor(z)
    x_frac, z_frac = x-x_floor, z-z_floor
    # Calculate weights for each of the four neighboring grid points
    w_tl = (1 - x_frac) * (1 - z_frac)  # top-left
    w_tr = x_frac * (1 - z_frac)        # top-right
    w_bl = (1 - x_frac) * z_frac        # bottom-left
    w_br = x_frac * z_frac              # bottom-right

    z_floor = Int(clamp(z_floor,1,nx-1))
    x_floor = Int(clamp(x_floor,1,nx-1))

    # Distribute the amplitude to the four nearest grid points
    snapshot[z_floor, x_floor] += amp * w_tl    # top-left
    snapshot[z_floor, x_floor + 1] += amp * w_tr  # top-right
    snapshot[z_floor + 1, x_floor] += amp * w_bl  # bottom-left
    snapshot[z_floor + 1, x_floor + 1] += amp * w_br  # bottom-right

    return snapshot
end

function get_number_of_rays_used(simulation::Simulation)
    data = simulation.data 
    for i = 1:size(data,3)
        data_slice = data[:,:,i]
        if maximum(.!isnan.(data_slice)) == false
            return i 
        end
    end
    return size(data,3)
end
