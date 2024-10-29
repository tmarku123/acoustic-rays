include("initialisations.jl")


function find_polar_zone(r::Float64,domain::PolarDomain)
    r_max = domain.r_max
    no_layers = domain.no_layers
    radii = collect(range(0,r_max,no_layers+1))[2:end]

    for i = 1:size(radii,1)
        if r < radii[i]
            return i
        end
    end
    
end

function is_in_polar_domain(ray::Ray,domain::PolarDomain)
    rad = norm(ray.position)
    rmax = domain.r_max
    if rad < rmax 
        return true
    end
end;

function check_polar_boundary(ray::Ray,domain::PolarDomain,dt::Float64)
    r = norm(ray.position)
    r_max = domain.r_max
    no_layers = domain.no_layers
    zone = find_polar_zone(r,domain)
    velocity = domain.properties[zone].velocity
    new_position = ray.position .+ ray.direction.*velocity.*dt 
    new_r = norm(new_position)
    new_zone = find_polar_zone(new_r,domain)
    if zone != new_zone
        return true
    end
    return false 
end;

function reflect_polar_ray(ray::Ray,domain::PolarDomain,dt::Float64)
    pos = ray.position 
    dir = ray.direction
    amp = ray.amplitude
    dir = dir ./ norm(dir)
    n = pos ./ norm(pos)
    outgoing = dir - 2*dot(dir,n)*n 
    outgoing = outgoing ./ norm(outgoing)
    angle = atand(outgoing[2]/outgoing[1])
    r = norm(pos)
    id = find_polar_zone(r,domain)
    vel = domain.properties[id].velocity
    new_pos = pos .+vel.*dir.*dt 
    new_id = find_polar_zone(norm(new_pos),domain)
    if !isnothing(new_id)
        velocity_2 = domain.properties[new_id].velocity
    else
        velocity_2 = 0
    end
    distance = outgoing.*vel.*dt 
    pos2 = pos + distance 
    IC = InitialConditions(pos2,angle,amp)
    length = norm(distance)
    damping = 0.8
    new_amp = amp*damping
    reflected_ray = Ray(pos2,outgoing,amp,length,IC,"active")
    return reflected_ray
end

function update_polar_ray(ray::Ray,domain::PolarDomain,dt::Float64)
    r = norm(ray.position)
    zone = find_polar_zone(r,domain)
    velocity = domain.properties[zone].velocity
    ds = ray.direction.*velocity.*dt
    new_pos = ray.position .+ ds 
    length = ray.length + norm(ds)
    dir = ray.direction ./ norm(ray.direction)
    IC = ray.initial_conditions
    a0 = IC.amplitude
    new_amp = a0*sqrt(norm(ds)/length)
    return Ray(new_pos,dir,new_amp,length,IC,"active")
end

function polar_snells_law(ray::Ray,domain::PolarDomain,dt::Float64)
    pos = ray.position
    dir = ray.direction
    amp = ray.amplitude
    r1 = norm(pos)
    id1 = find_polar_zone(r1,domain)
    velocity_1 = domain.properties[id1].velocity
    updated_ray = update_polar_ray(ray,domain,dt)
    r2 = norm(updated_ray.position)
    id2 = find_polar_zone(r2,domain)
    if !isnothing(id2)
        velocity_2 = domain.properties[id2].velocity
    else
        velocity_2 = 10^10
    end
    hyp = norm(dir)
    I = dir./ norm(dir)
    N = pos ./ norm(pos)
    ratio= velocity_2/velocity_1
    cosθ1 = dot(I,N)
    sinθ2_squared = ratio*(1-cosθ1^2)
    if sinθ2_squared <= 1
        cosθ2 = sqrt(1 - sinθ2_squared)
        cosθ2 = clamp(cosθ2, -1, 1)
        refracted_dir = ratio*I + (ratio*cosθ1-cosθ2)*N
        refracted_dir = refracted_dir / norm(refracted_dir)
        θ2 = acosd(cosθ2)
        new_position = pos + velocity_2.*refracted_dir.*dt
        a0 = ray.initial_conditions.amplitude
        IC = InitialConditions(new_position,θ2,a0)
        length = updated_ray.length
        new_amp = updated_ray.amplitude
        snells_ray = Ray(new_position,refracted_dir,new_amp,length,IC,"active")
        return snells_ray
    end
    return complete_ray(ray)
end

function compute_polar_rays(simulation::Simulation)
    t = simulation.parameters.time
    dt = t[2] - t[1]
    nt = size(t,1)
    domain = simulation.domain
    min_amp = domain.minimum_amplitude
    rays = simulation.rays
    data = simulation.data

    ray = rays[1]
    data = save_ray(ray,data,1,1)

    for i = 1:nt-1
        status = get_active_rays(rays)
        for j = 1:size(status)[1]
            if status[j] == 0
                continue
            end
            ray = rays[j]
            if is_in_polar_domain(ray,domain) == true
            
                amp = ray.amplitude
                if amp < min_amp
                    completed_ray = complete_ray(ray)
                    continue 
                end
                if check_polar_boundary(ray,domain,dt)
                    empty_indc = get_empty_ray(rays)
                    if !isnothing(empty_indc)
                        reflected_ray = reflect_polar_ray(ray,domain,dt)
                        data = save_ray(reflected_ray,data,empty_indc,i+1)
                        rays[empty_indc] = reflected_ray
                    end
                    snells_ray = polar_snells_law(ray,domain,dt)
                    data = save_ray(snells_ray,data,j,i+1)
                    rays[j] = snells_ray 
                else
                    updated_ray = update_polar_ray(ray,domain,dt)
                    data = save_ray(updated_ray,data,j,i+1)
                    rays[j] = updated_ray
                end
            else
                completed_ray = complete_ray(ray)
                rays[j] = completed_ray
            end
        end
    end
    return simulation
end

function obtain_polar_wavefronts(simulations::Vector{Simulation},scaling::Int64;blur::Bool=false)
    r_max = simulations[1].domain.r_max
    time = simulations[1].parameters.time
    nt = size(time,1)
    nx = Int(round(r_max*2*scaling))
    wavefronts = Vector{Wavefront}(undef,nt)  
    snapshot_template = zeros(nx,nx)
    for i = 1:nx,j = 1:nx
        pos = [i,j]
        i_unscaled = ((i-1)/scaling)-r_max
        j_unscaled = ((j-1)/scaling)-r_max

        r = norm([i_unscaled,j_unscaled])
        if r>r_max 
            snapshot_template[i,j] = NaN
        end
    end
    for i = 1:nt
        t = time[i]
        snapshot = copy(snapshot_template)
        for j = 1:size(simulations)[1]
            data = simulations[j].data[i,:,:]
            amp = data[3,:]
            z = data[1,:]
            x = data[2,:]
            indc = .!isnan.(amp)
            for k = 1:size(amp,1)
                if isnan(amp[k]) == true
                    continue
                end
                if norm([z[k],x[k]]) < r_max
                    x1 = get_scaled_points(x[k],r_max,nx)
                    z1 = get_scaled_points(z[k],r_max,nx)
                    if blur == true
                        snapshot = blur_amps(z1,x1,amp[k],snapshot,r_max)
                    else
                        x1 = Int(round(x1))
                        z1 = Int(round(z1))
                        x1 = clamp(x1, 1, nx);
                        z1 = clamp(z1, 1, nx);
                        snapshot[z1,x1] += amp[k]
                    end
                end
            end
        end

        # Define the standard deviation for the Gaussian blur (adjust for desired smoothness)
        sigma = 1  # e.g., 2.0 is a common starting value; increase for more blur
        # Apply Gaussian blur to the amplitude matrix
        #snapshot = imfilter(snapshot, Kernel.gaussian(sigma))
        wavefront = Wavefront(snapshot,t)
        wavefronts[i] = wavefront
    end
    return wavefronts 
end
