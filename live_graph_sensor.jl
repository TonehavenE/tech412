using LinearAlgebra, GLMakie, GLMakie.GLFW, SerialPorts
using GLMakie: to_native
using Statistics
using Base.Threads

struct Packet
    temperature_sht30::Float64
    light_lux::Float64
    co2_ppm::Float64
    temperature_scd4x::Float64
    relative_humidity_scd4x::Float64
    humidity_sht30::Float64
end

function parse_packet(data::Vector{Float64})
    if length(data) == 6
        return Packet(data...)
    else
        error("Invalid packet length: expected 6 values, got $(length(data))")
    end
end

function read_packet(port)
    buffer = ""
    
    while true
        data = readavailable(port)  # Read incoming data as a string
        buffer *= data  # Append new data to the buffer
    
        # Look for a full packet
        m = match(r"START, (.*?), END", buffer)
        
        if m !== nothing
            extracted_data = split(m.captures[1], ", ")  # Split values by ", "
            parsed_data = tryparse.(Float64, extracted_data)  # Convert to numbers
            return parse_packet(parsed_data)
        end
    
        sleep(0.01)  # Small delay to prevent excessive CPU usage
    end
end

function read_serial(port::SerialPort, sensor_data::Observable{Packet})
    buffer = ""

    while true
        data = readavailable(port)  # Read serial data
        buffer *= data  # Append to buffer

        # Look for a full packet
        m = match(r"START, (.*?), END", buffer)
        
        if m !== nothing
            extracted_data = split(m.captures[1], ", ")
            parsed_data = tryparse.(Float64, extracted_data)

            if length(parsed_data) == 6
                sensor_data[] = Packet(parsed_data...)  # Update Observable
                println("Updated packet: ", sensor_data[])
            end

            # Remove processed data from buffer
            buffer = buffer[findfirst("END", buffer)[end]+1:end]  
        end
        
        sleep(0.05)  # Avoid CPU overuse
    end
end


"""
    rotation(θ)

return 2 x 2 rotation matrix that rotates plane by angle θ
"""
function rotation(θ)
    c, s = cos(θ), sin(θ)
    return [c -s; s c]
end

"""
    dihedralgroup(n, flip=true)

Generate symmetries of the n-gon that rotate and flip the n-gon in the plane.
The return value is an array of n or 2n matrices representing the elements of the group.
"""
function dihedralgroup(n, flip=true)    
    S = [-1 0; 0 1]                 # a reflection about y axis
    I = [1.0 0.0; 0.0 1.0]          # the identity
    Dn = fill(I, flip ? 2n : n) # allocate an array of 2n or n matrices
    
    for k=1:n
        @views Dn[k] = rotation(2(k-1)π/n) # set Dn[k] to rotation by θ = 2(k-1)π/n
        if flip
            Dn[k+n] = S*Dn[k]       # set Dn[k+n] to reflection of Dn[k]
        end
    end
    Dn
end

"""
    symmetrize(X, G)

Symmetrize a set of data points X by symmetry group G. The return value
is a matrix containing all columns of X mapped by all matrices in G
"""
function symmetrize(X, G)
    m, nX = size(X)
    nG = length(G)
    GX = similar(X, m, nX * nG)  # Preallocate memory

    for j in 1:nX
        @views Xj = X[:, j]  # Avoid unnecessary indexing
        for k in 1:nG
            GX[:, (j-1) * nG + k] = G[k] * Xj
        end
    end
    return GX
end


"""
    f(x, X, a=1, k=1)

return 1/N sum_j cos(k|x-xj|) exp(-a|x-xj|^2) where 
  x is a 2d vector 
  xj is the jth column of 2 x N matrix X
"""
function f(x, X, a=1, k=1)
	norms = vec(norm.(eachcol(X .- x)))
	return sum(cos.(k .* norms) .* exp.(-a .* norms.^2))/ size(X, 2)
end

function generate_zgrid(Xsymm, a, k, x1grid, x2grid)
	N1 = length(x1grid)
	N2 = length(x2grid)
	zgrid = zeros(N1, N2)

	@threads for i in 1:N1
		for j in 1:N2
			x = [x1grid[i]; x2grid[j]]
			zgrid[i, j] = f(x, Xsymm, a, k)
		end
	end
	return zgrid
end

function compute_pattern(X, a, k, w, n, flip)	
    G = dihedralgroup(n, flip)
    Xsymm = symmetrize(X, G)

    x1grid = range(-w, w, length=100)
    x2grid = range(-w, w, length=100)

    # Generate zgrid based on the input data and current parameters
    zgrid = generate_zgrid(Xsymm, a, k, x1grid, x2grid)
    zscale = maximum(abs.(zgrid))

    return zgrid / zscale
end


function map_range(x, a_s, b_s, a_d, b_d)
    return (b_d - a_d) * (x - a_s) / (b_s - a_s) + a_d
end

"""
Create an animation with sliders.
"""
function animation_with_sliders()
    w = 400
    fps = 60
    nframes = 240
    time = 0
    polygon = 5
    a0 = 3    # scale of blobs (larger a, narrower blobs)
    k0 = 5    # scale of ripples (larger k, more rapid ripples)
    s = 3     # scale of data points (larger s, further spread out)
	Npts = 3
	width = 3
	flip = true
	levels = -1:.3:1  # number or values of contour levels

	X = randn(2, Npts)
    Nx = size(X,2)
    Xt = copy(X)
	dt = pi/512

	# some choices for rotation rates
	#ω = 1.5*(2*rand(Npts).-1) # random rotation rates for data points 
	#ω = ω .- sum(ω)/Npts      # remove mean rotation
	#ω = rand(Int, Npts) .% 6 .+ 1
	ω = (rand(Int,Npts) .% 4 .+ 1).* (-1).^(rand(Int,Npts) .% 2)
	if !flip
	    ω = ω .- sum(ω)/Npts      # remove mean rotation
	end
	println(ω)
	R = rotation.(ω*dt)   # vector of incremental rotation matrices, one for each data point


	a(a0, t) = a0*(1 .- 2/3*cos.(t))
	k(k0, t) = k0*(1 .- 2/3*cos.(2*t));

    zgrid = Observable(fill(0.0, 100, 100))  # Placeholder, will be updated with actual data

	# create figure
	fig = Figure(size=(1920, 1080))
    ax = Axis(fig[1, 1], aspect = 1, xgridvisible = false, ygridvisible = false)

	port = SerialPort("/dev/ttyACM0", 9600)
	sensor_data = Observable(Packet(0, 0, 0, 0, 0, 0))
	serial_task = Threads.@spawn read_serial(port, sensor_data)

    selected_colormap = "viridis"
    contourf!(ax, range(-w, w, length=100), range(-w, w, length=100), zgrid, colormap=selected_colormap, levels=levels)

    glfw_window = to_native(display(fig))
    on(events(fig).keyboardbutton) do event
        if event.key == Keyboard.q
            GLFW.SetWindowShouldClose(glfw_window, true)
        end
    end
    # for i = 1:nframes
	history_packets = []
	window_size = 4
	pkt = 0
	polygon = 5
	a0 = 0
	k0 = 0
	while !GLFW.WindowShouldClose(glfw_window)
		if pkt != sensor_data[]
			pkt = sensor_data[]
			push!(history_packets, pkt)
			if length(history_packets) > window_size
				popfirst!(history_packets)
			end
			temp_sht30 = mean(history_packets[i].temperature_sht30 for i in eachindex(history_packets))
			light = mean(history_packets[i].light_lux for i in eachindex(history_packets))
			co2 = mean(history_packets[i].co2_ppm for i in eachindex(history_packets))
			temp_scd4x = mean(history_packets[i].temperature_scd4x for i in eachindex(history_packets))
			hum_scd4x = mean(history_packets[i].relative_humidity_scd4x for i in eachindex(history_packets))
			hum_sht30 = mean(history_packets[i].humidity_sht30 for i in eachindex(history_packets))

			a0 = floor(map_range(hum_sht30, 1, 100, 1, 10))
			k0 = floor(map_range(light, 1, 4000, 1, 10))
			ω_0 = floor(map_range(hum_sht30, 1, 100, 0, pi/2)) + floor(map_range(light, 1, 1500, 0, pi/2))
			ω = repeat([ω_0], 3)

			# polygon = floor(Int, map_range(light, 1, 250, 2, 10))
			println("history_packets: $history_packets")
			println("Sensor Data: $pkt")
			println("Polygon: $polygon")
			println("A0: $a0")
			println("K0: $k0")
			println("ω: $ω")
		end
        new_zgrid = compute_pattern(Xt, a(a0, time), k(k0, time), width, polygon, flip)

        # Update zgrid observable
        zgrid[] = new_zgrid

		sleep(1/fps)
		time += 1/fps

		for j=1:Nx
			Xt[:, j] = rotation.(ω * dt)[j]*Xt[:, j]
		end
    end

	# Prevent script from exiting
	wait(serial_task)
    return fig
end

animation_with_sliders()

