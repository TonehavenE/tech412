using SerialPorts, GLMakie
using GLMakie.GLFW
using GLMakie: to_native

s = SerialPort("/dev/ttyACM0", 9600)

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


# Observable to store sensor data
sensor_data = Observable(Packet(0, 0, 0, 0, 0, 0))

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
        
        sleep(0.25)  # Avoid CPU overuse
    end
end

function plot_live_data()
    # Observable to store sensor data
    port = SerialPort("/dev/ttyACM0", 9600)
	init_packet = read_packet(port)
    sensor_data = Observable(init_packet)
    Threads.@spawn read_serial(port, sensor_data)

    println("Finished spawning thread")
    
    fig = Figure()

    # Axes for different sensors

    ax_temp_sht30 = Axis(fig[1, 1], title="Temperature SHT30 (°C)")
    ax_light = Axis(fig[1, 2], title="Light Lux")
    ax_co2 = Axis(fig[2, 1], title="CO2 ppm")
    ax_temp_scd4x = Axis(fig[2, 2], title="Temperature SCD4x (°C)")
    ax_humidity_scd4x = Axis(fig[3, 1], title="Humidity SCD4x (%)")
    ax_humidity_sht30 = Axis(fig[3, 2], title="Humidity SHT30 (%)")

    # Time series storage
    start_time = time()
	times = [start_time]
	temp_sht30 = [init_packet.temperature_sht30]
	light_lux = [init_packet.light_lux]
	co2_ppm = [init_packet.co2_ppm]
	temp_scd4x = [init_packet.temperature_scd4x]
	humidity_scd4x = [init_packet.relative_humidity_scd4x]
	humidity_sht30 = [init_packet.humidity_sht30]

	on(sensor_data) do pkt
		println("New Packet recieved: $pkt")
		push!(times, time() - start_time)
		push!(temp_sht30, pkt.temperature_sht30)
        push!(light_lux, pkt.light_lux)
        push!(co2_ppm, pkt.co2_ppm)
        push!(temp_scd4x, pkt.temperature_scd4x)
        push!(humidity_scd4x, pkt.relative_humidity_scd4x)
        push!(humidity_sht30, pkt.humidity_sht30)
	end
    display(fig)

	glfw_window = to_native(display(fig))

	on(events(fig).keyboardbutton) do event
    	if event.key == Keyboard.q
      		GLFW.SetWindowShouldClose(glfw_window, true) # this will close the window after all callbacks are finished
		end
    end

    lines!(ax_temp_sht30, times, temp_sht30)
    lines!(ax_light, times, light_lux)
    lines!(ax_co2, times, co2_ppm)
    lines!(ax_temp_scd4x, times, temp_scd4x)
    lines!(ax_humidity_scd4x, times, humidity_scd4x)
    lines!(ax_humidity_sht30, times, humidity_sht30)
end

plot_live_data()
