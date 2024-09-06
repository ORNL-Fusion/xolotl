#!/usr/bin/julia

using Plots

struct MemEvent
    time::Float64
    ptr::String
    size::Int64
    space::String
    name::String
end

function MemEvent(ser::Vector{<:AbstractString})
    name = "no_name"
    if length(ser) == 6
        name = ser[6]
    end
    MemEvent(parse(Float64, ser[1]), ser[2], parse(Int64, ser[3]), ser[4], name)
end

function readMemEvents(filename::AbstractString; memspace = "")
    ret = MemEvent[]
    lines = readlines(filename)
    sizehint!(ret, length(lines))
    for line in lines
        if startswith(line, "#")
            continue
        end
        find_f = occursin(line)
        if find_f("PushRegion") || find_f("PopRegion")
            continue
        end
        event = MemEvent(split(line))
        if isempty(memspace) || event.space == memspace
            push!(ret, event)
        end
    end
    return ret
end

struct NamedSeries
    name::String
    events::Vector{MemEvent}
    time::Vector{Float64}
    total::Vector{Float64}
end

function NamedSeries(name::String, events::Vector{MemEvent}, size::Int64)
    ret = NamedSeries(name, events, [0.0], [0.0])
    sizehint!(ret.time, size + 1)
    sizehint!(ret.total, size + 1)
    return ret
end

NamedSeries(name::String, size::Int64) = NamedSeries(name, MemEvent[], size)

@inline Base.isequal(a::NamedSeries, b::NamedSeries) = isequal(a.name, b.name)

function getTotalSeries(events::Vector{MemEvent})
    ret = NamedSeries("Total", events, [0.0], [0.0])
    sizehint!(ret.time, length(events) + 1)
    sizehint!(ret.total, length(events) + 1)
    for e in events
        push!(ret.time, e.time)
        push!(ret.total, last(ret.total) + e.size)
    end
    return ret
end

function getAllNamedSeries(totalSeries::NamedSeries)
    map = Dict{String,NamedSeries}()
    events = totalSeries.events
    numEvents = length(events)
    ret = Vector{NamedSeries}()
    for i in eachindex(events)
        e = events[i]
        ser = get!(map, e.name, NamedSeries(e.name, numEvents))
        if findfirst(isequal(ser), ret) == nothing
            push!(ret, ser)
        end
        push!(ser.events, e)
        append!(ser.time, totalSeries.time[length(ser.time)+1 : i+1])
        append!(ser.total, repeat([last(ser.total)], i - length(ser.total)))
        push!(ser.total, last(ser.total) + e.size)
    end

    for ser in ret
        i = numEvents + 1
        append!(ser.time, totalSeries.time[length(ser.time)+1 : end])
        append!(ser.total, repeat([last(ser.total)], i - length(ser.total)))
    end

    return ret
end

@userplot StackedArea

# a simple "recipe" for Plots.jl to get stacked area plots
# usage: stackedarea(xvector, datamatrix, plotsoptions)
@recipe function f(pc::StackedArea)
    x, y = pc.args
    n = length(x)
    y = cumsum(y, dims=2)
    seriestype := :shape

    # create a filled polygon for each item
    for c=1:size(y,2)
        sx = vcat(x, reverse(x))
        sy = vcat(y[:,c], c==1 ? zeros(n) : reverse(y[:,c-1]))
        @series (sx, sy)
    end
end

function main()
    if length(ARGS) != 1
        error("Need mem_events file")
    end

    events = readMemEvents(ARGS[1], memspace = "Cuda")

    totalSeries = getTotalSeries(events)
    allSeries = getAllNamedSeries(totalSeries)

    sNames = map(s -> s.name, allSeries)
    nSeries = length(sNames)
    series = mapreduce(s -> s.total, hcat, allSeries)

    # plotlyjs()
    stackedarea(
        totalSeries.time,
        series,
        labels=reshape(sNames, (1,nSeries)),
        legend=:outertopleft,
        scalefonts=0.5,
        wsize=(1200,400)
    )
    gui()
    # plot(totalSeries.time, totalSeries.total, show = true)
    readline()
end

main()
