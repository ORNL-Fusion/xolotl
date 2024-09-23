#!/usr/bin/julia

using Plots

mutable struct RegionLabel
    str::String
end

RegionLabel() = RegionLabel("")
Base.convert(::Type{RegionLabel}, s::AbstractString) = RegionLabel(s)

Base.string(rl::RegionLabel) = rl.str
Base.isempty(rl::RegionLabel) = isempty(rl.str)

function Base.push!(rl::RegionLabel, name::AbstractString)
    if isempty(rl.str)
        rl.str = name
    else
        rl.str = rl.str * " > " * name
    end
end

function Base.pop!(rl::RegionLabel)
    words = split(rl.str)
    if length(words) == 1
        rl.str = ""
    else
        rl.str = join(words[1:(end - 2)], " ")
    end
end

function Base.last(rl::RegionLabel)
    return split(rl.str)[end]
end

function popfront!(rl::RegionLabel)
    words = split(rl.str)
    if words[1] == ">"
        rl.str = " " * join(words[3:end], " ")
    else
        rl.str = " " * join(words[2:end], " ")
    end
end

struct MemEvent
    time::Float64
    ptr::String
    size::Int64
    space::String
    name::String
    region::String
end

function MemEvent(ser::Vector{<:AbstractString}, region::AbstractString)
    name = "no_name"
    if length(ser) >= 6
        name = join(ser[6:end], " ")
    end
    time = parse(Float64, ser[1])
    ptr = ser[2]
    size = parse(Int64, ser[3])
    space = ser[4]
    MemEvent(time, ptr, size, space, name, region)
end

mutable struct RegionStack
    stack::Vector{String}
    label::String
    indent::String
    newLine::Any

    RegionStack() = new(Vector{String}(), "", "", rs -> "\n" * rs.indent)
end

function setNewLine!(rs::RegionStack, newLine::Function)
    rs.newLine = newLine
end

function Base.push!(rs::RegionStack, name::AbstractString)
    push!(rs.stack, name)
    if isempty(rs.label)
        rs.label = name
    else
        rs.label = rs.label * rs.newLine(rs) * " > " * name
        rs.indent = rs.indent * "    "
    end
    return rs
end

function Base.pop!(rs::RegionStack)
    name = pop!(rs.stack)
    if isempty(rs.stack)
        rs.indent = ""
        rs.label = ""
    else
        rs.indent = chopsuffix(rs.indent, "    ")
        rs.label = chopsuffix(rs.label, rs.newLine(rs) * " > " * name)
    end
    return rs
end

Base.isempty(rs::RegionStack) = isempty(rs.stack)

function readMemEvents(filename::AbstractString; memspace = "", skip = 0, depth = 0)
    startLevel = 1 + skip
    ret = MemEvent[]
    lines = readlines(filename)
    sizehint!(ret, length(lines))
    rlabel = RegionLabel()
    regStack = Vector{String}()
    level = 1
    checkdepth = level -> true
    if depth > 0
        checkdepth = level -> (level - startLevel) < depth
    end
    for line in lines
        if startswith(line, "#")
            continue
        end
        if occursin("PushRegion", line)
            if (level >= startLevel) && checkdepth(level)
                region = split(line)[3]
                push!(rlabel, region)
            end
            push!(regStack, split(line)[3])
            level += 1
            continue
        end
        if occursin("PopRegion", line)
            if (level >= startLevel) && checkdepth(level)
                pop!(rlabel)
            end
            pop!(regStack)
            level -= 1
            continue
        end
        label = isempty(rlabel) ? "unlabeled" : string(last(rlabel))
        # if label == "unlabeled"
        #     @show rlabel.str, label
        #     @show line
        #     @show regStack
        #     readline()
        # end
        event = MemEvent(split(line), label)
        if isempty(memspace) || occursin(event.space, memspace)
            push!(ret, event)
        end
    end
    return ret
end

# function readMemEvents(filename::AbstractString; memspace = "")
#     ret = MemEvent[]
#     lines = readlines(filename)
#     sizehint!(ret, length(lines))
#     label = "unlabeled"
#     regStack = Vector{String}()
#     # regStack = RegionStack()
#     # setNewLine!(regStack, rs -> "")
#     for line in lines
#         if startswith(line, "#")
#             continue
#         end
#         find_f = occursin(line)
#         if find_f("PushRegion")
#             region = split(line)[3]
#             push!(regStack, region)
#             continue
#         end
#         if find_f("PopRegion")
#             pop!(regStack)
#             continue
#         end
#         label = isempty(regStack) ? "unlabeled" : last(regStack)
#         event = MemEvent(split(line), label)
#         if isempty(memspace) || occursin(event.space, memspace)
#             push!(ret, event)
#         end
#     end
#     return ret
# end

mutable struct NamedSeries
    name::String
    events::Vector{MemEvent}
    time::Vector{Float64}
    total::Vector{Float64}
    curEvtIdx::Int64
end

function NamedSeries(name::String, events::Vector{MemEvent}, size::Int64)
    ret = NamedSeries(name, events, [0.0, 0.0], [0.0, 0.0], 1)
    sizehint!(ret.time, size + 1)
    sizehint!(ret.total, size + 1)
    return ret
end

NamedSeries(name::String, size::Int64) = NamedSeries(name, MemEvent[], size)

@inline Base.isequal(a::NamedSeries, b::NamedSeries) = isequal(a.name, b.name)

function getTotalSeries(events::Vector{MemEvent})
    ret = NamedSeries("Total", events, [0.0, 0.0], [0.0, 0.0], 1)
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
        ser = get!(map, e.region, NamedSeries(e.region, numEvents))
        if findfirst(isequal(ser), ret) == nothing
            push!(ret, ser)
        end
        push!(ser.events, e)
        for ei = (ser.curEvtIdx + 1):(i + 1)
            push!(ser.time, totalSeries.time[ei])
            push!(ser.time, totalSeries.time[ei])
        end
        append!(ser.total, repeat([last(ser.total)], (i + 1 - ser.curEvtIdx) * 2 - 1))
        push!(ser.total, last(ser.total) + e.size)
        ser.curEvtIdx = i+1
    end

    for ser in ret
        i = numEvents + 1
        for ei = (ser.curEvtIdx + 1):numEvents+1
            push!(ser.time, totalSeries.time[ei])
            push!(ser.time, totalSeries.time[ei])
        end
        append!(ser.total, repeat([last(ser.total)], (i - ser.curEvtIdx) * 2))
        if last(ser.total) < 0.0
            println(ser.name, ": ", last(ser.total))
            # for e in events
            #     if e.region != ser.name || e.size > 0
            #         continue
            #     end
            #     println("\t", e.name, ": ", e.size)
            # end
        end
    end

    return ret
end

function extractInstances(allSeries::Vector{NamedSeries}) end

@userplot StackedArea

# a simple "recipe" for Plots.jl to get stacked area plots
# usage: stackedarea(xvector, datamatrix, plotsoptions)
@recipe function f(pc::StackedArea)
    x, y = pc.args
    n = length(x)
    y = cumsum(y, dims = 2)
    seriestype := :shape
    seriescolor := :auto

    # create a filled polygon for each item
    for c = 1:size(y, 2)
        sx = vcat(x, reverse(x))
        sy = vcat(y[:, c], c == 1 ? zeros(n) : reverse(y[:, c - 1]))
        @series (sx, sy)
    end
end

function main()
    if length(ARGS) < 1
        error("Need mem_events file")
    end

    memspace = "Cuda"
    if length(ARGS) == 2
        memspace = ARGS[2]
    end

    events = readMemEvents(ARGS[1], memspace = memspace, skip = 0)
    # events = readMemEvents(ARGS[1], memspace = "Cuda")

    totalSeries = getTotalSeries(events)
    println("Max CUDA allocated: ", maximum(totalSeries.total))
    allSeries = getAllNamedSeries(totalSeries)

    sNames = map(s -> s.name, allSeries)
    nSeries = length(sNames)
    series = mapreduce(s -> s.total, hcat, allSeries)

    stackedarea(
        # totalSeries.time,
        allSeries[1].time,
        series,
        labels = reshape(sNames, (1, nSeries)),
        # legend = :outerleft,
        legend=:topleft,
        # scalefonts=0.5,
        # thickness_scaling=0.75,
        # wsize=(1200,400)
    )
    gui()
    # plot(allSeries[1].time, series, labels = reshape(sNames, (1, nSeries)), show = true)
    # plot(totalSeries.time, totalSeries.total, show = true,
    #     scalefonts=2)
    print("To finish, press ENTER: ")
    readline()
end

pyplot()
main()
