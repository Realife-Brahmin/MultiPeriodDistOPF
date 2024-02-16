
function cycles_finding(graph,numCycles,cycles)
    for i in 1:size(graph,1)
        for j in 1:2
            (graph,numCycles,cycles) = findNewCycles(graph[i,j],graph,numCycles,cycles)
        end
    end
    cycleList = cycles
return cycleList
end

function findNewCycles(path,graph,numCycles,cycles)

startNode = path[1]
nextNode = NaN
sub = []

# visit each edge and each node of each edge
    for i in 1:size(graph,1)
        node1 = graph[i,1]
        node2 = graph[i,2]
        if startNode == graph[i,1] || startNode == graph[i,2]
            if node1 == startNode
                nextNode = node2
            elseif node2 == startNode
                nextNode = node1
            end

            if ~(visited(nextNode,path))
                # neighbor node not on path yet
                sub = nextNode
                sub = [sub path]
                # explore extended path
                (graph,numCycles,cycles) = findNewCycles(sub,graph,numCycles,cycles)
            elseif size(path,2) > 2 && nextNode == path[end]
                # cycle found
                path_new = permutedims(path)
                p   = rotate_to_smallest(path_new)
                inv = invert(p)
                if isNew(p,cycles) && isNew(inv,cycles)
                    numCycles = numCycles + 1
                    # cycles[numCycles] = [p]
                    # push!([p],cycles)
                    push!(cycles,[p])
                end
            end
        end
    end
    return graph,numCycles,cycles
end

function invert(path)
    inv = rotate_to_smallest(path[end:-1:1])
return inv
end

# rotate cycle path such that it begins with the smallest node
function rotate_to_smallest(path)
    # path = permutedims(path)
    (value,idx)  = findmin(path)
    idx = idx[1]
    # new_path = [path[idx:end], path[1:idx-1]]
    new_path = vcat(path[idx:end],path[1:idx-1])
    return new_path
end


function isNew(path,cycles)
result = true
    for i in 1:size(cycles,1)
        if size(path,1) == size(cycles[i][1],1) && all(path == cycles[i][1])
            result = false
            break
        end
    end
return result
end


function visited(node,path)
    result = false
    if isnan(node) && any(isnan.(path))
        result = true
        return result
    end
    for i in 1:size(path,2)
        if node == path[i]
            result = true
            break
        end
    end
return result
end
