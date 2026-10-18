module BoxRecords

export record_box, box_record_path

"""Name a graph transcript without depending on graph6 punctuation."""
function box_record_path(root::AbstractString, stage::AbstractString,
                         graph6::AbstractString; leaf::Union{Nothing,Int}=nothing)
    hex = uppercase(bytes2hex(codeunits(graph6)))
    suffix = leaf === nothing ? "" : "-leaf-$leaf"
    return joinpath(root, stage, "$hex$suffix.tsv")
end

"""Record one DFS box depth and disposition; its rational bounds are derived by replay."""
function record_box(io::IO, box, decision::Symbol)
    println(io, box.depth, '\t', decision)
end

end
