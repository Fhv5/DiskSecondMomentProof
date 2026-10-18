using DiskSecondMomentProof

function option_value(arguments, name)
    index = findfirst(==(name), arguments)
    index === nothing && return nothing
    index < length(arguments) || throw(ArgumentError("missing value after $name"))
    return arguments[index + 1]
end

function main(arguments)
    length(arguments) % 2 == 0 || throw(ArgumentError("options require values"))
    all(arguments[index] in ("--n", "--input", "--output")
        for index in 1:2:length(arguments)) ||
        throw(ArgumentError("unknown census filter option"))
    n_text = option_value(arguments, "--n")
    input_path = option_value(arguments, "--input")
    output_path = option_value(arguments, "--output")
    n_text === nothing && throw(ArgumentError("--n 2, --n 3, --n 4, --n 5 or --n 6 is required"))
    input_path === nothing && throw(ArgumentError("--input PATH is required"))
    n = parse(Int, n_text)
    n in (2, 3, 4, 5, 6) || throw(ArgumentError("the census supports n=2,3,4,5,6"))

    codes = filter(!isempty, strip.(readlines(input_path)))
    expected_planar = Dict(2 => 1, 3 => 2, 4 => 6, 5 => 20, 6 => 99)[n]
    length(codes) == expected_planar ||
        throw(ArgumentError("expected $expected_planar planar graph classes, got $(length(codes))"))
    if n == 4
        coverage = Census.verify_four_disk_planar_census(codes)
        println("Four-disk census: complete=$(coverage.connected_classes) " *
                "excluded=$(coverage.k4_graph6) Gram_rank=$(coverage.gram_rank)>$(coverage.ambient_dimension)")
    end
    selected = Census.filter_planar_graphs(codes, n)
    expected_outerplanar = Dict(2 => 1, 3 => 2, 4 => 5, 5 => 13, 6 => 46)[n]
    length(selected) == expected_outerplanar ||
        throw(ArgumentError("expected $expected_outerplanar outerplanar classes, got $(length(selected))"))
    Census.compare_with_stored_census(selected, n) ||
        throw(ArgumentError("filtered graph classes differ from the stored census"))

    if output_path !== nothing
        isfile(output_path) && throw(ArgumentError("output already exists: $output_path"))
        abspath(output_path) == abspath(connected_graph_path(n)) &&
            throw(ArgumentError("cannot overwrite the stored proof census"))
        open(output_path, "w") do io
            for code in selected
                println(io, code)
            end
        end
    end
    println("n=$n planar=$(length(codes)) outerplanar=$(length(selected)) stored_classes_match=true")
    output_path === nothing || println("Filtered graph6: $output_path")
end

main(ARGS)
