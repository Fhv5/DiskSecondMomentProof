module EvidenceArchive

using JSON3
using SHA

export write_run_record, sha256_file, sha256_text_lf, text_file_hashes

sha256_file(path::AbstractString) = bytes2hex(open(sha256, path))

const TEXT_EXTENSIONS = Set([".jl", ".py", ".ps1", ".sh", ".toml", ".md", ".json",
                             ".jsonl", ".g6", ".txt", ".tsv", ".csv", ".tex",
                             ".yaml", ".yml"])

"""Hash UTF-8 text after changing CRLF to LF; preserve all other bytes."""
function sha256_text_lf(path::AbstractString)
    bytes = read(path)
    text = String(copy(bytes))
    isvalid(text) || throw(ArgumentError("normalized source hash requires valid UTF-8 text"))
    return bytes2hex(sha256(codeunits(replace(text, "\r\n" => "\n"))))
end

function text_file_hashes(root::AbstractString)
    hashes = Dict{String,String}()
    for (directory, _, files) in walkdir(root), filename in files
        lowercase(splitext(filename)[2]) in TEXT_EXTENSIONS || continue
        path = joinpath(directory, filename)
        hashes[replace(relpath(path, root), '\\' => '/')] = sha256_text_lf(path)
    end
    return hashes
end

function file_hashes(root::AbstractString; exclude=Set{String}())
    hashes = Dict{String,String}()
    for (directory, _, files) in walkdir(root), filename in files
        path = joinpath(directory, filename)
        relative = replace(relpath(path, root), '\\' => '/')
        relative in exclude && continue
        hashes[relative] = sha256_file(path)
    end
    return hashes
end

"""Freeze parameters, version information, file inventory, and file digests."""
function write_run_record(output_directory::AbstractString, n::Int,
                          parameters, connected_census, search_results, analysis)
    isfile(joinpath(output_directory, "search", "conventions.json")) ||
        error("search conventions must be written before the run inventory")
    package_root = normpath(joinpath(@__DIR__, "..", ".."))
    revision = readchomp(`git -C $package_root rev-parse HEAD`)
    status = readchomp(`git -C $package_root status --porcelain`)
    sources = Dict{String,String}()
    normalized_sources = Dict{String,String}()
    for directory in ("src", "scripts", "test", "data", "docs")
        source_root = joinpath(package_root, directory)
        merge!(sources, Dict("$directory/$name" => digest
            for (name, digest) in file_hashes(source_root)))
        merge!(normalized_sources, Dict("$directory/$name" => digest
            for (name, digest) in text_file_hashes(source_root)))
    end
    for filename in ("README.md", ".gitattributes", ".gitignore")
        path = joinpath(package_root, filename)
        isfile(path) || continue
        sources[filename] = sha256_file(path)
        normalized_sources[filename] = sha256_text_lf(path)
    end
    manifest = joinpath(package_root, "Manifest.toml")
    record = Dict(
        "schema_version" => 3,
        "n" => n,
        "julia_version" => string(VERSION),
        "git_revision" => revision,
        "git_status_porcelain" => status,
        "project_sha256" => sha256_file(joinpath(package_root, "Project.toml")),
        "manifest_sha256" => isfile(manifest) ? sha256_file(manifest) : nothing,
        "source_sha256" => sources,
        "source_text_sha256_lf" => normalized_sources,
        "source_text_normalization" => "UTF-8; CRLF replaced by LF; all other bytes preserved",
        "project_text_sha256_lf" => sha256_text_lf(joinpath(package_root, "Project.toml")),
        "manifest_text_sha256_lf" => isfile(manifest) ? sha256_text_lf(manifest) : nothing,
        "parameters" => parameters,
        "search_conventions_file" => "search/conventions.json",
        "connected_census" => connected_census,
        "interval_graphs" => length(search_results),
        "leaf_graphs" => count(result -> result.chart.selected_leaf !== nothing, search_results),
        "leaf_excluded" => count(result -> result.chart.selected_leaf !== nothing &&
            result.status == :complete && isempty(result.candidates), search_results),
        "isolated_boxes" => sum((length(result.candidates) for result in search_results); init=0),
        "analysis_objects" => analysis.processed_candidates,
        "equality_family_objects" => count(result -> result.conclusion == :target_family_exactly_verified, analysis.results),
        "equality_endpoint_objects" => count(result -> result.conclusion == :target_endpoint_exactly_verified, analysis.results),
        "four_disk_scope" => n == 4 ? Dict(
            "raw_connected_classes" => 6, "excluded_graph6" => "C~",
            "k4_gram_rank" => 3, "ambient_dimension" => 2,
            "family_interior" => "1<a^2<3; a>0; b>0; a^2+b^2=4",
            "family_boundaries" => [1, 3],
            "global_and_family_coverage_argument" => "docs/proof.md",
        ) : nothing,
        "analysis_conclusions" => Dict(string(key) => value
            for (key, value) in analysis.conclusion_counts),
        "files_sha256" => file_hashes(output_directory;
            exclude=Set(["run.json"])),
    )
    path = joinpath(output_directory, "run.json")
    open(path, "w") do io
        JSON3.write(io, record)
        println(io)
    end
    return path
end

end
