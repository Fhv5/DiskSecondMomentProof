using JSON3

"""Convert immutable JSON3 containers into ordinary Julia dictionaries and arrays."""
function _materialize_json(value)
    value isa JSON3.Object && return Dict{String,Any}(
        string(key) => _materialize_json(item) for (key, item) in pairs(value)
    )
    value isa JSON3.Array && return Any[_materialize_json(item) for item in value]
    return value
end

"""Read one UTF-8 JSON file into ordinary Julia containers."""
function read_json_file(filepath::AbstractString)
    return _materialize_json(JSON3.read(read(filepath, String)))
end
