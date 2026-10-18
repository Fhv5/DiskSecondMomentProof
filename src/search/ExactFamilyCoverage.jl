module ExactFamilyCoverage

using ..DiskSecondMomentProof: ContactGraph
import ..DiskSecondMomentProof

export verify_edg_coverage, verify_family_coverage

const Q = Rational{BigInt}

"""Dispatch only to a coverage argument tied to the declared family model."""
function verify_family_coverage(entry::Dict{String,Any})
    model = string(entry["family"]["model"])
    model == "four_disk_rhombus" &&
        return DiskSecondMomentProof.verify_four_disk_rhombus_coverage(entry)
    return verify_edg_coverage(entry)
end

"""Check that the declared E`dg arc is precisely the analytically derived family."""
function verify_edg_coverage(entry::Dict{String,Any})
    graph = ContactGraph(string(entry["graph6"]))
    graph.n == 6 || throw(ArgumentError("the quadratic arc must have six vertices"))
    edges = Set((edge.src, edge.dst) for edge in graph.edges)
    edges == Set(((1, 2), (1, 5), (2, 6), (3, 4), (3, 6), (4, 5), (5, 6))) ||
        throw(ArgumentError("the declared arc has the wrong contact graph"))
    family = entry["family"]::Dict{String,Any}
    string(family["model"]) == "quadratic_arc" ||
        throw(ArgumentError("the E`dg coverage proof requires the quadratic arc model"))
    domain = family["domain"]::Dict{String,Any}
    string(domain["relation"]) == "s^2=4-a^2" ||
        throw(ArgumentError("the E`dg arc relation changed"))
    interval = domain["a_interval"]::Dict{String,Any}
    DiskSecondMomentProof._parse_rational(string(interval["lower"])) == -1//1 &&
        DiskSecondMomentProof._parse_rational(string(interval["upper"])) == 1//1 &&
        interval["lower_closed"] === false && interval["upper_closed"] === false ||
        throw(ArgumentError("the E`dg domain must be exactly -1<a<1"))

    a = DiskSecondMomentProof._arc(Q[0//1, 1//1], Q[0//1])
    s = DiskSecondMomentProof._arc(Q[0//1], Q[1//1])
    one_arc = DiskSecondMomentProof._arc(1)
    zero_arc = DiskSecondMomentProof._arc(0)
    expected = [
        (-one_arc - a, s),
        ( one_arc - a, s),
        ( one_arc + a, -s),
        (-one_arc + a, -s),
        (-one_arc, zero_arc),
        ( one_arc, zero_arc),
    ]
    centers = DiskSecondMomentProof._parse_arc_centers(family["centers"], 6)
    all(centers[1, vertex] == expected[vertex][1] &&
        centers[2, vertex] == expected[vertex][2] for vertex in 1:6) ||
        throw(ArgumentError("the supplied arc does not equal the exhaustive critical family"))
    return :complete_critical_family
end

end
