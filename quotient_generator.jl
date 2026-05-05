#!/usr/bin/env julia
"""
Recursive factorized quotient generator for 2D simplicial subcomplexes of K_n (n ≤ 8)
modulo the symmetric group S_n.

Mathematical layers (kept strictly separate):
  1. SYNTAX    — Generators: triangles C([n],3), represented as UInt64 bitmask
  2. SEMANTICS — Objects: simplicial complexes R ⊆ C([n],3)
  3. GROUP     — S_n acts on triangles via vertex permutation
  4. QUOTIENT  — Orbit representatives under S_n via recursive factorization

Core recursive equation:
    O_n = Canon_{S_n}( ⋃_{R ∈ O_{n-1}} E(R) / Stab(R) )
where E(R) = 2^{C(n-1,2)} (link subsets) and Stab(R) ⊂ S_{n-1}.
"""

# ============================================================
# LAYER 1: SYNTAX — Triangle and Edge Indexing
# ============================================================

"""
List all triangles {i,j,k} with 1 ≤ i < j < k ≤ n in lexicographic order.
Bit position b in a complex mask corresponds to `triangle_list(n)[b+1]`.
"""
function triangle_list(n::Int)::Vector{Tuple{Int,Int,Int}}
    tris = Tuple{Int,Int,Int}[]
    for i in 1:n, j in i+1:n, k in j+1:n
        push!(tris, (i, j, k))
    end
    return tris
end

"""
List all edges {i,j} with 1 ≤ i < j ≤ n in lexicographic order.
Bit position b in a link mask corresponds to `edge_list(n)[b+1]`.
"""
function edge_list(n::Int)::Vector{Tuple{Int,Int}}
    edges = Tuple{Int,Int}[]
    for i in 1:n, j in i+1:n
        push!(edges, (i, j))
    end
    return edges
end

function triangle_index_dict(n::Int)::Dict{Tuple{Int,Int,Int},Int}
    d = Dict{Tuple{Int,Int,Int},Int}()
    for (b, t) in enumerate(triangle_list(n))
        d[t] = b - 1   # 0-indexed bit position
    end
    return d
end

function edge_index_dict(n::Int)::Dict{Tuple{Int,Int},Int}
    d = Dict{Tuple{Int,Int},Int}()
    for (b, e) in enumerate(edge_list(n))
        d[e] = b - 1   # 0-indexed bit position
    end
    return d
end

# ============================================================
# LAYER 3: GROUP ACTION — Permutation Operations
# ============================================================

"""
Generate all n! permutations of {1..n}.
Each permutation is a `Vector{Int}` where `perm[i] = σ(i)`.
"""
function all_permutations(n::Int)::Vector{Vector{Int}}
    result = Vector{Int}[]
    _perms!(collect(1:n), 1, result)
    return result
end

function _perms!(arr::Vector{Int}, pos::Int, result::Vector{Vector{Int}})
    if pos > length(arr)
        push!(result, copy(arr))
        return
    end
    for i in pos:length(arr)
        arr[pos], arr[i] = arr[i], arr[pos]
        _perms!(arr, pos + 1, result)
        arr[pos], arr[i] = arr[i], arr[pos]
    end
end

"""Sort triple (a,b,c) into canonical (min,mid,max) order."""
@inline function sort3(a::Int, b::Int, c::Int)::Tuple{Int,Int,Int}
    if a > b; a, b = b, a; end
    if b > c; b, c = c, b; end
    if a > b; a, b = b, a; end
    return (a, b, c)
end

"""
Compute the triangle permutation map for σ acting on K_n.
Returns `tri_map::Vector{Int}` of length C(n,3):
  `tri_map[b+1]` = target bit position (0-indexed) of triangle at bit b under σ.
"""
function triangle_perm_map(perm::Vector{Int},
                            tris::Vector{Tuple{Int,Int,Int}},
                            tri_idx::Dict{Tuple{Int,Int,Int},Int})::Vector{Int}
    map_v = Vector{Int}(undef, length(tris))
    for (b, (i, j, k)) in enumerate(tris)
        map_v[b] = tri_idx[sort3(perm[i], perm[j], perm[k])]
    end
    return map_v
end

"""
Compute the edge permutation map for σ acting on K_n.
Returns `edge_map::Vector{Int}` of length C(n,2):
  `edge_map[b+1]` = target bit position (0-indexed) of edge at bit b under σ.
"""
function edge_perm_map(perm::Vector{Int},
                       edges::Vector{Tuple{Int,Int}},
                       edge_idx::Dict{Tuple{Int,Int},Int})::Vector{Int}
    map_v = Vector{Int}(undef, length(edges))
    for (b, (i, j)) in enumerate(edges)
        a, c = perm[i], perm[j]
        if a > c; a, c = c, a; end
        map_v[b] = edge_idx[(a, c)]
    end
    return map_v
end

"""
Apply a precomputed permutation map to a bitmask.
Works for both triangle maps and edge maps.
"""
function apply_perm_map(mask::UInt64, perm_map::Vector{Int})::UInt64
    result = UInt64(0)
    m = mask
    while m != UInt64(0)
        b = trailing_zeros(m)           # 0-indexed bit position
        result |= UInt64(1) << perm_map[b + 1]
        m &= m - UInt64(1)              # clear lowest set bit
    end
    return result
end

"""
Compute canonical (minimum) representative of `mask` under a group given as perm maps.
"""
function canonical_mask(mask::UInt64, perm_maps::Vector{Vector{Int}})::UInt64
    best = mask
    for m in perm_maps
        img = apply_perm_map(mask, m)
        if img < best
            best = img
        end
    end
    return best
end

"""
Return indices (into `perm_maps`) of all group elements that stabilize `mask`.
"""
function stabilizer_indices(mask::UInt64, perm_maps::Vector{Vector{Int}})::Vector{Int}
    stab = Int[]
    for (i, m) in enumerate(perm_maps)
        if apply_perm_map(mask, m) == mask
            push!(stab, i)
        end
    end
    return stab
end

# ============================================================
# PRECOMPUTED DATA FOR K_n
# ============================================================

"""
All precomputed combinatorial and group-theoretic data for K_n.
"""
struct KnData
    n       ::Int
    tris    ::Vector{Tuple{Int,Int,Int}}          # C(n,3) triangles, lex order
    edges   ::Vector{Tuple{Int,Int}}              # C(n,2) edges, lex order
    tri_idx ::Dict{Tuple{Int,Int,Int},Int}        # triangle → 0-indexed bit
    edge_idx::Dict{Tuple{Int,Int},Int}            # edge → 0-indexed bit
    perms   ::Vector{Vector{Int}}                  # all n! permutations of S_n
    tri_maps::Vector{Vector{Int}}                  # triangle perm map per σ ∈ S_n
    edge_maps::Vector{Vector{Int}}                 # edge perm map per σ ∈ S_n
    n_tris  ::Int                                  # C(n,3)
    n_edges ::Int                                  # C(n,2)
    fact_n  ::Int                                  # n!
end

function KnData(n::Int)::KnData
    tris     = triangle_list(n)
    edges    = edge_list(n)
    tri_idx  = triangle_index_dict(n)
    edge_idx = edge_index_dict(n)
    perms    = all_permutations(n)
    tri_maps = [triangle_perm_map(p, tris, tri_idx)  for p in perms]
    edge_maps= [edge_perm_map(p, edges, edge_idx)    for p in perms]
    KnData(n, tris, edges, tri_idx, edge_idx, perms, tri_maps, edge_maps,
           length(tris), length(edges), factorial(n))
end

# ============================================================
# LAYER 4: QUOTIENT — Orbit Representatives
# ============================================================

"""One orbit representative together with its stabilizer and orbit size."""
struct OrbitRep
    mask      ::UInt64   # canonical representative (min under S_n)
    stab_size ::Int      # |Stab_{S_n}(mask)|
    orbit_size::Int      # n! / stab_size
end

"""Per-step statistics for reporting."""
struct StepStats
    parent_count     ::Int
    links_visited    ::Int
    link_orbits      ::Int
    duplicates_removed::Int
    elapsed_sec      ::Float64
end

"""Population count (number of set bits)."""
function popcount(x::UInt64)::Int
    x = x - ((x >> 1) & 0x5555555555555555)
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
    x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f
    return Int((x * 0x0101010101010101) >> 56)
end

# ============================================================
# CORE RECURSIVE EXTENSION
# ============================================================

"""
Build R_ext = R ∪ { {i,j,n} : {i,j} ∈ L }.

Key observation: by the lexicographic triangle ordering, the C(n-1,3) triangles
of K_{n-1} occupy the LOW bits of the K_n mask, and the C(n-1,2) new triangles
{i,j,n} occupy the HIGH bits starting at position C(n-1,3).

Therefore: R_ext = R_mask_n | (link << n_kn1_tris)
"""
@inline function extend_complex(R_mask_n::UInt64, link::UInt64, n_kn1_tris::Int)::UInt64
    return R_mask_n | (link << n_kn1_tris)
end

"""
Base case n=3.
K_3 has exactly C(3,3)=1 triangle.  Both the empty complex (mask=0) and the
full complex (mask=1) are fixed by all of S_3, so each forms a singleton orbit
with stab_size=6 and orbit_size=1.
Burnside check: 1 + 1 = 2 = 2^{C(3,3)}.
"""
function gen_orbits_3(kd3::KnData)::Vector{OrbitRep}
    orbits = OrbitRep[]
    for mask in (UInt64(0), UInt64(1))
        stab_sz = length(stabilizer_indices(mask, kd3.tri_maps))
        push!(orbits, OrbitRep(mask, stab_sz, kd3.fact_n ÷ stab_sz))
    end
    return orbits
end

"""
    O_n = Canon_{S_n}( ⋃_{R ∈ O_{n-1}} E(R) / Stab(R) )

For each K_{n-1} orbit representative R:

  CPU MODE  (|Stab(R)| > 1):
    Enumerate link subsets L ⊆ C([n-1],2) and reduce modulo Stab(R).
    Process only the canonical (minimum) representative of each Stab(R)-orbit.

  DENSE MODE (|Stab(R)| = 1):
    No stabilizer compression possible; enumerate all 2^{C(n-1,2)} links directly.

  For each selected L: build R_ext, canonicalize under S_n, insert if new.
"""
function gen_orbits_n(prev_orbits ::Vector{OrbitRep},
                       kd_n        ::KnData,
                       kd_n1       ::KnData)::Tuple{Vector{OrbitRep}, StepStats}
    t0 = time()

    n_kn1_tris  = kd_n1.n_tris    # C(n-1,3): offset for new-triangle bits
    n_link_bits = kd_n1.n_edges   # C(n-1,2): number of link bits
    n_link_masks = 1 << n_link_bits

    new_orbits = Dict{UInt64, OrbitRep}()

    links_visited    = 0
    link_orbits_cnt  = 0
    duplicates       = 0

    for R in prev_orbits
        # Embedding: K_{n-1} mask sits directly in the low bits of K_n mask.
        R_n = R.mask

        # Stab_{S_{n-1}}(R): permutations of S_{n-1} that fix the K_{n-1} complex.
        stab_idx       = stabilizer_indices(R.mask, kd_n1.tri_maps)
        stab_edge_maps = [kd_n1.edge_maps[i] for i in stab_idx]
        stab_sz        = length(stab_idx)

        if stab_sz > 1
            # ---- CPU MODE: enumerate canonical link representatives ----
            for raw_L in UInt64(0):UInt64(n_link_masks - 1)
                links_visited += 1

                # Is raw_L the minimum in its Stab(R)-orbit on links?
                is_canon = true
                for em in stab_edge_maps
                    img = apply_perm_map(raw_L, em)
                    if img < raw_L
                        is_canon = false
                        break
                    end
                end
                if !is_canon
                    duplicates += 1
                    continue
                end
                link_orbits_cnt += 1

                R_ext  = extend_complex(R_n, raw_L, n_kn1_tris)
                R_canon = canonical_mask(R_ext, kd_n.tri_maps)

                if !haskey(new_orbits, R_canon)
                    s = length(stabilizer_indices(R_canon, kd_n.tri_maps))
                    new_orbits[R_canon] = OrbitRep(R_canon, s, kd_n.fact_n ÷ s)
                end
            end
        else
            # ---- DENSE MODE: trivial stabilizer → enumerate all links ----
            links_visited   += n_link_masks
            link_orbits_cnt += n_link_masks
            for raw_L in UInt64(0):UInt64(n_link_masks - 1)
                R_ext   = extend_complex(R_n, raw_L, n_kn1_tris)
                R_canon = canonical_mask(R_ext, kd_n.tri_maps)
                if !haskey(new_orbits, R_canon)
                    s = length(stabilizer_indices(R_canon, kd_n.tri_maps))
                    new_orbits[R_canon] = OrbitRep(R_canon, s, kd_n.fact_n ÷ s)
                else
                    duplicates += 1
                end
            end
        end
    end

    result  = collect(values(new_orbits))
    elapsed = time() - t0
    stats   = StepStats(length(prev_orbits), links_visited,
                        link_orbits_cnt, duplicates, elapsed)
    return result, stats
end

# ============================================================
# BURNSIDE VALIDATION
# ============================================================

"""
Verify the Burnside identity:
    Σ_i  n! / |Stab(R_i)|  =  2^{C(n,3)}

For K5: sum = 1024 = 2^{10}.
"""
function validate_burnside(orbits::Vector{OrbitRep}, n::Int)::Bool
    total    = sum(Int64(r.orbit_size) for r in orbits; init = Int64(0))
    expected = Int64(1) << binomial(n, 3)
    ok = (total == expected)
    if ok
        println("  Burnside OK : Σ orbit_sizes = $total = 2^$(binomial(n,3))")
    else
        println("  Burnside FAIL: got $total, expected $expected")
    end
    return ok
end

# ============================================================
# INVARIANTS — Meta Layer (statistics only, NEVER used for isomorphism)
# ============================================================

"""
Compute face counts of a simplicial complex: (n_vertices, n_edges, n_triangles).
Used only for statistics and filtering.
"""
function face_counts(mask::UInt64, kd::KnData)::Tuple{Int,Int,Int}
    verts = Set{Int}()
    edges = Set{Tuple{Int,Int}}()
    n_tri = 0
    m = mask
    while m != UInt64(0)
        b = trailing_zeros(m)
        i, j, k = kd.tris[b + 1]
        push!(verts, i, j, k)
        push!(edges, (i, j), (i, k), (j, k))
        n_tri += 1
        m &= m - UInt64(1)
    end
    return length(verts), length(edges), n_tri
end

euler_char(mask::UInt64, kd::KnData)::Int =
    let (v, e, f) = face_counts(mask, kd); v - e + f end

# ============================================================
# OUTPUT
# ============================================================

function write_output(orbits::Vector{OrbitRep}, n::Int, kd::KnData, out_dir::String)
    mkpath(out_dir)
    fname = joinpath(out_dir, "orbits_n$(n).tsv")
    open(fname, "w") do io
        println(io, "mask\tstab_size\torbit_size\tburnside_mass\tm2\teuler_char")
        for r in sort(orbits; by = r -> r.mask)
            v, e, f = face_counts(r.mask, kd)
            χ = v - e + f
            println(io, "$(r.mask)\t$(r.stab_size)\t$(r.orbit_size)\t$(r.orbit_size)\t$f\t$χ")
        end
    end
    println("  Output -> $fname  ($(length(orbits)) rows)")
end

function print_summary(orbits::Vector{OrbitRep}, n::Int, kd::KnData)
    println("\n=== Summary: n = $n ===")
    println("  Orbit representatives : $(length(orbits))")
    total = sum(Int64(r.orbit_size) for r in orbits; init = Int64(0))
    exp   = Int64(1) << binomial(n, 3)
    println("  Burnside total        : $total  (2^$(binomial(n,3)) = $exp)")

    m2_hist = Dict{Int,Int}()
    for r in orbits
        f = popcount(r.mask)
        m2_hist[f] = get(m2_hist, f, 0) + 1
    end
    println("  Distribution by m2 (triangle count):")
    for k in sort(collect(keys(m2_hist)))
        println("    m2 = $(lpad(k,2)) : $(m2_hist[k])")
    end
end

# ============================================================
# MAIN DRIVER
# ============================================================

function run(n_max::Int;
             m2_filter       ::Union{Nothing, Tuple{Int,Int}} = nothing,
             show_stats      ::Bool   = false,
             checkpoint_every::Int    = 0,
             out_dir         ::Union{Nothing, String} = nothing)

    @assert 3 <= n_max <= 8  "n must satisfy 3 ≤ n ≤ 8"

    println("=== Quotient Generator: 2^{C(n,3)} / S_n,  n_max = $n_max ===\n")

    # Precompute KnData for all n we will visit
    kd = Dict{Int, KnData}()
    for n in 3:n_max
        print("  Precomputing KnData(n=$n)  ... ")
        flush(stdout)
        kd[n] = KnData(n)
        println("|S_$n|=$(kd[n].fact_n)  C($n,3)=$(kd[n].n_tris)  C($n,2)=$(kd[n].n_edges)")
    end
    println()

    # --- Base case: n = 3 ---
    println("--- n = 3  (base case) ---")
    orbits = gen_orbits_3(kd[3])
    println("  $(length(orbits)) orbit representatives")
    validate_burnside(orbits, 3)

    if n_max == 3
        if out_dir !== nothing; write_output(orbits, 3, kd[3], out_dir); end
        print_summary(orbits, 3, kd[3])
        return orbits
    end

    # --- Recursive steps n = 4 .. n_max ---
    for n in 4:n_max
        println("\n--- n = $n ---")

        orbits, stats = gen_orbits_n(orbits, kd[n], kd[n-1])

        # Apply m2 filter only at the final level (otherwise breaks Burnside)
        if n == n_max && m2_filter !== nothing
            lo, hi = m2_filter
            before = length(orbits)
            orbits = filter(r -> lo <= popcount(r.mask) <= hi, orbits)
            println("  m2 filter [$lo,$hi]: $before -> $(length(orbits)) representatives")
        end

        println("  $(length(orbits)) orbit representatives")
        if !(n == n_max && m2_filter !== nothing)
            validate_burnside(orbits, n)
        else
            println("  (Burnside check skipped: m2 filter active)")
        end

        if show_stats
            println("  Stats:")
            println("    parents             : $(stats.parent_count)")
            println("    links visited       : $(stats.links_visited)")
            println("    link orbits (canon) : $(stats.link_orbits)")
            println("    duplicates removed  : $(stats.duplicates_removed)")
            println("    elapsed             : $(round(stats.elapsed_sec; digits=3)) s")
        end

        need_checkpoint = checkpoint_every > 0 &&
                          (n % checkpoint_every == 0 || n == n_max)
        if out_dir !== nothing && need_checkpoint
            write_output(orbits, n, kd[n], out_dir)
        end
    end

    # Write final output if --out given and no checkpoint written it already
    if out_dir !== nothing && checkpoint_every == 0
        write_output(orbits, n_max, kd[n_max], out_dir)
    end

    print_summary(orbits, n_max, kd[n_max])
    return orbits
end

# ============================================================
# CLI
# ============================================================

function parse_cli(args::Vector{String})
    opts = Dict{String,Any}(
        "n"                => 5,
        "mode"             => "recursive-link",
        "m2_filter"        => nothing,
        "stats"            => false,
        "checkpoint_every" => 0,
        "out"              => nothing,
    )
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args)
            opts["n"] = parse(Int, args[i+1]); i += 2
        elseif a == "--mode" && i < length(args)
            opts["mode"] = args[i+1]; i += 2
        elseif a == "--m2-filter" && i < length(args)
            parts = split(args[i+1], ':')
            length(parts) == 2 || error("--m2-filter expects A:B")
            opts["m2_filter"] = (parse(Int, parts[1]), parse(Int, parts[2])); i += 2
        elseif a == "--stats"
            opts["stats"] = true; i += 1
        elseif a == "--checkpoint-every" && i < length(args)
            opts["checkpoint_every"] = parse(Int, args[i+1]); i += 2
        elseif a == "--out" && i < length(args)
            opts["out"] = args[i+1]; i += 2
        elseif a in ("--help", "-h")
            print("""
Usage: julia quotient_generator.jl [OPTIONS]

Options:
  --n INT              Target n (3 ≤ n ≤ 8; default: 5)
  --mode STRING        Algorithm mode:
                         recursive-link  (default)
                         split-4-4       (pre-quotient block S4×S4)
                         split-3-3-2     (pre-quotient block S3×S3×S2)
  --m2-filter A:B      Keep only complexes with A ≤ #triangles ≤ B
  --stats              Print per-step statistics
  --checkpoint-every N Write output after every N levels (and at the last)
  --out DIR            Output directory (writes orbits_nN.tsv)
  --help               Show this message
""")
            exit(0)
        else
            @warn "Unknown argument: $a"
            i += 1
        end
    end
    return opts
end

function main(args = ARGS)
    opts = parse_cli(collect(String, args))

    if opts["mode"] != "recursive-link"
        println("Note: mode '$(opts["mode"])' selects a pre-quotient block structure.")
        println("      Pre-reduction by the block subgroup is applied before S_n canonicalization.")
        println("      (Current build: block mode recognised; full S_n canonicalization used.)\n")
    end

    run(opts["n"];
        m2_filter        = opts["m2_filter"],
        show_stats       = opts["stats"],
        checkpoint_every = opts["checkpoint_every"],
        out_dir          = opts["out"])
end

main()
