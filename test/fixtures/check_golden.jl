# Phase 0 golden checker：对 dump 的 port_golden JSON 做解析与自洽校验。
#
# 用法（Phase 0，旧实现自洽）：
#   julia --project=. test/fixtures/check_golden.jl
# Phase 2+：同一文件在换轨后的实现上重跑（rebuild_case 换新 API 路径），
# 断言对象 = 物理量（occupation/gf 1e-10 档）；K/IF 稠密系数只 dump 不跨栈
# 断言（变量语义不同——vertex-pair vs Chrono cell，见 plan Phase 0 铁律）。
#
# JSON 方言：dump_golden.jl 的极简写出器（键不带引号、复数 [re,im]、
# 非有限数编码为字符串）——本文件自带同构解析器（U1GTEMPO2 check_golden 同款）。

using Test, LinearAlgebra
using GTEMPO
using GTEMPO: branches
using Z2Tensors: NoTruncation, truncdimcutoff
using ImpurityModelBase

const GOLDEN_DIR = joinpath(@__DIR__, "port_golden")

# ---------------------------------------------------------------- JSON 解析
function skipws(s, i)
    while i <= ncodeunits(s) && isspace(s[i]); i = nextind(s, i); end
    return i
end

function parse_value(s, i)
    i = skipws(s, i)
    c = s[i]
    if c == '{'
        i = skipws(s, i + 1); dict = Dict{String,Any}()
        if s[i] == '}'; return dict, nextind(s, i); end
        while true
            kstart = i; while s[i] != ':'; i = nextind(s, i); end
            key = strip(s[kstart:prevind(s, i)], '"')
            v, i = parse_value(s, nextind(s, i))
            dict[key] = v
            i = skipws(s, i)
            if s[i] == ','; i = skipws(s, nextind(s, i)); continue; end
            s[i] == '}' && return dict, nextind(s, i)
            error("bad object at $i")
        end
    elseif c == '['
        i = skipws(s, i + 1); arr = Any[]
        if s[i] == ']'; return arr, nextind(s, i); end
        while true
            v, i = parse_value(s, i)
            push!(arr, v)
            i = skipws(s, i)
            if s[i] == ','; i = skipws(s, nextind(s, i)); continue; end
            s[i] == ']' && return arr, nextind(s, i)
            error("bad array at $i")
        end
    elseif c == '"'
        j = nextind(s, i)
        while s[j] != '"'; j = nextind(s, j); end
        return s[nextind(s, i):prevind(s, j)], nextind(s, j)
    else
        j = i
        while j <= ncodeunits(s) && (s[j] != ',' && s[j] != ']' && s[j] != '}'); j = nextind(s, j); end
        tok = strip(s[i:prevind(s, j)])
        return parse_token(tok), j
    end
end

function parse_token(tok)
    if startswith(tok, "\"")
        return strip(tok, '"')
    elseif tok == "null"; return nothing
    elseif tok == "true"; return true
    elseif tok == "false"; return false
    elseif occursin(",", tok)  # 复数 [re,im]
        parts = strip.(split(strip(tok, ['[', ']']), ","))
        re = parts[1] in ("NaN", "Inf", "-Inf") ? parse_literal(parts[1]) : parse(Float64, parts[1])
        im = parts[2] in ("NaN", "Inf", "-Inf") ? parse_literal(parts[2]) : parse(Float64, parts[2])
        return ComplexF64(re, im)
    else
        return parse_literal(tok)
    end
end
parse_literal(t) = t == "NaN" ? NaN : t == "Inf" ? Inf : t == "-Inf" ? -Inf :
    occursin(".", t) || occursin("e", t) ? parse(Float64, t) : parse(Int, t)

function load_golden(name)
    path = joinpath(GOLDEN_DIR, "port_golden_$name.json")
    isfile(path) || return nothing
    s = read(path, String)
    v, _ = parse_value(s, 1)
    return v
end

golden_names() = [replace(f, r"^port_golden_|\.json$" => "")
                  for f in readdir(GOLDEN_DIR) if endswith(f, ".json")]

# ---------------------------------------------------------------- 重建算例
# Phase 0（旧实现自洽）：dump 与重跑同路径，无需重建——只验 JSON 健全性。
# Phase 2+：rebuild_case(name) 用换轨后的 API 重建算例并断言（见各 Phase 验收）。

# ---------------------------------------------------------------- 断言
function check(label, ok, detail="")
    println(ok ? "PASS  " : "FAIL  ", label, detail == "" ? "" : "  [$detail]")
    ok || global failures += 1
end

global failures = 0

# dump 方言里复数 = [re, im] 数组（parse_value 走数组分支）；取数处归一
ascomplex(v::Complex) = v
ascomplex(v::AbstractVector) = ComplexF64(v[1], v[2])

@testset "golden fixture 自洽（Phase 0）" begin
    names = golden_names()
    @test length(names) >= 10
    for name in names
        ref = load_golden(name)
        @test ref !== nothing
        @test haskey(ref, "Z")
        @test haskey(ref, "indexmap")
        @test haskey(ref, "corr_a_abar")
        # JSON 往返数值健全性
        z = ascomplex(ref["Z"])
        @test z isa Complex
        @test isfinite(abs(z)) || iszero(abs(z))
        # indexmap 是全量双射
        sites = [e["site"] for e in ref["indexmap"]]
        @test length(sites) == length(unique(sites)) == ref["L"]
    end
end

if failures == 0
    println("check_golden Phase 0: all green")
else
    println("check_golden Phase 0: $failures FAILURES")
    exit(1)
end
