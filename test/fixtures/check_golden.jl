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
include(joinpath(@__DIR__, "golden_json.jl"))

# ---------------------------------------------------------------- 重建算例
# Phase 0（旧实现自洽）：dump 与重跑同路径，无需重建——只验 JSON 健全性。
# Phase 2+：rebuild_case(name) 用换轨后的 API 重建算例并断言（见各 Phase 验收）。

# ---------------------------------------------------------------- 断言
function check(label, ok, detail="")
    println(ok ? "PASS  " : "FAIL  ", label, detail == "" ? "" : "  [$detail]")
    ok || global failures += 1
end

global failures = 0

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
