# golden JSON 方言解析器（dump_golden.jl / dump_golden_weyl.jl 的极简写出器：
# 键不带引号、复数 [re,im]、非有限数编码为字符串——U1GTEMPO2 check_golden 同款）。

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


const GOLDEN_DIR = joinpath(@__DIR__, "port_golden")

# dump 方言里复数 = [re, im] 数组；取数处归一
ascomplex(v::Complex) = v
ascomplex(v::AbstractVector) = ComplexF64(v[1], v[2])
