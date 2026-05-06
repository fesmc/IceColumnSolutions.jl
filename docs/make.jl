using Documenter, IceColumnSolutions

DocMeta.setdocmeta!(IceColumnSolutions, :DocTestSetup,
                    :(using IceColumnSolutions); recursive=true)

makedocs(
    sitename = "IceColumnSolutions.jl",
    modules  = [IceColumnSolutions],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://fesmc.github.io/IceColumnSolutions.jl",
        mathengine = Documenter.MathJax2(),
    ),
    pages = [
        "Home"   => "index.md",
        "Theory" => "theory.md",
        "Benchmark Experiments" => [
            "Overview"                        => "benchmarks/overview.md",
            "Exp 1 — Pure diffusion"          => "benchmarks/exp1.md",
            "Exp 2 — Advection"               => "benchmarks/exp2.md",
            "Exp 3 — Strain heating"          => "benchmarks/exp3.md",
            "Exp 4 — Full case"               => "benchmarks/exp4.md",
        ],
        "Model Verification" => "verification.md",
        "API Reference"      => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo       = "github.com/fesmc/IceColumnSolutions.jl.git",
    devbranch  = "main",
    push_preview = true,
)
