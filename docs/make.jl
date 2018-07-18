using Documenter, MathProgComplex, DataStructures

makedocs(
    modules = [MathProgComplex],
    format = :html,
    sitename = "MathProgComplex.jl",
    pages = [
        "Home" => "index.md",
        "Polynomial Optimization" => Any[
            "PolynomialOptim/polynomialoptim_structures.md"
        ],
        "SDP hierarchy" => Any[
            "SDPhierarchy/intro.md",
            "SDPhierarchy/mathprinciple.md",
            "SDPhierarchy/relax_settings.md",
            "SDPhierarchy/problems.md",
            "SDPhierarchy/solving.md"
            ]
        ]
)