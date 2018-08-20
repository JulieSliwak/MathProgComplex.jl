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
            "SDPhierarchy/guide.md",
            "SDPhierarchy/implementation_details.md",
            "SDPhierarchy/relax_settings.md",
            # "SDPhierarchy/mathprinciple.md",
            "SDPhierarchy/state_of_implementation.md",
            ]
        ]
)