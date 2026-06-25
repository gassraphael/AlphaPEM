"""
    make.jl

Documenter.jl configuration for AlphaPEM documentation.
This script generates the static HTML documentation and prepares it for GitHub Pages deployment.

The documentation is organized hierarchically with the following structure:
- Home page with overview and quick links
- Getting Started section with installation and quick start guides
- User Guide for web interface and CLI usage
- Advanced section for calibration and parameter analysis
- API Reference for programmatic usage
- About section including roadmap, contributing guidelines, and publications
"""

using Documenter
using AlphaPEM

# Configure doctests to use AlphaPEM environment
DocMeta.setdocmeta!(
    AlphaPEM,
    :DocTestSetup,
    :(using AlphaPEM);
    recursive=true
)

# Build documentation with Documenter.jl
makedocs(
    # Project metadata
    sitename="AlphaPEM Documentation",
    authors="Raphaël Gass",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        size_threshold=nothing,
        assets=String[],
    ),

    # Documentation structure and pages
    pages=[
        "Home" => "index.md",
        "Getting Started" => [
            "Installation" => "getting_started/installation.md",
            "Quick Start" => "getting_started/quickstart.md",
        ],
        "User Guide" => [
            "Web Interface" => "user_guide/web_interface.md",
            "Command Line Usage" => "user_guide/cli_usage.md",
        ],
        "Advanced" => [
            "Parameter Calibration" => "advanced/calibration.md",
            "Valid Parameter Region Analysis" => "advanced/validity_analysis.md",
        ],
        "About" => [
            "Roadmap" => "about/roadmap.md",
            "Publications" => "about/publications.md",
        ],
    ],

    # Documentation source and build directories
    source="src",
    build="build",

    # Build options
    clean=true,
    doctest=false,
    linkcheck=false,
)

# Deploy documentation to GitHub Pages
deploydocs(
    repo="github.com/gassraphael/AlphaPEM.git",
    devbranch="main",
    branch="gh-pages",
    push_preview=false,
)
