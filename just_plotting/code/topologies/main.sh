cd ./flow_analysis/
julia main_flow.jl
cd ../MadrasSokal/
julia "scripts/table.jl"
if command -v latex > /dev/null 2>&1; then
    echo "LaTeX is installed."
    julia "scripts/plaquette.jl"
    julia "scripts/topological_charge.jl"  
else
    echo "LaTeX is not installed."  
fi
cd ..

