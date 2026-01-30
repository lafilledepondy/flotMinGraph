# Min-Cost Flow Solver

This C++ program implements a successive shortest augmenting path algorithm (with Bellman-Ford potentials) to compute minimum-cost flows on capacitated directed graphs. Instance files follow the DIMACS `p min` format with `n` supply lines and `a` arc lines; ready-to-use samples live under `instances/`.

## Build & Run

1. Compile: `g++ -std=c++17 -O2 main_FlotCoutMin.cpp -o min_cost_flow`
2. Execute: `./min_cost_flow`

By default `main()` runs an internal regression suite (`checkCoutMin()`) that loads every sample instance and validates both optimal cost and flow feasibility. To inspect a single instance, uncomment the `printSolution(filename);` call in `main_FlotCoutMin.cpp` and point `filename` to the desired `.dat` file.

## Output

Successful runs report whether the resulting flow is feasible and print the optimal cost. Additional arc-level flow details can be shown by uncommenting the block in `printSolution()`.
