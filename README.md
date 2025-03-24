# Heat Conduction in a 2D Composite Rod

## Description
This project presents a numerical study of transient heat conduction in a two-dimensional rod composed of four different materials. The heat equation is solved using the implicit finite difference method (Backward Euler) with Gauss-Seidel iteration for numerical convergence.

The study includes:
- Definition of boundary conditions.
- Numerical method implementation.
- Analysis of the temperature distribution over time.
- Visualization of results through Gnuplot.

## Code Structure
- **`Heat_Conduction_2D.cpp`** — Main C++ code implementing the heat conduction model.
- **`Heat_Conduction_2D.h`** — Header file containing function declarations and constants.
- **`Exercise2_Results.txt`** — Temperature data at selected points over time.
- **`Exercise2_TemperatureMap_Data.txt`** — Temperature distribution data for visualization.
- **`Exercise2_TemperatureDistribution_Data.txt`** — Temperature values for each node in the mesh.
- **`TemperatureMap_Plot.plt`** — Gnuplot script for visualizing temperature maps and contours.
- **`Es1_Heat_Conduction_2D_AlessiGiada.pdf`** — Report presenting the results and analysis in detail.

## Results
The detailed results, including temperature evolution at selected points and temperature distribution maps, are presented in the report: **`Es1_Heat_Conduction_2D_AlessiGiada.pdf`**.

## How to Run
1. Ensure **Gnuplot** is installed on your system.
2. Compile the code using a C++ compiler (e.g., `g++ Heat_Conduction_2D.cpp -o HeatConduction`).
3. Run the executable: `./HeatConduction`.
4. To visualize results, run the Gnuplot script: `gnuplot TemperatureMap_Plot.plt`.

## Author
**Giada Alessi**  
Master in Thermal Engineering  
Universitat Politècnica de Catalunya



