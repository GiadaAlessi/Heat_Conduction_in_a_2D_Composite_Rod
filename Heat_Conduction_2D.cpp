// 2-D Heat Conduction Problem in a Composite Rod

#include "Heat_Conduction_2D.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace std::chrono;

// Geometry
double const xM1 = 0.5, yM1 = 0.4;
double const xM2 = 0.6, yM2 = 0.7;
double const xM3 = 0.5, yM3 = 0.4;
double const xM4 = 0.6, yM4 = 0.1;

// Target Points
double xP1 = 0.65;
double yP1 = 0.56;
double xP2 = 0.74;
double yP2 = 0.72;

// Boundary Conditions
double const Tf = 33 + 273.15, alpha = 9, Tb = 23 + 273.15, Qf = 60;

// Thermo-Physical Properties
double const rho1 = 1500, rho2 = 1600, rho3 = 1900, rho4 = 2500;
double const cp1 = 750, cp2 = 770, cp3 = 810, cp4 = 930;
double const lambda1 = 170, lambda2 = 140, lambda3 = 200, lambda4 = 140;

// Numerical Data
const int Nx = 55, Ny = 40;
double Dx = (0.5 + 0.6) / Nx, Dy = 0.8 / Ny;
double dt = 1.0;
double t_start = 0.0;
double Interval = 2306;
int maxIte = 1e6;
double maxRes = 1e-6;

// Vectors Definition
int Mat[Ny][Nx]; // Material matrix
double T[Ny][Nx], T1[Ny][Nx], Tg[Ny][Nx]; // Temperature matrix
double aP[Ny][Nx], aE[Ny][Nx], aW[Ny][Nx],
aN[Ny][Nx], aS[Ny][Nx], bP[Ny][Nx];  // Coefficients

// Function to Compute Harmonic Mean
double HarmonicMean(double lambdaA, double lambdaB) {
    return 2.0 / ((1 / lambdaA) + (1 / lambdaB));
}

// Function to Fill the Material Matrix and Initialize the Temperature
void MaterialMatrix() {
    int jM1 = static_cast<int>(round(xM1 / Dx));
    int iM3 = static_cast<int>(round(yM3 / Dy));
    int iM4 = static_cast<int>(round(yM4 / Dy));

    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            if (i <= iM3 && j < jM1)
                Mat[i][j] = 3;  // M3
            else if (i > iM3 && j < jM1)
                Mat[i][j] = 1;  // M1
            else if (i < iM4 && j >= jM1)
                Mat[i][j] = 4;  // M4
            else
                Mat[i][j] = 2;  // M2

            cout << Mat[i][j] << " ";
        }
        cout << endl;
    }
}


// Function to Evaluate Coefficients
void InternalNodesCoefficients() {
    for (int i = 1; i < Ny-1; i++) {
        for (int j = 1; j < Nx-1; j++) {
            double lambda_mat[5] = {0, lambda1, lambda2, lambda3, lambda4};
            double rho_cp[5] = {0, rho1 * cp1, rho2 * cp2, rho3 * cp3, rho4 * cp4};

            double rhocpP = rho_cp[Mat[i][j]];
            double lambdaP = lambda_mat[Mat[i][j]];

            double lambdaE = HarmonicMean(lambdaP, lambda_mat[Mat[i][j+1]]);
            double lambdaW = HarmonicMean(lambdaP, lambda_mat[Mat[i][j-1]]);
            double lambdaN = HarmonicMean(lambdaP, lambda_mat[Mat[i-1][j]]);
            double lambdaS = HarmonicMean(lambdaP, lambda_mat[Mat[i+1][j]]);

            aE[i][j] = lambdaE / (Dx * Dx);
            aW[i][j] = lambdaW / (Dx * Dx);
            aN[i][j] = lambdaN / (Dy * Dy);
            aS[i][j] = lambdaS / (Dy * Dy);
            aP[i][j] = rhocpP / dt + aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j];
        }
    }
}

// Function for Boundary Conditions and Source Term
void BoundaryConditions(double t) {
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            double lambda_mat[5] = {0, lambda1, lambda2, lambda3, lambda4};
            double rho_cp[5] = {0, rho1 * cp1, rho2 * cp2, rho3 * cp3, rho4 * cp4};

            double rhocpP = rho_cp[Mat[i][j]];
            double lambdaP = lambda_mat[Mat[i][j]];

            if (j == 0)
                aE[i][j] = lambdaP / (Dx * Dx),
                aW[i][j] = 0,
                aN[i][j] = lambdaP / (Dy * Dy),
                aS[i][j] = lambdaP / (Dy * Dy),
                aP[i][j] = rhocpP / dt + aE[i][j] + aN[i][j] + aS[i][j] + alpha / Dx,
                bP[i][j] = T[i][j] * rhocpP / dt + alpha * Tf / Dx;  // Convection at Left
            else if (j == Nx-1)
                aE[i][j] = 0,
                aW[i][j] = lambdaP / (Dx * Dx),
                aN[i][j] = lambdaP / (Dy * Dy),
                aS[i][j] = lambdaP / (Dy * Dy),
                aP[i][j] = rhocpP / dt + aW[i][j] + aN[i][j] + aS[i][j] + lambdaP / (Dx * Dx),
                bP[i][j] = T[i][j] * rhocpP / dt + (8 + 0.005 * t + 273.15) * lambdaP / (Dx * Dx);  // Variable T at Right
            else if (i == 0)
                aE[i][j] = lambdaP / (Dx * Dx),
               aW[i][j] = lambdaP / (Dx * Dx),
               aN[i][j] = 0,
               aS[i][j] = lambdaP / (Dy * Dy),
               aP[i][j] = rhocpP / dt + aE[i][j] + aW[i][j] + aS[i][j],
               bP[i][j] = T[i][j] * rhocpP / dt + Qf / Dy;  // Constant Flux at the Top
            else if (i == Ny-1)
                aE[i][j] = lambdaP / (Dx * Dx),
                aW[i][j] = lambdaP / (Dx * Dx),
                aN[i][j] = lambdaP / (Dy * Dy),
                aS[i][j] = 0,
                aP[i][j] = rhocpP / dt + aE[i][j] + aW[i][j] + aN[i][j] + lambdaP / (Dy * Dy),
                bP[i][j] = T[i][j] * rhocpP / dt + Tb * lambdaP / (Dy * Dy); // Constant T at the Base
            else
                bP[i][j] = T[i][j] * rhocpP / dt;  // No Heat Generation
        }
    }
}

// Gauss-Seidel Solver
void GaussSeidelSolver() {
    // Time Loop
    double res = maxRes + 1;  // condition to enter the loop
    int ite = 0;
    while (res > maxRes && ite < maxIte) {
        double maxDiff = 0.0;
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                T1[i][j] = (aE[i][j] * T1[i][j+1] + aW[i][j] * T1[i][j-1] +
                               aN[i][j] * T1[i-1][j] + aS[i][j] * T1[i+1][j] + bP[i][j]) / aP[i][j];
                double diff = fabs(T1[i][j] - Tg[i][j]);
                if (diff > maxDiff) {
                    maxDiff = diff;
                }
            }
        }
        res = maxDiff;
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                Tg[i][j] = T1[i][j];
            }
        }
        ite++;
    }
}

int main() {
    // Initial Temperature at t = 0;
    double T0 = 8 + 273.15;

    // Start Timer
    auto start = high_resolution_clock::now();

    // Find Indexes of Target Points
    int jP1 = static_cast<int>(round(xP1 / Dx));
    int iP1 = static_cast<int>(round((0.8 - yP1) / Dy));
    int jP2 = static_cast<int>(round(xP2 / Dx));
    int iP2 = static_cast<int>(round((0.8 - yP2) / Dy));

    // Initial Map
    for (auto& row : T) {
        for (auto& elem : row) {
            elem = T0;
        }
    }

    MaterialMatrix(); // Compute Matrix
    InternalNodesCoefficients(); // Compute Internal Nodes

    // File to Print Results
    ofstream TargetPoints("Exercise2_Results.txt");
    ofstream TemperatureMap("Exercise2_TemperatureMap_Data.txt");
    ofstream TemperatureDistribution("Exercise2_TemperatureDistribution_Data.txt");

    for (double t = 0.0; t <= Interval; t += dt) {

        BoundaryConditions(t); // Compute Boundary Nodes and Source Term
        GaussSeidelSolver(); // Evaluate Temperature Profile

        // Output results
        if (t == 0 or t == 1000 or t == 2000 or t == 3000 or t == 4000 or t == 5000) {
            TargetPoints << t << " " << T[iP1][jP1] - 273.15 << " " << T[iP2][jP2] - 273.15<< endl;
        }

        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                T[i][j] = T1[i][j];
            }
        }
    }

    // Temperature Map at t = Interval
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            TemperatureMap <<Dx/2 + j * Dx << " " << (Ny - i - 1) * Dy + Dy/2<< " " << T[i][j] - 273.15 << endl;
            TemperatureDistribution << T[i][j] - 273.15 << " ";
        }
        TemperatureDistribution << endl;
    }

    TargetPoints.close();
    TemperatureMap.close();
    TemperatureDistribution.close();

    //Gnuplot
    ofstream GnuPlot("TemperatureMap_Plot.plt");

    GnuPlot << "set terminal pngcairo size 1000,700 enhanced font 'Arial,12'\n";
    GnuPlot << "set output 'TemperatureMap_Plot.png'\n";

    GnuPlot << "set pm3d at b\n";
    GnuPlot << "set dgrid3d 55,40 splines\n";
    GnuPlot << "set cbrange [8:40]\n";
    GnuPlot << "set palette defined ("
            << "0 'blue', "
            << "0.2 'cyan', "
            << "0.4 'green', "
            << "0.5 'light-green', "
            << "0.6 'yellow', "
            << "0.75 'orange', "
            << "0.9 'red', "
            << "1 'dark-red')\n";
    GnuPlot << "set colorbox\n";

    GnuPlot << "set xlabel 'X (m)'\n";
    GnuPlot << "set ylabel 'Y (m)'\n";
    GnuPlot << "set title 'Temperature Distribution at t = 5000 s' font 'Arial,24'\n";
    GnuPlot << "set xrange [0:1.1]\n";

    GnuPlot << "set contour base\n";
    GnuPlot << "set cntrparam levels incremental 10,1,40\n";
    GnuPlot << "set cntrlabel start 5 interval 5 format '%2.0f°C' font 'Arial,10'\n";
    GnuPlot << "unset surface\n";
    GnuPlot << "set view map\n";

    GnuPlot << "set table 'contour_data.txt'\n";
    GnuPlot << "splot 'Exercise2_TemperatureMap_Data.txt' using 1:2:3 with lines\n";
    GnuPlot << "unset table\n";

    GnuPlot << "unset pm3d\n";
    GnuPlot << "set multiplot\n";

    GnuPlot << "plot 'Exercise2_TemperatureMap_Data.txt' using 1:2:3 with image notitle\n";

    GnuPlot << "replot 'contour_data.txt' with lines lc rgb 'black' lw 1 notitle\n";

    GnuPlot << "unset multiplot\n";

    GnuPlot.close();

    int result = system("TemperatureMap_Plot.plt");
    if (result != 0) {
        cerr << "Error during the execution of Gnuplot!" << endl;
        return 1;
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);  // Total Duration

    cout << "Transient simulation complete. Results saved to Exercise2_Results.txt" << endl;
    cout << "Total Execution Time: " << duration.count() << " seconds" << endl;

    return 0;
}