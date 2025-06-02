#include <iostream>
#include <cmath>
#include "HENS_math_core.h"
#include "HENS_matrix.h"
#include <fstream>
#include <iomanip>

void testMatrix(){
    HENS::Matrix m(4,4,2.0);
    m.setValue(0, 0, 2);
    m.setValue(0, 1, 3);
    m.setValue(0, 2, 4);
    m.setValue(1, 0, 4);
    m.setValue(1, 1, 1);
    m.setValue(1, 2, 6);
    m.setValue(2, 0, 3);
    m.setValue(2, 1, 3);
    m.setValue(2, 2, 7);
    m.setValue(0, 3, 3.1);
    m.setValue(2, 3, 5.5);
    m.setValue(3, 2, 6.8);
    m.print();
    std::cout << "-------------------" << std::endl;

    auto mInverse = m.inverseMatrix();
    mInverse.print();
    std::cout << "-------------------" << std::endl;
    auto I = mInverse * m;
    I.print();
    std::cout << "-------------------" << std::endl;
}
double func(double x){
    return std::sin(PI * x);
}
int main(){
    std::ofstream outfile("output.txt");
    outfile << std::fixed << std::setprecision(10);
    for(int N = 1; N <=7; N +=2){
        HENS::Solver mainSolver(N, HENS::ResdualMethod::GALERKIN);
        auto X = mainSolver.solve();
        for(auto& x : X){
            outfile << x <<"\t";
        }
        outfile << std::endl;
    }
    for(int N = 1; N <=7; N +=2){
        HENS::Solver mainSolver(N, HENS::ResdualMethod::SUBDOMAIN);
        auto X = mainSolver.solve();
        for(auto& x : X){
            outfile << x <<"\t";
        }
        outfile << std::endl;
    }
    for(int N = 1; N <=7; N +=2){
        HENS::Solver mainSolver(N, HENS::ResdualMethod::COLLOCATION);
        auto X = mainSolver.solve();
        for(auto& x : X){
            outfile << x <<"\t";
        }
        outfile << std::endl;
    }

    
    



}