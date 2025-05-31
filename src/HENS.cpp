#include <iostream>
#include <cmath>
#include "HENS_math_core.h"
#include "HENS_matrix.h"

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
int main(){
    HENS::Solver mainSolver(3);
    auto vec = mainSolver.residualIntegral(0.0, 1.0);
    for(auto& col : vec){
        for(auto& number : col){
            std::cout << number << " ";
        }
        std::cout << std::endl;
    }
    mainSolver.solve();
    mainSolver.printMatrices();

    
    
}