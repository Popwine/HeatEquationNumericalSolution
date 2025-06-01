#ifndef HNE_MATH_CORE_H
#define HNE_MATH_CORE_H

#include <vector>
#include "HENS_matrix.h"
#define PI 3.141592653589793238462643

namespace HENS{

class Solver{
    
private:
    size_t N_;
    Matrix M_;
    Matrix K_;
    Matrix F_;
public:
    
    Solver(size_t number);
    size_t getN() const;
    std::vector<std::vector<double>> residualIntegral(const double x1, const double x2) const;
    void solve();
    void printMatrices();
    /**
     * @brief Runge-Kutta
     */
    std::vector<double> RK4(
        std::vector<double>& XInitial,
        std::vector<std::function<double(std::vector<double>)>>& funcs,
        double stepSize,
        double endTime
    );

};

void printVec(std::vector<double> vec);

}





#endif // HNE_MATH_CORE_H