#ifndef HNE_MATH_CORE_H
#define HNE_MATH_CORE_H

#include <vector>
#include <functional>
#include "HENS_matrix.h"
#define PI 3.141592653589793238462643

namespace HENS{
enum ResdualMethod{
    GALERKIN,
    SUBDOMAIN,
    COLLOCATION

};
class Solver{
    
private:
    size_t N_;
    Matrix M_;
    Matrix K_;
    Matrix F_;
    ResdualMethod method_;
public:
    
    Solver(size_t number, ResdualMethod method);
    size_t getN() const;
    std::vector<std::vector<double>> residualIntegral(const int i) const;
    std::vector<double> solve();
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


double simpsonIntegral(
    std::function<double(double)> func,
    double x1, double x2,
    int N
);

}





#endif // HNE_MATH_CORE_H