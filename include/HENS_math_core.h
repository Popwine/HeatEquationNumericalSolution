#ifndef HNE_MATH_CORE_H
#define HNE_MATH_CORE_H

#include <vector>
#define PI 3.141592653589793238462643

namespace HENS{

class Solver{
private:
    int N_;
public:
    Solver(int number);
    int getN() const;
    std::vector<std::vector<double>> residualIntegral(const double x1, const double x2) const;
};



}





#endif // HNE_MATH_CORE_H