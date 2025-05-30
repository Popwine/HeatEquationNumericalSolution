#include "HENS_math_core.h"
#include <cmath>
namespace HENS{

Solver::Solver(int number): N_(number){

}


int Solver::getN() const{
    return N_;
}


/**
 * @brief 返回残差的在区间x1到x2的积分
 * @return 二维vector vec[i][j]，i=0代表对时间的导数aj点的系数，i=1代表aj的系数，j=0时，不管i等于几，都代表常数项。
 */
std::vector<std::vector<double>> Solver::residualIntegral(const double x1, const double x2) const{

    std::vector<std::vector<double>> parameters(2, std::vector<double>(getN() + 1, 0.0));
    parameters[0][0] = parameters[1][0] = -PI * std::cos(PI * x2) - (-PI * std::cos(PI * x1));

    //处理
    for(int j = 1; j <= getN(); j++){
        parameters[0][j] = 1.0 / (j + 1) * std::pow(x2, j + 1.0) - 1.0 / (j + 2) * std::pow(x2, j + 2) - (1.0 / (j + 1.0) * std::pow(x1, j + 1.0) - 1.0 / (j + 2.0) * std::pow(x1, j + 2.0));
    }

    parameters[1][1] = 2.0 * x2 - 2.0 * x1;

    for(int j = 2; j <= getN(); j++){
        parameters[1][j] = -(j * std::pow(x2, j - 1.0) - (j + 1) * std::pow(x2, j)) + (j * std::pow(x1, j - 1) - (j + 1.0) * std::pow(x1, j));
    }
    return parameters;
}

}
