#include "HENS_math_core.h"
#include <cmath>
namespace HENS{

Solver::Solver(size_t number): 
N_(number), 
M_(number, ZERO_MATRIX),
K_(number, ZERO_MATRIX),
F_(number, 1, 0.0)
{

}


size_t Solver::getN() const{
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
    for(size_t j = 1; j <= getN(); j++){
        parameters[0][j] = 1.0 / (j + 1) * std::pow(x2, j + 1.0) - 1.0 / (j + 2) * std::pow(x2, j + 2) - (1.0 / (j + 1.0) * std::pow(x1, j + 1.0) - 1.0 / (j + 2.0) * std::pow(x1, j + 2.0));
    }

    parameters[1][1] = 2.0 * x2 - 2.0 * x1;

    for(size_t j = 2; j <= getN(); j++){
        parameters[1][j] = -(j * std::pow(x2, j - 1.0) - (j + 1) * std::pow(x2, j)) + (j * std::pow(x1, j - 1) - (j + 1.0) * std::pow(x1, j));
    }
    return parameters;
}

void Solver::solve(){
    // Integrate in N child domains.
    for(size_t i = 0; i < getN(); i++){
        auto parameters = residualIntegral(
            static_cast<double>(i) / static_cast<double>(getN()),
            static_cast<double>(i + 1) / static_cast<double>(getN())
        );
        for(size_t j = 0; j < getN(); j++){
            M_.setValue(i, j, parameters[0][j + 1]);
            K_.setValue(i, j, parameters[1][j + 1]);
            
        }
        F_.setValue(i, 0, -parameters[0][0]);
    }

    //

}
void Solver::printMatrices(){
    std::cout<< "M:" << std::endl;
    M_.print();
    std::cout<< "K:" << std::endl;
    K_.print();
    std::cout<< "F:" << std::endl;
    F_.print();

}


}
