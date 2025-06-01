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

    // Find the inverse matrix of M_
    auto M_inversed = M_.inverseMatrix();
    auto M_inversed_multiply_F = M_inversed * F_;
    auto M_inversed_multiply_K = M_inversed * K_;

    // std::cout << "M_inversed_multiply_F" << std::endl;
    // M_inversed_multiply_F.print();
    // std::cout << "M_inversed_multiply_K" << std::endl;
    // M_inversed_multiply_K.print();

    std::vector<std::function<double(std::vector<double>)>> funcs;
    funcs.push_back([](std::vector<double> X){
        if(X[0] == 0){

        }
        return 1.0;
    });
    // get funtions for runge-kutta, 
    // func[0] stands for Time derivative of time, which is 1 
    for(size_t i = 0; i < getN(); i++){
        funcs.push_back([
            &M_inversed_multiply_F,
            &M_inversed_multiply_K,
            i
        ](std::vector<double> X){
            if(X.size() != M_inversed_multiply_F.getRowNumber() + 1){
                throw std::runtime_error("X doesn't match M^-1*F.");
            }
            if(X.size() != M_inversed_multiply_K.getRowNumber() + 1){
                throw std::runtime_error("X doesn't match M^-1*K.");
            }
            double Row_i_M_inversed_multiply_K_multiply_a = 0.0;
            for(size_t j = 0; j < M_inversed_multiply_K.getRowNumber(); j++){
                //std::cout << j <<std::endl;
                Row_i_M_inversed_multiply_K_multiply_a += 
                M_inversed_multiply_K.getValue(i, j) * X[j + 1];
            }
            //X = M^-1*F - M^-1*K*a
            return (
                M_inversed_multiply_F.getValue(i, 0) 
                - Row_i_M_inversed_multiply_K_multiply_a
            );
        });
    }
    
    // set original value of X
    // XInitial[0] stand for time
    // XInitial[1] to X[N] stand for a_1 to a_N
    // All of X's initial value is 0.0

    std::vector<double> XInitial(getN() + 1, 0.0);
    
    auto XEnd = RK4(XInitial, funcs, 0.2*0.01, 0.2);



}

void Solver::printMatrices(){
    std::cout<< "M:" << std::endl;
    M_.print();
    std::cout<< "K:" << std::endl;
    K_.print();
    std::cout<< "F:" << std::endl;
    F_.print();

}



std::vector<double> Solver::RK4(
    std::vector<double>& XInitial,
    std::vector<std::function<double(std::vector<double>)>>& funcs,
    double stepSize,
    double endTime
){
    std::vector<double> X = XInitial;
    while(X[0] < endTime){
        // k1 = f(x)
        std::vector<double> K1(X.size());
        for(size_t i = 0; i < X.size(); i++){
            K1[i] = funcs[i](X);
        }
        // k2 = f(x + stepSize / 2.0 * k1)
        std::vector<double> K2(X.size());
        for(size_t i = 0; i < X.size(); i++){
            std::vector<double> XInput(X.size());
            for(size_t j = 0; j < X.size(); j++){
                XInput[j] = X[j] + stepSize / 2.0 * K1[j];
            }
            K2[i] = funcs[i](XInput);
        }
        // k3 = f(x + stepSize / 2.0 * k2)
        std::vector<double> K3(X.size());
        for(size_t i = 0; i < X.size(); i++){
            std::vector<double> XInput(X.size());
            for(size_t j = 0; j < X.size(); j++){
                XInput[j] = X[j] + stepSize / 2.0 * K2[j];
            }
            K3[i] = funcs[i](XInput);
        }
        // k4 = f(x + stepSize * k3)
        std::vector<double> K4(X.size());
        for(size_t i = 0; i < X.size(); i++){
            std::vector<double> XInput(X.size());
            for(size_t j = 0; j < X.size(); j++){
                XInput[j] = X[j] + stepSize * K3[j];
            }
            K4[i] = funcs[i](XInput);
        }
        // x_n+1 = x_n + 1/6*(k1 + 2k2 + 2k3 + k4)
        for(size_t i = 0; i < X.size(); i++){
            X[i] = X[i] + 1.0 / 6.0 * stepSize * (K1[i] + 2.0 * K2[i] + 2.0 * K3[i] + K4[i]);
            
        }
        std::cout << "K1:" << std::endl;
        printVec(K1);
        std::cout << "K2:" << std::endl;
        printVec(K2);
        std::cout << "K3:" << std::endl;
        printVec(K3);
        std::cout << "K4:" << std::endl;
        printVec(K4);
        std::cout << "X:" << std::endl;
        printVec(X);

    }
    return X;
}

void printVec(std::vector<double> vec){
    for(auto& v : vec){
        std::cout << v << " ";
    }
    std::cout << std::endl;
}

}
