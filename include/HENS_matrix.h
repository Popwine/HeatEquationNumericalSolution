#ifndef HENS_MATRIX_H
#define HENS_MATRIX_H
#include <vector>
#include <iostream>
#include <stdexcept>
#include <functional>

namespace HENS{
enum MatrixType{
    UNIT_MATRIX,
    ZERO_MATRIX
};


class Matrix{
private:
    std::vector<std::vector<double>> matrix_;
    const size_t rowNumber_;
    const size_t colNumber_;
public:
    Matrix(size_t rowNum, size_t colNum);
    Matrix(size_t rowNum, size_t colNum, double value);
    /**
     * @brief 生成N阶方阵，可选类型
     * @param size 方阵大小
     * @param type 方阵类型，UNIT_MATRIX：单位矩阵, ZERO_MATRIX：零矩阵
     */
    Matrix(size_t size, MatrixType type);
    size_t getRowNumber() const;
    size_t getColNumber() const;
    double getValue(size_t i, size_t j) const;
    void setValue(size_t i, size_t j, double value);
    void print() const;

    /**
    * @brief 把第i行的k倍加到第j行
    */
    void addRowMultiple(size_t i, size_t j, double k);
    Matrix operator*(const Matrix& rhs) const;
    /**
    * @brief 把第i行变为原来的k倍
    */
    void scaleRow(size_t i, double k);

    /**
    * @brief 交换两行
    */
    void swap(size_t i, size_t j);

    /**
    * @brief 求逆矩阵
    */
    Matrix inverseMatrix();


};





}

#endif //HENS_MATRIX_H