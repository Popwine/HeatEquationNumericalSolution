#include "HENS_matrix.h"
namespace HENS{

Matrix::Matrix(const size_t rowNum, const size_t colNum)
: matrix_(rowNum, std::vector<double>(colNum, 0.0)), rowNumber_(rowNum), colNumber_(colNum)
{
    
}

Matrix::Matrix(const size_t rowNum, const size_t colNum, double value)
: matrix_(rowNum, std::vector<double>(colNum, value)), rowNumber_(rowNum), colNumber_(colNum)
{
    
}

Matrix::Matrix(size_t size, MatrixType type)
: matrix_(size, std::vector<double>(size, 0.0)), rowNumber_(size), colNumber_(size)
{
    if(type == UNIT_MATRIX){
        for(size_t i = 0; i < size; i++){
            setValue(i, i, 1.0);
        }
    }
    if(type == ZERO_MATRIX){
        /*Doing Nothing.*/
    }


}

size_t Matrix::getRowNumber() const{
    return rowNumber_;
}

size_t Matrix::getColNumber() const{
    return colNumber_;
}

double Matrix::getValue(size_t i, size_t j) const{
    
    if (i >= getRowNumber() || j >= getColNumber()) {
        
        throw std::out_of_range("Row index out of range in getValue.\n");
    }
    return matrix_[i][j];
}

void Matrix::setValue(size_t i, size_t j, double value) {
    //std::cout << i << getRowNumber() << j << getColNumber() << std::endl;
    if (i >= getRowNumber() || j >= getColNumber()) {
        throw std::out_of_range("Row index out of range in setValue");
    }
    matrix_[i][j] = value;
}

void Matrix::print() const{
    for(const auto& row : matrix_){
        for(const auto& number : row){
            std::cout << number << " ";
        }
        std::cout << std::endl;
    }
}

void Matrix::addRowMultiple(size_t i, size_t j, double k){
    if (i >= getRowNumber() || j >= getRowNumber()) {
        throw std::out_of_range("Row index out of range in addRowMultiple");
    }
    for(size_t idx = 0; idx < getColNumber(); idx++){
        setValue(j, idx, getValue(j, idx) + k * getValue(i, idx));
    }


}

Matrix Matrix::operator*(const Matrix& rhs) const{
    if(rhs.getRowNumber() != getColNumber()){
        throw std::runtime_error("Matrix size mismatch!");
    }
    Matrix result(getRowNumber(), rhs.getColNumber());
    for(size_t i = 0; i < result.getRowNumber(); i++){
        for(size_t j = 0; j < result.getColNumber(); j++){
            double value = 0;
            for(size_t k = 0; k < getColNumber(); k++){
                value += getValue(i, k) * rhs.getValue(k, j);
            }
            result.setValue(i, j, value);
            
        }
    }
    return result;
}

void Matrix::scaleRow(size_t i, double k){
    if (i >= getRowNumber()) {
        throw std::out_of_range("Row index out of range in scaleRow");
    }
    for(size_t j = 0; j < getColNumber(); j++){
        setValue(i, j, getValue(i, j) * k);
    }
}

void Matrix::swap(size_t i, size_t j){
    if (i >= getRowNumber() || j >= getRowNumber()) {
        throw std::out_of_range("Row index out of range in scaleRow");
    }
    for(size_t k = 0; k < getColNumber(); k++){
        double iValue = getValue(i, k);
        double jValue = getValue(j, k);
        setValue(i, k, jValue);
        setValue(j, k, iValue);
    }

}

Matrix Matrix::inverseMatrix(){
    // 高斯-约旦消元法求逆
    // 1. 对每一列，从当前行到最后一行找最大值作为主元，避免数值不稳定
    // 2. 交换当前行和主元行
    // 3. 把当前行缩放使主元为1
    // 4. 对所有其他行消去当前列
    // 完成后 right 即为原矩阵的逆
    if(getColNumber() != getRowNumber()){
        throw std::runtime_error("The number of rows and columns does not match. Can not inverse matrix.");
    }
    Matrix left = *this;
    Matrix right(left.getColNumber(), UNIT_MATRIX);

    for(size_t i = 0; i < left.getColNumber(); i++){
        //step 1： Find the max abs value's index in colum i.
        double maxAbsValue = 0;
        size_t maxAbsIndex = i;
        for(size_t j = i; j < left.getColNumber(); j++){
            if(std::abs(left.getValue(j, i)) > maxAbsValue){
                maxAbsValue = std::abs(left.getValue(j, i));
                maxAbsIndex = j;
            }
        }

        //step 2: if maxAbsValue is 0 or very small, the matrix can't be inversed
        if(std::abs(maxAbsValue) < 1e-12){
            throw std::runtime_error("Matrix is singular or nearly singular.");
        }

        //step 3: swap the maxAbs row and the i row for both left and right
        right.swap(i, maxAbsIndex);
        left.swap(i, maxAbsIndex);
        

        //step 4: Scale the main value to 1;
        right.scaleRow(i, 1.0 / left.getValue(i, i));
        left.scaleRow(i, 1.0 / left.getValue(i, i));
        

        //step 5: use the row to minus other row, make this colum is 0 except i, i
        for(size_t j = 0; j < left.getColNumber(); j++){
            if( j == i) continue;
            right.addRowMultiple(i, j, -left.getValue(j, i));
            left.addRowMultiple(i, j, -left.getValue(j, i));
        }






    }
    return right;
    


}

}