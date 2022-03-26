#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <stdexcept>

namespace KLinear {
    /** Matrix*/
    struct matrix {
        int row;
        int column;
        std::vector<std::vector<std::pair<int, int>>> element;
    };

    /** Determinant*/
    struct determinant {
        int order;
        std::vector<std::vector<std::pair<int, int>>> element;
    };

    /** Transposed Matrix*/
    struct TransMatrix {
        matrix origin;
    };

    /** Transposed Determinant*/
    struct TransDeterminant {
        determinant origin;
    };

    /** Fraction Reduction */
    std::pair<int, int> reduction(std::pair<int, int> fraction) {
        if (fraction.second == 0)throw std::invalid_argument("Invalid fraction expression.");
        if (fraction.first == 0)return fraction;
        int gcd = std::__gcd(fraction.first, fraction.second);
        if (gcd == 1)return fraction;
        fraction.first /= gcd;
        fraction.second /= gcd;
        return fraction;
    }

    /** Fraction Add A Number*/
    std::pair<int, int> add(std::pair<int, int> fraction, int num) {
        if (fraction.second == 0)throw std::invalid_argument("Invalid fraction expression.");
        fraction.first += num * fraction.second;
        return reduction(fraction);
    }

    /** Fraction Add A Fraction*/
    std::pair<int, int> add(std::pair<int, int> fraction, std::pair<int, int> num) {
        if (fraction.second == 0 || num.second == 0)throw std::invalid_argument("Invalid fraction expression.");
        (fraction.first *= num.second) += num.first * fraction.second;
        fraction.second *= num.second;
        return reduction(fraction);
    }

    /** Fraction Minus A Number*/
    std::pair<int, int> minus(std::pair<int, int> fraction, int num) {
        if (fraction.second == 0)throw std::invalid_argument("Invalid fraction expression.");
        fraction.first -= num * fraction.second;
        return reduction(fraction);
    }

    /** Fraction Minus A Fraction*/
    std::pair<int, int> minus(std::pair<int, int> fraction, std::pair<int, int> num) {
        if (fraction.second == 0 || num.second == 0)throw std::invalid_argument("Invalid fraction expression.");
        (fraction.first *= num.second) -= num.first * fraction.second;
        fraction.second *= num.second;
        return reduction(fraction);
    }

    /** Fraction Multiply A Number*/
    std::pair<int, int> multiply(std::pair<int, int> fraction, int num) {
        if (fraction.second == 0)throw std::invalid_argument("Invalid fraction expression.");
        if (num != 0 || fraction.first == 0) {
            int gcd = abs(std::__gcd(fraction.second, num));
            fraction.second /= gcd;
            num /= gcd;
            fraction.first *= num;
        } else return std::make_pair(0, 1);
        return fraction;
    }

    /** Fraction Multiply A Fraction*/
    std::pair<int, int> multiply(std::pair<int, int> fraction, std::pair<int, int> num) {
        if (fraction.second == 0 || num.second == 0)throw std::invalid_argument("Invalid fraction expression.");
        if (num.first != 0 || fraction.first == 0) {
            std::pair<int, int> R_fraction = reduction(fraction);
            std::pair<int, int> R_num = reduction(num);
            R_fraction.first *= R_num.first;
            R_fraction.second *= R_num.second;
            fraction = reduction(R_fraction);
        } else return std::make_pair(0, 1);
        return fraction;
    }

    /** Fraction Divide A Number*/
    std::pair<int, int> divide(std::pair<int, int> fraction, int num) {
        return multiply(fraction, std::make_pair(1, num));
    }

    /** Fraction Divide A Fraction*/
    std::pair<int, int> divide(std::pair<int, int> fraction, std::pair<int, int> num) {
        return multiply(fraction, std::make_pair(num.second, num.first));
    }

    /** Turn Square Matrix Into Determinant*/
    determinant transform(matrix square) {
        if (square.row != square.column) throw std::invalid_argument("Invalid determinant expression.");
        determinant res;
        res.order = square.row;
        res.element = square.element;
        return res;
    }

    /** Turn Transposed Square Matrix Into Transposed Determinant*/
    TransDeterminant transform(TransMatrix square) {
        TransDeterminant res;
        res.origin = transform(square.origin);
        return res;
    }

    /** Calculate The Value Of Determinant*/
    std::pair<int, int> det(determinant square) {
        for (int i = 0; i < square.order; ++i) {
            for (int j = square.order; j > i; --j) {
                std::pair<int, int> mid = divide(square.element[j][0], square.element[i][0]);
                for (int k = square.order - 1; k > j; --k) {
                    square.element[j][k] = minus(square.element[j][k], multiply(square.element[i][k], mid));
                }
            }
        }
        std::pair<int, int> res(1, 1);
        for (int i = 0; i < square.order; ++i) {
            if (square.element[i][i].first == 0) {
                res = std::make_pair(0, 1);
                break;
            }
            res = multiply(res, square.element[i][i]);
        }
        return res;
    }

    /** Calculate The Value Of The Determinant Of Transposed Matrix*/
    std::pair<int, int> det(TransDeterminant square){
        return det(square.origin);
    }

    /** Get A Minor*/
    determinant minor(determinant det, int i, int j) {
        determinant res;
        res.order = det.order - 1;
        res.element = det.element;
        std::swap(res.element[i], res.element[res.order]);
        res.element.pop_back();
        std::swap(res.element[i], res.element[res.order - 1]);
        for (int k = 0; k < res.order; ++k) {
            std::swap(res.element[k][j], res.element[k][res.order]);
            res.element[k].pop_back();
            std::swap(res.element[k][j], res.element[k][res.order - 1]);
        }
        return res;
    }

    /** Create A Null Matrix*/
    matrix nullMatrix(int row, int column) {
        matrix null;
        null.row = row;
        null.column = column;
        for (int j = 0; j < column; ++j) {
            null.element[0].push_back(std::make_pair(0, 1));
        }
        for (int i = 1; i < row; ++i) {
            null.element.push_back(null.element[0]);
        }
        return null;
    }

    /** Create A Square Null Matrix*/
    matrix nullMatrix(int order) {
        return nullMatrix(order, order);
    }

    /** Create A Identity Matrix*/
    matrix identityMatrix(int order) {
        matrix identity;
        std::vector<std::pair<int, int>> tempRow;
        tempRow.emplace_back(1, 1);
        for (int i = 2; i < order; ++i) {
            tempRow.emplace_back(0, 1);
        }
        identity.element.push_back(tempRow);
        for (int i = 1; i < order; ++i) {
            std::swap(tempRow[i], tempRow[i - 1]);
            identity.element.push_back(tempRow);
        }
        return identity;
    }

    /** Matrix Transpose*/
    matrix transpose(matrix A) {
        matrix At;
        At.row = A.column;
        At.column = A.row;
        if (At.row == At.column){
            At.element = A.element;
            for (int i = 0; i < A.row; ++i) {
                for (int j = 0; j < i; ++j) {
                    std::swap(At.element[i][j],At.element[j][i]);
                }
            }
        }
        else {
            for(int i = 0; i<A.row;++i){
                for (int j = 0; j < A.column; ++j) {
                    At.element[j][i] = A.element[i][j];
                }
            }
        }
        return At;
    }

}

