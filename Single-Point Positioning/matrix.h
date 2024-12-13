#pragma once
#include<iostream>
#include<math.h>
#include<vector>
#include <stdexcept>
using namespace std;

class Matrix {
public:
	int rows;//h
	int cols;//l
	vector<vector<double>>data;

	Matrix(int r = 1, int c = 1, const double* arr =NULL) : rows(r), cols(c), data(r, vector<double>(c, 0)) {
		if (arr != NULL) {
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					data[i][j] = arr[i * cols + j];
				}
			}
		}
	}

	int getRows() const {
		return rows;
	}
	int getCols() const {
		return cols;
	}	
	// 动态调整矩阵大小
	void resizeMatrix(int r, int c) {
		rows = r;
		cols = c;
		data.resize(r); // 调整行数
		for (int i = 0; i < r; ++i) {
			data[i].resize(c, 0); // 调整每行的列数，并初始化为 0
		}
	}
	//设置矩阵
	void setMatrix(int r, int c, double* arr) {
		rows = r; // 更新矩阵行数
		cols = c; // 更新矩阵列数
		resizeMatrix(r, c);
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				data[i][j] = arr[i * c + j];
			}
		}
	}
	// 清空矩阵，将所有元素设置为0
	void clear() {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				data[i][j] = 0;
			}
		}
	}
	// 修改矩阵中的某个元素
	void setElement(int r, int c, double value) {
		if (r >= 0 && r < rows && c >= 0 && c < cols) {
			data[r][c] = value;
		}
	}
	//输入矩阵
	void Matrix_Input() {
		for (int i = 0; i < this->getRows(); i++) {
			for (int j = 0; j < this->getCols(); j++) {
				double x;
				cin >> x;
				this->setElement(i, j, x);
			}
		}
	}
	//打印矩阵
	void display() const {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				cout << data[i][j] << " ";
			}
			cout << endl;
		}
	}
	//设置单位矩阵
	void setIdentity() {
		for (int i = 0; i < this->cols; i++) {
			this->setElement(i, i, 1.0);
		}
	}
	//设置对角矩阵
	void setDiagonal(double x) {
		for (int i = 0; i < this->cols; i++) {
			this->setElement(i, i, x);
		}
	}

	// 获取矩阵中的某个元素
	double getElement(int r, int c) const {
		if (r >= 0 && r < rows && c >= 0 && c < cols) {
			return data[r][c];
		}
		return 0; 
	}
	//矩阵加法
	Matrix operator+(const Matrix& other) const {
		if (rows != other.rows || cols != other.cols) {
			throw invalid_argument("error!");
		}

		Matrix result(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] + other.data[i][j];
			}
		}
		return result;
	}
	//矩阵减法
	Matrix operator-(const Matrix& other) const {
		if (rows != other.rows || cols != other.cols) {
			throw invalid_argument("error!");
		}

		Matrix result(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] - other.data[i][j];
			}
		}
		return result;
	}
	//实数与矩阵相乘
	Matrix operator*(double scalar) const {
		Matrix result(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				result.data[i][j] = data[i][j] * scalar;
			}
		}
		return result;
	}
	//矩阵与矩阵相乘(重载)
	Matrix operator*(const Matrix& other) const {
		if (cols != other.rows) {
			throw invalid_argument("error!");
		}

		Matrix result(rows, other.cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < other.cols; ++j) {
				result.data[i][j] = 0;
				for (int k = 0; k < cols; ++k) {
					result.data[i][j] += data[i][k] * other.data[k][j];
				}
			}
		}
		return result;
	}
	//矩阵转置
	Matrix transpose() const {
		Matrix result(cols, rows);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				result.setElement(j, i, data[i][j]);
			}
		}
		return result;
	}
	// 交换矩阵中的两行
	void swapRows(int row1, int row2) {
		if (row1 >= 0 && row1 < rows && row2 >= 0 && row2 <= rows) {
			for (int j = 0; j < cols; ++j) {
				double temp = data[row1-1][j];
				data[row1-1][j] = data[row2-1][j];
				data[row2-1][j] = temp;
			}
		}
		else {
			cout << "行索引超出范围，无法交换！" << endl;
		}
	}
	//矩阵求逆
	Matrix inverse() const {
		if (rows != cols) {
			cout << "行列数不等无法求逆" << endl;
		}
		//创建单位矩阵
		vector<double>m(rows*cols,0.0);
		Matrix identity;
		identity.setMatrix(rows, cols, m.data());
		identity.setIdentity();
		//复制原始矩阵
		Matrix copy;
		vector<double> flattened; 
		for (int i = 0; i < this->rows; ++i) {
			for (int j = 0; j < this->cols; ++j) {
				flattened.push_back(this->data[i][j]);
			}
		}
		copy.setMatrix(this->rows, this->cols, flattened.data());
		//高斯-约旦消元法
		for (int i = 0; i < rows; i++) {
			if (copy.getElement(i, i) == 0.0) {
				int j = i + 1;
				while (j < rows && copy.getElement(j, i) == 0.0) { j++; }
				if (j == rows) {cout<<"矩阵不可逆!";}
				copy.swapRows(i, j); identity.swapRows(i,j);
			}
			double pivot = copy.getElement(i, i);
			for (int j = 0; j < rows; j++) {
				copy.setElement(i, j, copy.getElement(i, j) / pivot);
				identity.setElement(i, j, identity.getElement(i, j) / pivot);
			}
			for (int j = 0; j < rows; j++) {
				if (j != i) {
					double factor = copy.getElement(j, i);
					for (int k = 0; k < rows; k++){
						copy.setElement(j, k, copy.getElement(j, k)-factor* copy.getElement(i, j));
						identity.setElement(j, k, identity.getElement(j, k) - factor * identity.getElement(i, j));
					}
				}
			}
		}
		return identity;
	}







	
	//正则化
	Matrix regularizeMatrix( double lambda = 1e-8) {
		

		Matrix regularized;
		vector<double> flattened;
		for (int i = 0; i < this->rows; ++i) {
			for (int j = 0; j < this->cols; ++j) {
				flattened.push_back(this->data[i][j]);
			}
		}
		regularized.setMatrix(this->rows, this->cols, flattened.data());
		for (int i = 0; i < regularized.getRows(); ++i) {
			regularized.setElement(i, i, regularized.getElement(i, i) + lambda);
		}
		return regularized;
	}






	// 计算矩阵的迹
	double trace() const {
		if (rows != cols) {
			cout << "矩阵不是方阵，无法计算迹。" << endl;
			return 0;
		}

		double sum = 0;
		for (int i = 0; i < rows; ++i) {
			sum += data[i][i];
		}
		return fabs(sum);
	}
};


//加权最小二乘
Matrix reweightedLeastSquares(Matrix B,Matrix L,Matrix P) {
	Matrix X_estimated;
	X_estimated = (B.transpose() * P * B).inverse() * (B.transpose() * P * L);//X最优估计量
	return X_estimated;
}
Matrix Identity(int row, int col) {//生成单位矩阵
	Matrix A(row, col);
	for (int i = 0; i < A.getCols(); i++) {
		A.setElement(i, i, 1.0);
	}
	return A;
}

