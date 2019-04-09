/*
#pragma once

#include <array>

template <class T, size_t rows, size_t columns>
class Matrix {
public:
	Matrix();
	~Matrix();

	T& Get(int i, int j);
	T Get(int i, int j) const;

	void Set()

	void Fill();
private:
	std::array<T, rows * columns> data_;
};


template <typename T, size_t rows, size_t columns>
Matrix<T, rows, columns>::Matrix() {
}

template <typename T, size_t rows, size_t columns>
Matrix<T, rows, columns>::~Matrix() {
}

template <typename T, size_t rows, size_t columns>
T& Matrix<T, rows, columns>::Get(int i, int j) {
	return data_[i * columns + j];
}

template <typename T, size_t rows, size_t columns>
T Matrix<T, rows, columns>::Get(int i, int j) const {
	return data_[i * columns + j];
}

template <typename T, size_t rows, size_t columns>
void Matrix<T, rows, columns>::Fill() {

}
*/