// matrixmultiplication.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "omp.h"
#include <iostream>
#include <random>
#include <stdio.h>
#include <vector>;
#include <chrono>;

using namespace std;

long long RunSequencialMultiplication(int n1, int n2, int n3);
long long RunBlockMultiplication(int n1, int n2, int n3, int block_side);
vector<int> MultipleMatrixParallel(vector<int> m1, vector<int> m2, int n1, int n2, int n3, int chunk);
long long RunParallelSequentialMultiplication(int n1, int n2, int n3, int chunk);
long long RunBlockMultiplicationParallel(int n1, int n2, int n3, int block_side, int chunk);
vector<vector<int>> MultipleMatrixBlocksParallel(vector<vector<int>> m1_blocks, vector<vector<int>> m2_blocks, int n1, int n2, int n3, int block_size, int chunk);
int GenerateRandomValue();
void GenerateMatrix(int row_number, int column_number, vector<int> &m);
void PrintMatrix(int row_number, int column_number, vector<int> m);
vector<int> MultipleMatrix(vector<int> m1, vector<int> m2, int n1, int n2, int n3);
vector<vector<int>> MultipleMatrixBlocks(vector<vector<int>> m1_blocks, vector<vector<int>> m2_blocks, int n1, int n2, int n3, int block_size);
vector<int> AddMatrixes(vector<int> m1, vector<int> m2, int row_number, int column_number);
void GenerateMatrixBlocks(vector<vector<int>> &m_blocks, int row_number, int column_number, int block_side);

int main() {
	
	const int n = 100, m = 100;
	const int block_size = 10;

	long long sequentialTime = RunSequencialMultiplication(n, n, n);
	long long blockTime = RunBlockMultiplication(n, n, n, block_size);
	long long sequentialParallelTime = RunParallelSequentialMultiplication(n, n, n, 10);
	long long blockParallelTime = RunBlockMultiplicationParallel(n, n, n, block_size, 10);

	cout << "sequential time " << sequentialTime << endl;
	cout << "block time " << blockTime << endl;
	cout << "sequential parallel time " << sequentialParallelTime << endl;
	cout << "block parallel time " << blockParallelTime << endl;
	
	system("pause");
    return 0;
}

long long RunSequencialMultiplication(int n1, int n2, int n3) {
	vector<int> m1 = {}, m2 = {}, m3 = {};
	m3.resize(n1, n3);

	GenerateMatrix(n1, n2, m1);
	GenerateMatrix(n2, n3, m2);

	auto start = chrono::steady_clock::now();
	m3 = MultipleMatrix(m1, m2, n1, n2, n3);
	auto duration = chrono::duration_cast<chrono::milliseconds>
                            (std::chrono::steady_clock::now() - start);
	return duration.count();
}

long long RunBlockMultiplication(int n1, int n2, int n3, int block_side) {
	vector<vector<int>> m1_blocks = {};
	vector<vector<int>> m2_blocks = {};
	vector<vector<int>> m3_blocks = {};

	GenerateMatrixBlocks(m1_blocks, n1, n2, block_side);
	GenerateMatrixBlocks(m2_blocks, n2, n3, block_side);

	auto start = chrono::steady_clock::now();
	m3_blocks = MultipleMatrixBlocks(m1_blocks, m2_blocks, n1/block_side, n2/block_side, n3/block_side, block_side);
	auto duration = chrono::duration_cast<chrono::milliseconds>
		(std::chrono::steady_clock::now() - start);
	
	return duration.count();
}

long long RunBlockMultiplicationParallel(int n1, int n2, int n3, int block_side, int chunk) {
	vector<vector<int>> m1_blocks = {};
	vector<vector<int>> m2_blocks = {};
	vector<vector<int>> m3_blocks = {};

	GenerateMatrixBlocks(m1_blocks, n1, n2, block_side);
	GenerateMatrixBlocks(m2_blocks, n2, n3, block_side);

	auto start = chrono::steady_clock::now();
	m3_blocks = MultipleMatrixBlocksParallel(m1_blocks, m2_blocks, n1 / block_side, n2 / block_side, n3 / block_side, block_side, chunk);
	auto duration = chrono::duration_cast<chrono::milliseconds>
		(std::chrono::steady_clock::now() - start);

	return duration.count();
}

void GenerateMatrixBlocks(vector<vector<int>> &m_blocks, int row_number, int column_number, int block_side) {
	m_blocks.resize(row_number * column_number * block_side);

	int m_blocks_number = row_number * column_number / (block_side * block_side);

	for (int i = 0; i < m_blocks_number; i++) {
		vector<int> block = {};
		GenerateMatrix(block_side, block_side, block);
		m_blocks[i] = block;
	}
}

long long RunParallelSequentialMultiplication(int n1, int n2, int n3, int chunk) {
	vector<int> m1 = {}, m2 = {}, m3 = {};
	m3.resize(n1, n3);

	GenerateMatrix(n1, n2, m1);
	GenerateMatrix(n2, n3, m2);

	auto start = chrono::steady_clock::now();

	m3 = MultipleMatrixParallel(m1, m2, n1, n2, n3, chunk);

	auto duration = chrono::duration_cast<chrono::milliseconds>
		(std::chrono::steady_clock::now() - start);
	return duration.count();
}

int GenerateRandomValue() {
	int sign = (rand() % 2) == 0 ? -1 : 1;
	return (rand() % 100) * sign;
}

void GenerateMatrix(int row_number, int column_number, vector<int> &m) {
	m.resize(row_number * column_number);
	for (int i = 0; i < row_number; i++) {
		for (int j = 0; j < column_number; j++) {
			m[i* row_number + j] = GenerateRandomValue();
		}
	}
}

void PrintMatrix(int row_number, int column_number, vector<int> m) {
	for (int i = 0; i < row_number; i++) {
		for (int j = 0; j < column_number; j++) {
			std::cout << m[i * row_number + j] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

vector<int> MultipleMatrix(vector<int> m1, vector<int> m2, int n1, int n2, int n3) {
	vector<int> result_matrix = {};
	result_matrix.resize(n1 * n3);
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n3; j++) {
			int row_sum = 0;
			for (int k = 0; k < n2; k++) {
				row_sum += m1[i * n1 + k] * m2[k * n2 + j];
			}
			result_matrix[i * n1 + j] = row_sum;
		}
	}
	return result_matrix;
}

vector<int> MultipleMatrixParallel(vector<int> m1, vector<int> m2, int n1, int n2, int n3, int chunk) {
	int i, j, k, row_sum;
	vector<int> result_matrix = {};
	result_matrix.resize(n1 * n3);

	#pragma omp parallel shared(result_matrix, m1, m2, n1, n2, n3, chunk) private(i, j, k, row_sum) 
	#pragma omp for schedule (static, chunk)
		for (i = 0; i < n1; i++) {
			for (j = 0; j < n3; j++) {
				row_sum = 0;
				for (k = 0; k < n2; k++) {
					row_sum += m1[i * n1 + k] * m2[k * n2 + j];
				}
				result_matrix[i * n1 + j] = row_sum;
			}
		}
	return result_matrix;
}

vector<vector<int>> MultipleMatrixBlocks(vector<vector<int>> m1_blocks, vector<vector<int>> m2_blocks, int n1, int n2, int n3, int block_size) {
	vector<vector<int>> result_matrix = {};
	result_matrix.resize(n1 * n3 * block_size);
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n3; j++) {
			vector<int> result_matrix_block = {};
			result_matrix_block.resize(block_size * block_size);
			for (int k = 0; k < n2; k++) {
				vector<int> matrix_production = {};
				matrix_production.resize(block_size * block_size);
				matrix_production = MultipleMatrix(m1_blocks[i * n1 + k], m2_blocks[k * n2 + j], block_size, block_size, block_size);
				result_matrix_block = AddMatrixes(result_matrix_block, matrix_production, 2, 2);
			}
			result_matrix[i * n1 + j] = result_matrix_block;
		}
	}
	return result_matrix;
}

vector<vector<int>> MultipleMatrixBlocksParallel(vector<vector<int>> m1_blocks, vector<vector<int>> m2_blocks, int n1, int n2, int n3, int block_size, int chunk) {
	int i, j, k;
	vector<int> result_matrix_block, matrix_production;
	vector<vector<int>> result_matrix = {};
	result_matrix.resize(n1 * n3 * block_size);

	#pragma omp parallel shared(result_matrix, m1_blocks, m2_blocks, block_size, n1, n2, n3, chunk) private(i, j, k, result_matrix_block, matrix_production) 
	#pragma omp for schedule (static, chunk)
	for (i = 0; i < n1; i++) {
		for (j = 0; j < n3; j++) {
			result_matrix_block = {};
			result_matrix_block.resize(block_size * block_size);
			for (k = 0; k < n2; k++) {
				matrix_production = {};
				matrix_production.resize(block_size * block_size);
				matrix_production = MultipleMatrix(m1_blocks[i * n1 + k], m2_blocks[k * n2 + j], block_size, block_size, block_size);
				result_matrix_block = AddMatrixes(result_matrix_block, matrix_production, 2, 2);
			}
			result_matrix[i * n1 + j] = result_matrix_block;
		}
	}
	return result_matrix;
}


vector<int> AddMatrixes(vector<int> m1, vector<int> m2, int row_number, int column_number) {
	vector<int> result_matrix = {};
	result_matrix.resize(row_number * column_number);
	for (int i = 0; i < row_number; i++) {
		for (int j = 0; j < column_number; j++) {
			result_matrix[i * row_number + j] = m1[i * row_number + j] + m2[i * row_number + j];
		}
	}
	return result_matrix;
}




