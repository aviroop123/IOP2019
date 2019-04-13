#include <bits/stdc++.h>
using namespace std;
 
typedef long double f80;

#define matrix vector<vector<f80>>

void print(matrix A){
	int n = A.size();
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			if(fabs(A[i][j]) < 1e-5)
				A[i][j] = 0;
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

matrix transpose(matrix A){
	int n = A.size();
	for(int i = 0; i < n; i++){
		for(int j = i+1; j < n; j++)
			swap(A[i][j], A[j][i]);
	}
	return A;
}

matrix multiply(matrix A, matrix B){
	int n = A.size();
	matrix C(n, vector<f80>(n, 0));
	for(int k = 0; k < n; k++){
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				C[i][j] += A[i][k]*B[k][j];
	}
	return C;
}

matrix Gram_Schmidt(matrix A){
	int n = A.size();
	for(int j = 0; j < n; j++){		//Gram-Schmidt's orthogonalisation
		vector<f80> v;
		f80 mag = 0;	
		for(int i = 0; i < n; i++){
			v.push_back(A[i][j]);	// v initialised as column of A
			mag += v[i]*v[i];
		}
		for(int i = 0; i < j; i++){
			f80 dot = 0;
			for(int k = 0; k < n; k++)
				dot += A[k][j]*A[k][i]; 	//dot product
			for(int k = 0; k < n; k++)
				v[k] -= dot*A[k][i];	// subtract the projection
		}
		mag = 0;
		for(int i = 0; i < n; i++)
			mag += v[i]*v[i];
		mag = sqrt(mag);
		for(int i = 0; i < n; i++)
			A[i][j] = -1*v[i]/mag;		//normalize v
	}
	return A;
}

signed main() {
//	freopen("input.txt","r",stdin);
	int n;
	matrix A;
	cout << "Enter size of matrix:\n";
	cin >> n;
	A.resize(n);
	cout << "Enter matrix:\n";
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			f80 p;	cin >> p;
			A[i].push_back(p);
		}
	}
	int k = 30;
	matrix Q, R, eigenvectors(n, vector<f80>(n, 0));
	for(int i = 0; i < n; i++)
		eigenvectors[i][i] = 1;
	while(k--){
		Q = Gram_Schmidt(A);
		R = multiply(transpose(Q), A);
		eigenvectors = multiply(eigenvectors, Q);
		A = multiply(R, Q);
	}
	print(A);
	print(eigenvectors);
	return 0;
}
