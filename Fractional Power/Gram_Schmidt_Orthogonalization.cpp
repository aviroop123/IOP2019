#include <bits/stdc++.h>
using namespace std;
 
typedef long double f80;

void print(vector<vector<f80> > A){
	int n = A.size();
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++)
			cout << A[i][j] << " ";
		cout << endl;
	}
}

vector<vector<f80> > Gram_Schmidt(vector<vector<f80> > &A){
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
			A[i][j] = v[i]/mag;		//normalize v
	}
	return A;
}

signed main() {
	freopen("input.txt","r",stdin);
	int n;
	vector<vector<f80> > A;
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
	Gram_Schmidt(A);
	cout << "Matrix after Gram_Schmidt's orthogonalization:\n";
	print(A);
	return 0;
}
