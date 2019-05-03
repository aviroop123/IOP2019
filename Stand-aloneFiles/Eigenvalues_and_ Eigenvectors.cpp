#include <bits/stdc++.h>
using namespace std;

typedef complex<long double> cmx;

const long double eps = 1e-7;

typedef vector<vector<cmx>> matrix;

void print(matrix A){
    int n = A.size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(abs(A[i][j]) < eps)
                A[i][j] = 0;
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

matrix transpose(matrix A){
    int n = A.size();
    for(int i = 0; i < n; i++)
        for(int j = i+1; j < n; j++)
            swap(A[i][j], A[j][i]);
    return A;
}

matrix multiply(matrix A, matrix B){
    int n = A.size();
    matrix C(n, vector<cmx>(n, 0));
    for(int k = 0; k < n; k++)
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

matrix Gram_Schmidt(matrix A){
    int n = A.size();
    for(int j = 0; j < n; j++){     //Gram-Schmidt's orthogonalisation
        vector<cmx> v;
        cmx mag = 0;    
        for(int i = 0; i < n; i++){
            v.push_back(A[i][j]);   // v initialised as column of A
            mag += v[i] * v[i];
        }
        for(int i = 0; i < j; i++){
            cmx dot = 0;
            for(int k = 0; k < n; k++)
                dot += A[k][j] * A[k][i];   //dot product
            for(int k = 0; k < n; k++)
                v[k] -= dot * A[k][i];  // subtract the projection
        }
        mag = 0;
        for(int i = 0; i < n; i++)
            mag += v[i] * v[i];
        mag = sqrt(mag);
        for(int i = 0; i < n; i++)
            A[i][j] = cmx(-1, 0) * v[i] / mag;      //normalize v
    }
    return A;
}

vector<cmx> get_eigenvalues(matrix A){
    int n = A.size();
    int k = 80; // number of iterations
    matrix Q, R, eigenvectors(n, vector<cmx>(n, 0));
    vector<cmx> eigenvalues;
    while(k--){
        Q = Gram_Schmidt(A);
        R = multiply(transpose(Q), A);
        eigenvectors = multiply(eigenvectors, Q);
        A = multiply(R, Q);
    }
    for(int i = 0; i < n; i++)
        eigenvalues.push_back(A[i][i]);
    return eigenvalues;
}

matrix forward_elimination(matrix A){
    int n = A.size();
    for(int i = 0; i < n; i++){
        int piv = i; long double val = abs(A[i][i]);
        for(int j = i + 1; j < n; j++){ // pivoting
            if(abs(A[j][i]) > val)
                val = abs(A[j][i]), piv = j;
        }
        if(val < eps) continue;
        swap(A[i], A[piv]);
        for(int j = i + 1; j < n; j++){
            cmx val = A[j][i] / A[i][i];
            for(int k = i; k < n; k++)
                A[j][k] -= val * A[i][k];
        }
    }
    return A;
}

matrix get_eigenvectors(matrix A, vector<cmx> lambda){
    int n = A.size();
    matrix B(n, vector<cmx>(n));
    int c = 1, t = 0;
    for(auto x : lambda){
        matrix C = A;
        for(int i = 0; i < n; i++)
            C[i][i] -= x;
        C = forward_elimination(C);
        vector<cmx> ans(n, 0);
        for(int i = n - 1; i >= 0; i--){
            if(abs(C[i][i]) < eps)
                ans[i] = c++;
            else{
                for(int j = i + 1; j < n; j++)
                    ans[i] -= ans[j] * C[i][j];
                ans[i] /= C[i][i];
            }
        }
        for(int i = 0; i < n; i++)
            B[i][t] = ans[i];
        t++;
    }
    return B;
}
signed main(){
    matrix A = {{0.7, 0.2, 0.1}, {0.2, 0.5, 0.3}, {0, 0, 1}};
    vector<cmx> lambda = get_eigenvalues(A);
    for(auto x : lambda){
        cout << x << " ";
    }
    cout << endl;
    print(get_eigenvectors(A, lambda));
    return 0;
}
