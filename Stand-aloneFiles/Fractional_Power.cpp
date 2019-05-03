#include <bits/stdc++.h>
using namespace std;

typedef complex<long double> cmx;

const long double eps = 1e-7;

typedef vector<vector<cmx>> matrix;

vector<cmx> eigenvalues;
matrix D;
int N;

void print(matrix A){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(abs(A[i][j]) < eps)
                A[i][j] = 0;
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

matrix transpose(matrix A){
    for(int i = 0; i < N; i++)
        for(int j = i+1; j < N; j++)
            swap(A[i][j], A[j][i]);
    return A;
}

matrix multiply(matrix A, matrix B){
    matrix C(N, vector<cmx>(N, 0));
    for(int k = 0; k < N; k++)
        for(int i = 0; i < N; i++)
            for(int j = 0; j < N; j++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

matrix Gram_Schmidt(matrix A){		//To find matrix Q
    int n = A.size();
    for(int j = 0; j < N; j++){     //Gram-Schmidt's orthogonalisation
        vector<cmx> v;
        cmx mag = 0;    
        for(int i = 0; i < N; i++){
            v.push_back(A[i][j]);   // v initialised as column of A
            mag += v[i] * v[i];
        }
        for(int i = 0; i < j; i++){
            cmx dot = 0;
            for(int k = 0; k < N; k++)
                dot += A[k][j] * A[k][i];   //dot product
            for(int k = 0; k < N; k++)
                v[k] -= dot * A[k][i];  // subtract the projection
        }
        mag = 0;
        for(int i = 0; i < N; i++)
            mag += v[i] * v[i];
        mag = sqrt(mag);
        for(int i = 0; i < N; i++)
            A[i][j] =  -v[i] / mag;      //normalize v
    }
    return A;
}

vector<cmx> get_eigenvalues(matrix A){
    int k = 80; // number of iterations
    matrix Q, R, eigenvectors(N, vector<cmx>(N, 0));
    while(k--){
        Q = Gram_Schmidt(A);
        R = multiply(transpose(Q), A);
        eigenvectors = multiply(eigenvectors, Q);
        A = multiply(R, Q);
    }
    D = matrix(N, vector<cmx> (N, 0));
    for(int i = 0; i < N; i++){
    	D[i][i] = A[i][i];
        eigenvalues.push_back(A[i][i]);
    }
    return eigenvalues;
}

matrix forward_elimination(matrix A){
    for(int i = 0; i < N; i++){
        int piv = i; long double val = abs(A[i][i]);
        for(int j = i + 1; j < N; j++){ // pivoting
            if(abs(A[j][i]) > val)
                val = abs(A[j][i]), piv = j;
        }
        if(val < eps) continue;
        swap(A[i], A[piv]);
        for(int j = i + 1; j < N; j++){
            cmx val = A[j][i] / A[i][i];
            for(int k = i; k < N; k++)
                A[j][k] -= val * A[i][k];
        }
    }
    return A;
}

matrix get_eigenvectors(matrix A, vector<cmx> lambda){
    matrix B(N, vector<cmx>(N));
    int c = 1, t = 0;
    for(auto x : lambda){
        matrix C = A;
        for(int i = 0; i < N; i++)
            C[i][i] -= x;
        C = forward_elimination(C);
        vector<cmx> ans(N, 0);
        for(int i = N - 1; i >= 0; i--){
            if(abs(C[i][i]) < eps)
                ans[i] = c++;
            else{
                for(int j = i + 1; j < N; j++)
                    ans[i] -= ans[j] * C[i][j];
                ans[i] /= C[i][i];
            }
        }
        for(int i = 0; i < N; i++)
            B[i][t] = ans[i];
        t++;
    }
    return B;
}

void sub(matrix &a, int r2, int r1, cmx k){ // r2 = r2 - k * r1
    for(int i = 0; i < N; i++) {
        a[r2][i] -= k * a[r1][i];
    }
}

void rowmul(vector<vector<cmx>> &a, int r, cmx k){ // r = r * k;
    for(int i = 0; i < N; i++) {
        a[r][i] *= k;
    }
}

matrix L, U;

void LU_Dolittle(matrix A) { // L[i][i] = 1 for 0 <= i < n
    for(int i = 0; i < N; i++) {
        L[i][i] = 1;
        for(int j = i; j < N; j++) {
            U[i][j] = A[i][j];
            for(int k = 0; k < i; k++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
        for(int j = i + 1; j < N; j++) {
            L[j][i] = A[j][i];
            for(int k = 0; k < i; k++) {
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
    }
}

matrix inverse(matrix A) { // Call any LU - decomposition method first, otherwise DoLittle by default
   	LU_Dolittle(A);
    matrix IA(N, vector<cmx>(N, 0));
    for(int i = 0; i < N; i++) 
        IA[i][i] = 1;
   	matrix T = L;
    for(int i = 0; i < N; i++){ // forward substitution
        cmx val = conj(T[i][i])/(abs(T[i][i])*abs(T[i][i]));
        rowmul(T, i, val);
        rowmul(IA, i, val);
        for(int j = i + 1; j < N; j++) {
            cmx val = T[j][i];
            sub(T, j, i, val);
            sub(IA, j, i, val);
        }
    }
    T = U;
    for(int i = N - 1; i >= 0; i--) { // back substitution
        cmx val = conj(T[i][i])/(abs(T[i][i])*abs(T[i][i]));
        rowmul(T, i, val);
        rowmul(IA, i, val);
        for(int j = i - 1; j >= 0; j--) {
            cmx val = T[j][i];
            sub(T, j, i, val);
            sub(IA, j, i, val);
        }
    }
    return IA;
}

signed main(){
    matrix A = {{0.7, 0.2, 0.1}, {0.2, 0.5, 0.3}, {0, 0, 1}};
    N = A.size();
    cmx power;
    cout << "Enter fractional power to be raised:";
    cin >> power;
    L = matrix(N, vector<cmx>(N, 0));
    U = matrix(N, vector<cmx>(N, 0));
    vector<cmx> lambda = get_eigenvalues(A);
    matrix P = get_eigenvectors(A, lambda);
    matrix Q = inverse(P);
    for(int i = 0; i < N; i++)
    	D[i][i] = pow(D[i][i], power);
    A = multiply(P, D);
    A = multiply(A, Q);
    print(A);
    return 0;
}
