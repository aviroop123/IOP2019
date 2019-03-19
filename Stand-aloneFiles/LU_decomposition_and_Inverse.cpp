// Performs LU decomposition - Dolittle, Crout, Cholesky and computes inverse of a matrix using these decompositions.
#include<bits/stdc++.h>
using namespace std;
typedef long double f80;
void print(vector<vector<f80>> A) {
    int n = A.size();
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
vector<vector<f80>> operator *(vector<vector<f80>> A, vector<vector<f80>> B){
    int n = A.size();
    vector<vector<f80>> temp(n);
    for(int i = 0; i < n; i++) {
        temp[i].resize(n, 0);
    }
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
                temp[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return temp;
}
class Matrix{
public:
    vector<vector<f80>> A, L, U, IA;
    bool decomposed;
    int n;
    Matrix(int n) : n(n) {
        decomposed = 0;
        A.resize(n);
        for(int i = 0; i < n; i++){
            A[i].resize(n, 0);
        }
    }
    void reset(){
        L.clear();
        U.clear();
        L.resize(n);
        U.resize(n);
        for(int i = 0; i < n; i++){
            L[i].resize(n, 0);
            U[i].resize(n, 0);
        }
    }
    void LU_Crout() { // U[i][j] = 1 for 0 <= i < n
        reset();
        decomposed = 1;
       for (int j = 0; j < n; j++) {
        U[j][j] = 1;
            for (int i = j; i < n; i++) {
                L[i][j] = A[i][j];
                for (int k = 0; k < j; k++) {
                    L[i][j] -= L[i][k] * U[k][j];   
                }
            }
            for (int i = j; i < n; i++) {
                U[j][i] = A[j][i];
                for(int k = 0; k < j; k++) {
                    U[j][i] -= L[j][k] * U[k][i];
                }
                U[j][i] /= L[j][j];
            }
        }
    }
    void LU_Dolittle() { // L[i][i] = 1 for 0 <= i < n
        reset();
        decomposed = 1;
        for(int i = 0; i < n; i++) {
            L[i][i] = 1;
            for(int j = i; j < n; j++) {
                U[i][j] = A[i][j];
                for(int k = 0; k < i; k++) {
                    U[i][j] -= L[i][k] * U[k][j];
                }
            }
            for(int j = i + 1; j < n; j++) {
                L[j][i] = A[j][i];
                for(int k = 0; k < i; k++) {
                    L[j][i] -= L[j][k] * U[k][i];
                }
                L[j][i] /= U[i][i];
            }
        }
    }
    void LU_Cholesky() { // U = L^T, possible only if A is symmetric
        reset();
        decomposed = 1;
        for(int i = 0; i < n; i++) {
            L[i][i] = A[i][i];
            for(int j = 0; j < i; j++) {
                L[i][i] -= L[i][j] * L[i][j];
            }
            assert(L[i][i] >= 0);
            L[i][i] = sqrtl(L[i][i]);
            for(int j = i + 1; j < n; j++) {
                L[j][i] = A[j][i];
                for(int k = 0; k < i; k++) {
                    L[j][i] -= L[j][k] * L[i][k];
                }
                L[j][i] /= L[i][i];
            }
        }
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                U[i][j] = L[j][i];
            }
        }
    }
    void operator = (const Matrix &M){
        A = M.A, L = M.L, U = M.U, IA = M.IA, n = M.n;
    }
    void sub(vector<vector<f80>> &a, int r2, int r1, f80 k){ // r2 = r2 - k * r1
        for(int i = 0; i < n; i++) {
            a[r2][i] -= k * a[r1][i];
        }
    }
    void rowmul(vector<vector<f80>> &a, int r, f80 k){ // r = r * k;
        for(int i = 0; i < n; i++) {
            a[r][i] *= k;
        }
    }
    void inverse() { // Call any LU - decomposition method first, otherwise DoLittle by default
        if(!decomposed) LU_Dolittle(); // by default
        IA.clear();
        IA.resize(n);
        for(int i = 0; i < n; i++) {
            IA[i].resize(n, 0);
            IA[i][i] = 1;
        }
        auto T = L;
        for(int i = 0; i < n; i++){ // forward substitution
            f80 val = 1.0 / T[i][i];
            rowmul(T, i, val);
            rowmul(IA, i, val);
            for(int j = i + 1; j < n; j++) {
                f80 val = T[j][i];
                sub(T, j, i, val);
                sub(IA, j, i, val);
            }
        }
        T = U;
        for(int i = n - 1; i >= 0; i--) { // back substitution
            f80 val = 1.0 / T[i][i];
            rowmul(T, i, val);
            rowmul(IA, i, val);
            for(int j = i - 1; j >= 0; j--) {
                f80 val = T[j][i];
                sub(T, j, i, val);
                sub(IA, j, i, val);
            }
        }
    }
};

int main() {
    Matrix M(3);
    M.A = {{25, 15, -5}, {15, 18, 0}, {-5, 0, 11}};
    M.LU_Cholesky();
    cout << "L - Matrix: " << endl;
    print(M.L);
    cout << "U - Matrix " << endl;
    print(M.U);
    M.inverse();
    cout << "Inverse of Matrix : " << endl;
    print(M.IA);
    return 0;
}
