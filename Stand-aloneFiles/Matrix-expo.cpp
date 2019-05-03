#include <bits/stdc++.h>
using namespace std;
 
typedef complex<long double> cmx;

const int mod = 1e9 + 7;

int sz;
const int NN = 101;

class matrix{
public:
   	cmx mat[NN][NN];
    matrix(){
        for(int i = 0; i < NN; i++)
            for(int j = 0; j < NN; j++)
                mat[i][j] = 0;
    }
    inline matrix operator * (const matrix &a){
        matrix temp;
        for(int i = 0; i < sz; i++)
            for(int j = 0; j < sz; j++){
                for(int k = 0; k < sz; k++){
                temp.mat[i][j] += (mat[i][k] * a.mat[k][j]);
            }
        }
        return temp;
    }
    inline matrix operator + (const matrix &a){
        matrix temp;
        for(int i = 0; i < sz; i++)
            for(int j = 0; j < sz; j++)
                temp.mat[i][j] = mat[i][j] + a.mat[i][j];
        return temp;
    }
    inline matrix operator - (const matrix &a){
        matrix temp;
        for(int i = 0; i < sz; i++)
            for(int j = 0; j < sz; j++)
                temp.mat[i][j] = mat[i][j] - a.mat[i][j];
        return temp;
    }
    inline void operator = (const matrix &b){
        for(int i = 0; i < sz; i++)
            for(int j = 0; j < sz; j++)
                mat[i][j] = b.mat[i][j];
    }
    inline void print(){
        for(int i = 0; i < sz; i++){
            for(int j = 0; j < sz; j++){
                cout << setprecision(15) << mat[i][j] << " ";
            }
            cout << endl;
        }
    }
};

matrix pow(matrix a, int k){
    matrix ans;
    for(int i = 0; i < sz; i++)
        ans.mat[i][i] = 1;
    while(k){
        if(k & 1)
            ans = ans * a;
        a = a * a;
        k >>= 1;
    }
    return ans;
}

signed main() {
	matrix a;
	int p;
	cout << "Enter size of matrix:\n";
	cin >> sz;
	cout << "Enter elements of matrix:\n";
	for(int i = 0; i < sz; i++){
		for(int j = 0; j < sz; j++){
			long double x, y;
			cin >> x >> y;
			a.mat[i][j] = {x, y};
		}
	}
	cout << "Enter power to be raised:\n";
	cin >> p;
	a = pow(a, p);
	cout << "Resultant matrix:\n";
	a.print();
	return 0;
}
