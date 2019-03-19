#include <bits/stdc++.h>

using namespace std;

#define ll long long
#define N 10004

ll a[N][N], b[N][N], c[N][N] = {0};

int main(){
	int n, m1, m2, m, p;
	cout<<"Enter dimensions\n";
	cin>>n>>m1>>m2>>p;
	if(m1 != m2){
		cout<<"Matrix multiplication not possible!\n";
		return 0;
	}
	m = m1;
	cout<<"Enter Matrices in order, row first\n";
	int i, j, k;
	for(i=1; i<=n; i++){
		for(j=1; j<=m; j++){
			cin>>a[i][j];
		}
	}
	for(i=1; i<=m; i++){
		for(j=1; j<=p; j++){
			cin>>b[i][j];
		}
	}
	for(i=1; i<=n; i++){
		for(j=1; j<=m; j++){
			for(k=1; k<=p; k++){
				c[i][k] = c[i][k] + a[i][j]*b[j][k];
			}
		}
	}
	cout<<"Multiplication result\n";
	for(i=1; i<=n; i++){
		for(j=1; j<=p; j++){
			cout<<c[i][j]<<" ";
		}
		cout<<"\n";
	}
	return 0;
}
