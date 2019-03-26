#include <bits/stdc++.h>

using namespace std;

int seeds[] = {41661, 82263, 31243, 70828, 46412, 25622, 67921, 74947, 78627, 84862, 66061, 70145, 95453, 43355, 37710, 99532, 59418, 11642, 36212, 82328};

vector < vector <int> > add(vector < vector <int> > x, vector < vector <int> > y, int s, int f){
	vector < vector <int> > ans(f-s+1);
	for(int i=s; i<=f; i++){
		for(int j=s; j<=f; j++){
			ans[i-s].push_back(x[i][j] + y[i][j]);
		}
	}
	return ans;
}

void print(vector < vector <int> > &a){
	for(auto it : a){
		for(auto val : it){
			cout<<val<<" ";
		}
		cout<<"\n";
	}
	cout<<"\n";
	return;
}

vector < vector <int> > sub(vector < vector <int> > x, vector < vector <int> > y, int s, int f){
	vector < vector <int> > ans(f-s+1);
	for(int i=s; i<=f; i++){
		for(int j=s; j<=f; j++){
			ans[i-s].push_back(x[i][j] - y[i][j]);
		}
	}
	return ans;
}

vector < vector <int> > strassens_multiply(vector < vector <int> > x, vector < vector <int> > y, int s, int f){
	if(s + 1 == f){
		vector < vector <int> > c(2);
		c[0].push_back(x[s][s]*y[s][s] + x[s][f]*y[f][s]);
		c[0].push_back(x[s][s]*y[s][f] + x[s][f]*y[f][f]);
		c[1].push_back(x[f][s]*y[s][s] + x[f][f]*y[f][s]);
		c[1].push_back(x[f][s]*y[s][f] + x[f][f]*y[f][f]);
		return c;
	}
	int sz = (s + f)>>1, i, j;
	vector < vector <int> > a11(sz+1), a12(sz+1), a21(sz+1), a22(sz+1), b11(sz+1), b12(sz+1), b21(sz+1), b22(sz+1);
	for(i=s; i<=sz; i++){
		for(j=s; j<=sz; j++){
			a11[i-s].push_back(x[i][j]);
			b11[i-s].push_back(y[i][j]);
		}
	}
	for(i=s; i<=sz; i++){
		for(j=sz+1; j<=f; j++){
			a12[i-s].push_back(x[i][j]);
			b12[i-s].push_back(y[i][j]);
		}
	}
	for(i=sz+1; i<=f; i++){
		for(j=s; j<=sz; j++){
			a21[i-sz-1].push_back(x[i][j]);
			b21[i-sz-1].push_back(y[i][j]);
		}
	}
	for(i=sz+1; i<=f; i++){
		for(j=sz+1; j<=f; j++){
			a22[i-sz-1].push_back(x[i][j]);
			b22[i-sz-1].push_back(y[i][j]);
		}
	}
	// print(a11); print(a12);
	// print(a21); print(a22);
	// print(b11); print(b12);
	// print(b21); print(b22);
	vector < vector <int> > m1(sz+1), m2(sz+1), m3(sz+1), m4(sz+1), m5(sz+1), m6(sz+1), m7(sz+1);
	m1 = strassens_multiply(add(a11, a22, 0, sz), add(b11, b22, 0, sz), 0, sz);
	m2 = strassens_multiply(add(a21, a22, 0, sz), b11, 0, sz);
	m3 = strassens_multiply(a11, sub(b12, b22, 0, sz), 0, sz);
	m4 = strassens_multiply(a22, sub(b21, b11, 0, sz), 0, sz);
	m5 = strassens_multiply(add(a11, a12, 0, sz), b22, 0, sz);
	m6 = strassens_multiply(sub(a21, a11, 0, sz), add(b11, b12, 0, sz), 0, sz);
	m7 = strassens_multiply(sub(a12, a22, 0, sz), add(b21, b22, 0, sz), 0, sz);
	vector < vector <int> > c11(sz+1), c12(sz+1), c21(sz+1), c22(sz+1);
	c11 = add(m1, add(m4, sub(m7, m5, 0, sz), 0, sz), 0, sz);
	c12 = add(m3, m5, 0, sz);
	c21 = add(m2, m4, 0, sz);
	c22 = add(m1, add(m3, sub(m6, m2, 0, sz), 0, sz), 0, sz);
	vector < vector <int> > c(f - s + 1);
	for(i=s; i<=sz; i++){
		for(j=s; j<=sz; j++){
			c[i-s].push_back(c11[i-s][j-s]);
		}
	}
	for(i=s; i<=sz; i++){
		for(j=sz+1; j<=f; j++){
			c[i-s].push_back(c12[i-s][j-sz-1]);
		}
	}
	for(i=sz+1; i<=f; i++){
		for(j=s; j<=sz; j++){
			c[i-sz+((f-s)>>1)].push_back(c21[i-sz-1][j-s]);
		}
	}
	for(i=sz+1; i<=f; i++){
		for(j=sz+1; j<=f; j++){
			c[i-sz+((f-s)>>1)].push_back(c22[i-sz-1][j-sz-1]);
		}
	}
	return c;
}

int main(){
	clock_t clk;	
	int n, i, j;
	for(int t=1; t<=10; t++){
		n = (1<<t);
		srand(seeds[t-1]);
		vector < vector <int> > a(n), b(n);
		for(i=0; i<n; i++){
			for(j=0; j<n; j++){
				a[i].push_back(rand()%1000);
				b[i].push_back(rand()%1000);
			}
		}
		clk = clock();
		vector < vector <int> > c = strassens_multiply(a, b, 0, n-1);
		clk = clock() - clk;
		cout<<fixed<<setprecision(6)<<n<<" "<<((double)clk)/CLOCKS_PER_SEC<<"\n";
	}
	return 0;
}
