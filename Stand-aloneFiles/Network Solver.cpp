#include "bits/stdc++.h"
using namespace std;

typedef complex<long double> cmx;

const long double eps = 1e-9;
typedef vector<vector<cmx>> matrix;
class Matrix
{
	public:
		
		int N,M;
		matrix mat;
		bool inverse_exists;

		Matrix(int n,int m)
		{
			N=n;
			M=m;
			mat.clear();
			mat.resize(N,vector<cmx>(M,0));
			if(N!=M)
				inverse_exists=0;
		}

		void reset()
		{
			mat.clear();
			N=0;
		}

		void Print()
		{
			cout<<"Matrix:"<<endl;
			for(int i=0;i<N;i++)
			{
				for(int j=0;j<M;j++)
				{
					cout<<mat[i][j]<<" ";
				}
				cout<<endl;
			}
		}

		Matrix Gauss_Jordan_Elimination(Matrix X)
		{
			X.inverse_exists=0;
			for(int i=0;i<N;i++)
			{
				int select=i;
				for(int j=i;j<N;j++)
				{
					if(abs(X.mat[j][i])>abs(X.mat[select][i]))
					{
						select=j;
					}
				}
				if(abs(X.mat[select][i])<eps)
					return X;
				for(int j=i;j<2*N;j++)
					swap(X.mat[i][j],X.mat[select][j]);
				for(int j=0;j<N;j++)
				{
					if(j!=i)
					{
						cmx tmp=X.mat[j][i]/X.mat[i][i];
						for(int k=0;k<2*N;k++)
						{
							X.mat[j][k]-=tmp*X.mat[i][k];
						}
					}
				}
				cmx tmp=X.mat[i][i];
				for(int j=0;j<2*N;j++)
				{
					X.mat[i][j]/=tmp;
				}
			}
			X.inverse_exists=1;
			return X;
		}

		Matrix Inverse()
		{
			Matrix res(N,N);
			res.inverse_exists=0;
			if(N!=M)
				return res;
			Matrix B(N,N+N);
			for(int i=0;i<N;i++)
			{
				for(int j=0;j<N;j++)
				{
					B.mat[i][j]=mat[i][j];
				}
				B.mat[i][N+i]=1;
			}
			B=Gauss_Jordan_Elimination(B);
			if(B.inverse_exists==0)
				return res;
			res.inverse_exists=1;
			for(int i=0;i<N;i++)
				for(int j=0;j<N;j++)
					res.mat[i][j]=B.mat[i][j+N];
			return res;
		}
};

void print(matrix A){
	int n = A.size();
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

struct edges{
	int to, v, r, id, dir;
	edges() {}
	edges(int to, int v, int r, int id = 0, int dir = 0) : to(to), v(v), r(r), id(id), dir(dir) {}
};
vector<vector<edges>> g;
vector<bool> vis, taken;
vector<pair<int,edges>> E;
void dfs(int u){
	vis[u] = 1;
	for(auto x : g[u]){
		if(!vis[x.to]){
			taken[x.id] = 1;
			dfs(x.to);
		}
	}
}
bool dfs1(int u,int v, vector<cmx> &temp, cmx &tot){
	vis[u] = 1;
	if(u == v) return 1;
	for(auto x : g[u]){
		if(!vis[x.to] && taken[x.id]){
			if(dfs1(x.to, v, temp, tot)){
				tot += x.v;
				temp[x.id - 1] += x.dir * x.r;
				return 1;
			}
		}
	}
	return 0;
}
matrix T;
void solve(){
	cout << "Enter number of nodes: " << endl;
	int n;
	cin >> n;
	g.resize(n + 1);
	vis.resize(n + 1, 0);
	cout << "Enter number of branches: " << endl;
	int m;
	cin >> m;
	taken.resize(m + 1, 0);
	E.resize(m + 1);
	cout << "Enter the branches, voltages and resistance" << endl;
	for(int i = 1; i <= m; i++){
		int u, v, V, r;
		cin >> u >> v >> V >> r;
		g[u].push_back({v, V, r, i, 1});
		g[v].push_back({u, -V, r, i, -1});
		E[i] = {u, {v, V, r}};
	}
	dfs(1);
	vector<vector<cmx>> A(m);
	int cc = 0;
	vector<cmx> B(m, 0);
	for(int i = 1; i <= m; i++){
		if(!taken[i]){
			fill(vis.begin(), vis.end(), 0);
			vector<cmx> temp(m, 0);
			cmx tot = 0;
			dfs1(E[i].first, E[i].second.to, temp, tot);
			tot -= E[i].second.v;
			temp[i - 1] = -E[i].second.r;
			B[cc] = tot;
			A[cc++] = temp;
		}
	}
	T = matrix(m, vector<cmx>(m, 0));
	for(int i = 1; i <= n; i++){
		vector<cmx> temp(m);
		for(auto x : g[i])
			temp[x.id - 1] = x.dir;
		for(int j = 0; j < m; j++){
			if(abs(temp[j]) > eps){
				if(abs(T[j][j]) < eps){
					T[j] = temp;
					assert(cc < m);
					A[cc++] = temp;
					break;
				}
				else{
					cmx z = temp[j] / T[j][j];
					for(int k = j; k < m; k++){
						temp[k] -= z * T[j][k];
					}
				}
			}
		}
	}
	Matrix M(m, m);
	M.mat = A;
	matrix invA = M.Inverse().mat;
	vector<cmx> I(m);
	for(int i = 0; i < m; i++){
		 for(int j = 0; j < m; j++){
		 	I[i] += invA[i][j] * B[j];
		 }
	}
	for(int i = 0; i < m; i++){
		cout << "I(" << i + 1 << ") = " << I[i] << endl;
	}
	cout << endl;
}
signed main(){
    #ifdef LOCAL
        freopen("inp.txt","r", stdin);
        freopen("out.txt", "w", stdout);
    #endif
    cout << setprecision(10) << fixed;
    solve();
    return 0;
}
