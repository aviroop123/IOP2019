#include <bits/stdc++.h>
 
using namespace std;

const long double eps=1e-9; 
typedef long double ld;


class Matrix
{
	public:
		
		int N,M;
		vector<vector<ld>>mat;
		bool inverse_exists;

		Matrix(int n,int m)
		{
			N=n;
			M=m;
			mat.clear();
			mat.resize(N,vector<ld>(M,0));
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
					cout<<fixed<<setprecision(5)<<mat[i][j]<<" ";
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
					if(fabs(X.mat[j][i])>fabs(X.mat[select][i]))
					{
						select=j;
					}
				}
				if(fabs(X.mat[select][i])<eps)
					return X;
				for(int j=i;j<2*N;j++)
					swap(X.mat[i][j],X.mat[select][j]);
				for(int j=0;j<N;j++)
				{
					if(j!=i)
					{
						ld tmp=X.mat[j][i]/X.mat[i][i];
						for(int k=0;k<2*N;k++)
						{
							X.mat[j][k]-=tmp*X.mat[i][k];
						}
					}
				}
				ld tmp=X.mat[i][i];
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

int main()
{
	Matrix A(3,3);	
	A.mat={{5,7,9},{4,3,8},{7,5,6}};
	Matrix B=A.Inverse();
	if(B.inverse_exists==0)
	{
		cout<<"Inverse Does Not Exist!!"<<endl;
	}
	else
	{
		B.Print();
	}
	return 0;
}
