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
		ld det;

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

		ld Determinant(Matrix X)
		{
			ld res=1;
			for(int i=0;i<X.N;i++)
			{
				int select=i;
				while(select<X.N&&fabs(X.mat[select][i])<eps)
				{
					select++;
				}
				if(select==X.N)
					return 0;
				if(select!=i)
				{
					res*=-1;	
					for(int j=0;j<X.N;j++)
					{
						swap(X.mat[i][j],X.mat[select][j]);
					}
				}
				res*=X.mat[i][i];
				for(int j=i+1;j<X.N;j++)
				{
					ld tmp=X.mat[j][i]/X.mat[i][i];
					for(int k=0;k<X.N;k++)
					{
						X.mat[j][k]-=X.mat[i][k]*tmp;
					}
				}
			}
			return res;
		}

		Matrix GetCofactor(Matrix X,int I,int J)
		{
			Matrix res(X.N-1,X.N-1);
			int curi=0;
			for(int i=0;i<X.N;i++)
			{
				if(i==I)
					i++;
				if(i>=X.N)
					break;
				int curj=0;
				for(int j=0;j<X.N;j++)
				{
					if(j==J)
						j++;
					if(j>=X.N)
						break;
					res.mat[curi][curj]=X.mat[i][j];
					curj++;
				}
				curi++;
			}
			return res;
		}

		Matrix Adjoint(Matrix X)
		{
			Matrix res(X.N,X.N);
			for(int i=0;i<X.N;i++)
			{
				for(int j=0;j<X.N;j++)
				{
					res.mat[i][j]=Determinant(GetCofactor(X,i,j));
					if((i+j)%2==1)
						res.mat[i][j]*=-1;
				}
			}
			return res;
		}

		Matrix Inverse()
		{
			Matrix res(N,N);
			res.inverse_exists=0;
			if(N!=M)
				return res;
			Matrix B(N,N);
			B.mat=mat;
			ld det=Determinant(B);
			if(fabs(det)<eps)
				return res;
			res=Adjoint(B);
			res.inverse_exists=1;
			for(int i=0;i<N;i++)
			{
				for(int j=0;j<N;j++)
				{
					res.mat[i][j]/=det;
				}
			}
			for(int i=0;i<N;i++)
			{
				for(int j=i+1;j<N;j++)
				{
					swap(res.mat[i][j],res.mat[j][i]);
				}
			}
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
