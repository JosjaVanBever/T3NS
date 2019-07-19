
#include <iostream>
#include <math.h>
using namespace std;

// g++ -O3 -std=c++17 testing.cpp -o testing.exe -lblas

extern "C" {
void dgemm_(char & transa, char & transb, int & m, int & n, int & k,
		double & alpha, double * a, int & lda, double * b, int & ldb,
		double & beta, double * c, int & ldc);
}


void test1();
void mattest();
void objecttest();
void arraytest();
void functieptrtest();
void casttest();
void alloctest();
void doubleptrtest();

class Bar {
	int i;
};

class Foo {
	public:
		Foo(int i) : i(i) {};
		int get_i() { return i; }
		void set_i(int j) { i = j; }
	private:
		int i;
};

int main() {
	int test[5] = {4};
	cout << test[2];
	doubleptrtest();
	return 0;
}

void ptrarraymanip(int ** test) {
	// do nothing
}

void doubleptrtest() {
	int *test[4];
	ptrarraymanip(test);
}

void casttest() {
	cout << int(round(4.99)) % 2 << endl;
	// cout << 4.5 % 2 << endl;
}

void alloctest() {
	int * test = (int * ) malloc(3 * sizeof(int));
	test[0] = 5;
	cout << test[0] << endl;
	free(test);
}

int add(int i, int j) {
	return i + j;
}

void loop(int * i, int * j, int * r, int n, int (*operation)(int,int)) {
	for (int k=0; k<n; k++) {
		r[k] = operation(i[k],j[k]);
	}
}

void functieptrtest() {
	int r[4] = {};
	int i[4] = {0,1,2,3};
	int j[4] = {4,5,6,7};
	loop(i,j,r,4, &add);
	for (int k=0; k<4; k++) {
		cout << r[k] << endl;
	}

}

void arraytest() {
	int test[5] = {1,2,3,4,5};
	int copy[5];
	for (int i=0; i<5; i++) {
		copy[i] = test[i];
	}
	for (int i=0; i<5; i++) {
		cout << copy[i] << endl;
	}
}

void objecttest() {
	Foo instance(5);
	Foo test[3] = {instance, instance, instance};
	cout << test[1].get_i() << endl;
	test[2].set_i(10);
	cout << test[1].get_i() << endl;
}

void test1() {
	int input = 3;
	int boe[3] = {input%3,(input + 1)%3,(input + 2)%3};
	
	for (int i=0; i<3; i++) {
		cout << boe[i] << endl;
		//cout << i << endl;
	}
}

void mattest() {
	cout << "mattest:" << endl;

	char transA = 'N', transB = 'N';
	int M = 3, N = 3, K = 3;
	int Mz = 2, My = 2, Ky = 2;
	double alpha = 1, beta = 1;
	int lda = 3, ldb = 3, ldc = 3, ldz = 2, ldy = 2;

	// A = [1  2  0
	//      3  1  0
	//      0  0  0] = Z = Y
	double A[9] = {1,3,0,2,1,0,0,0,0};
	double Z[6] = {1,3,2,1,0,0};
	double Y[4] = {1,3,2,1};
	// B = [1  4  0
	//      0  1  3
	//      0  3  0]
	double B[9] = {1,0,0,4,1,3,0,3,0};
	// C_ij = A_ik * B_kj
	double C[9] = {0};

	// capture example:
	// auto printC = [&, nrows, ncols]()

	auto printMat = [](double * M, int nrows, int ncols) {
		for (int i=0; i<nrows; i++) {
			for (int j=0; j<ncols; j++) {
				cout << M[i + nrows * j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	};

	printMat(C,3,3);
	// works
	// dgemm_(transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
	// works
	// dgemm_(transA, transB, Mz, N, K, alpha, Z, ldz, B, ldb, beta, C, ldc);
	dgemm_(transA, transB, My, N, Ky, alpha, Y, ldy, B, ldb, beta, C, ldc);
	printMat(C,3,3);
}

// Result should be:
// mattest:
// 0 0 0 
// 0 0 0 
// 0 0 0 

// 1 6 6 
// 3 13 3 
// 0 0 0