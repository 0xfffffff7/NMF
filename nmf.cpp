#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <stdlib.h>

//-------------------------------------------
// マトリックスの和.
//-------------------------------------------
// n = 行数.
// m = 列数.
// c = 結果のマトリックス.
// a,b = 和を計算する対象の二つのマトリックス.
//-------------------------------------------
void matrix_add(double *c, const double *a, const double *b, const int n, const int m)
{
	int	i = 0;
	int j = 0;

	for (i = 0; i <n; i++){
		for (j = 0; j < m; j++){
			*(c + i * m + j) = *(a + i * m + j) + *(b + i * m + j);
		}
	}
}

//-------------------------------------------
// マトリックスの差.
//-------------------------------------------
// n = 行数.
// m = 列数.
//-------------------------------------------
int matrix_diff(const double *a, const double *b, const int n, const int m)
{
	int	i = 0;
	int j = 0;
	int diff = 0;

	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			diff = (int)pow(*(a + i * m + j) - *(b + i * m + j), 2);
		}
	}

	return diff;
}

//-------------------------------------------
// マトリックスの初期化.
// ランダムな値で初期化する.
//-------------------------------------------
void matrix_init(double *a, const int n, const int m)
{
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
			*(a + (i * m + j)) = (double)(rand() % 1000);
			srand((unsigned int)time(NULL) ^ i + j);
		}
	}
}

//-------------------------------------------
// マトリックスの初期化.
// ゼロ初期化.
//-------------------------------------------
void matrix_init_zero(double *a, const int n, const int m)
{
	memset(a, 0x00, n * m * sizeof(double));
}

//-------------------------------------------
// マトリックスのコピー.
// a = コピーされる行列.
// b = コピーする行列.
//-------------------------------------------
void matrix_copy(double *a, const double *b, const int n, const int m)
{
	memcpy(a, b, n * m * sizeof(double));
}

//-------------------------------------------
// マトリックスの表示.
// coefficient = 表示する時に乗ずる値.
// bDisplayPoint = 小数点以下を表示するかしないか.
//-------------------------------------------
void matrix_print(const double *a, const int n, const int m, int coefficient = 1, bool bDisplayPoint = false)
{
	if (bDisplayPoint){
		for (int i = 0; i < n; i++){
			for (int j = 0; j < m; j++){
				printf("%f ", *(a + (i * m + j)) * coefficient);
			}
			printf("\n");
		}
	}
	else{
		for (int i = 0; i < n; i++){
			for (int j = 0; j < m; j++){
				printf("%d ", (int)(*(a + (i * m + j)) * coefficient));
			}
			printf("\n");
		}

	}

	printf("\n");
}

//-------------------------------------------
// マトリックスの表示.
// 転置版.
//-------------------------------------------
void matrix_print_t(const double *a, const int n, const int m)
{
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			printf("%f ", *(a + (j * m + i)));
		}
		printf("\n");
	}

	printf("\n");
}

//-------------------------------------------
// マトリックスの積.
//-------------------------------------------
// c = 結果行列の先頭アドレス.
// a,b = 割る行列の先頭アドレス.
// l = a の行数.
// m = a の列数、b の行数.
// n = b の列数.
//-------------------------------------------
void matrix_x(double *c, const double *a, const double *b, const int l, const int m, const int n)
{
	double t = 0;
	for (int i = 0; i < l; i++){
		for (int j = 0; j < n; j++){
			for (int k = 0; k < m; k++){
				t += *(a + (i * m + k)) * *(b + (k * n + j));
			}
			*(c + (i * n + j)) = t;
			t = 0;
		}
	}
}

//-------------------------------------------
// マトリックスの商.
//-------------------------------------------
// c = 結果行列の先頭アドレス.
// a,b = 割る行列の先頭アドレス.
// l = a の行数.
// m = a の列数、b の行数.
// n = b の列数.
//-------------------------------------------
void matrix_divide(double *c, const double *a, const double *b, const int l, const int m, const int n)
{
	double t = 0;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < l; j++){
			for (int k = 0; k < m; k++){
				t += *(a + i * m + k) / *(b + k * n + j);
			}

			*(c + i * l + j) = t;
			t = 0;
		}
	}
}

//-------------------------------------------
// マトリックスの積.
// 対応する値同士の乗算.
//-------------------------------------------
// c = 結果行列の先頭アドレス.
// a,b = 掛ける行列の先頭アドレス.
// m = 行数.
// n = 列数.
//-------------------------------------------
void matrix_x(double *c, const double *a, const double *b, const int m, const int n)
{
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			*(c + (i * n + j)) = *(a + (i * n + j)) * *(b + (i * n + j));
		}
	}
}

//-------------------------------------------
// マトリックスの商.
// 対応する値同士の除算.
//-------------------------------------------
// c = 結果行列の先頭アドレス.
// a,b = 割る行列の先頭アドレス.
// m = 行数.
// n = 列数.
//-------------------------------------------
void matrix_divide(double *c, const double *a, const double *b, const int m, const int n)
{
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			*(c + (i * n + j)) = *(a + (i * n + j)) / *(b + (i * n + j));
		}
	}
}

//-------------------------------------------
// マトリックスの転置.
//-------------------------------------------
// a = 転置する行列.
// b = 転置した結果の行列.
// n = aの行数.
// m = aの列数.
//-------------------------------------------
void matrix_transpose(const double *a, double *b, const int n, const int m)
{
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
			*(b + j * n + i) = *(a + i * m + j);
		}
	}
}

//-------------------------------------------
// マトリックスの転置と積.
// 左集合転置版.
//-------------------------------------------
// c = 結果行列の先頭アドレス.
// a,b = 掛ける行列の先頭アドレス.
// l = a の列数.
// m = a の行数、b の行数.
// n = b の列数.
//-------------------------------------------
void matrix_left_transpose_x(double *c, const double *a, const double *b, const int l, const int m, const int n)
{
	double *t = new double[l * m];
	matrix_init_zero(t, l, m);

	matrix_transpose(a, t, l, m);

	matrix_x(c, t, b, l, m, n);

	delete[] t;
	t = NULL;
}

//-------------------------------------------
// マトリックスの転置と積.
// 右集合転置版.
//-------------------------------------------
// c = 結果行列の先頭アドレス.
// a,b = 掛ける行列の先頭アドレス.
// l = a の行数.
// m = a の列数、b の列数.
// n = b の行数.
//-------------------------------------------
void matrix_right_transpose_x(double *c, const double *a, const double *b, const int l, const int m, const int n)
{
	double *t = new double[n * m];
	matrix_init_zero(t, n, m);

	matrix_transpose(b, t, n, m);

	matrix_x(c, a, t, l, m, n);

	delete[] t;
	t = NULL;
}



//-------------------------------------------
// マトリックスの転置と積.
// 左集合転置版.
// 転置しながら乗算を行う.
//-------------------------------------------
// c = 結果行列の先頭アドレス.
// a,b = 掛ける行列の先頭アドレス.
// l = a の列数.
// m = a の行数、b の行数.
// n = b の列数.
//-------------------------------------------
void matrix_tx_left(double *c, const double *a, const double *b, const int l, const int m, const int n)
{
	double t = 0;
	for (int i = 0; i < l; i++){
		for (int j = 0; j < n; j++){
			for (int k = 0; k < m; k++){
				t += *(a + (k * l + i)) * *(b + (k * n + j));
			}
			*(c + (i * n + j)) = t;
			t = 0;
		}
	}
}


//-------------------------------------------
// マトリックスの転置と積.
// 右集合転置版.
// 転置しながら乗算を行う.
//-------------------------------------------
// c = 結果行列の先頭アドレス.
// a,b = 掛ける行列の先頭アドレス.
// l = a の行数.
// m = a の列数、b の列数.
// n = b の行数.
//-------------------------------------------
void matrix_tx_right(double *c, const double *a, const double *b, const int l, const int m, const int n)
{
	double t = 0;
	for (int i = 0; i < l; i++){
		for (int j = 0; j < n; j++){
			for (int k = 0; k < m; k++){
				t += *(a + (i * m + k)) * *(b + (j * m + k));
			}
			*(c + (i * n + j)) = t;
			t = 0;
		}
	}
}

class MatrixFactor{
public:
	MatrixFactor(){
		h = NULL;
		w = NULL;
	}
	double* h;
	double* w;

	~MatrixFactor(){
		h ? NULL : delete h;
		w ? NULL : delete h;
	}
};

void factorize(double *v, int row, int col, int features, unsigned int count, MatrixFactor *mf)
{

	// 重みの行列と特徴の行列をランダムな値で初期化.
	double * h = new double[features * col];
	double * hn = new double[features * col];
	double * hd = new double[features * col];
	double * hd_t = new double[features * features];
	double * w = new double[row * features];
	double * wn = new double[row * features];
	double * wd = new double[row * features];
	double * wd_t = new double[row * col];

	// 重みと特徴の乗算結果用行列.
	double * wh = new double[row * col];

	matrix_init_zero((double*)wh, row, col);

	matrix_init(h, features, col);
	matrix_init_zero(hn, features, col);
	matrix_init_zero(hd, features, col);
	matrix_init_zero(hd_t, features, features);
	matrix_init(w, row, features);
	matrix_init_zero(wn, row, features);
	matrix_init_zero(wd, row, features);
	matrix_init_zero(wd_t, row, col);


	//matrix_print((double*)h, features, col);
	//matrix_print((double*)w, raw, features);

	unsigned int i = 0;
	for (i = 0; i < count; i++){

		// 特徴の重みの行列の積から、元データとの差を計算する.
		matrix_x(wh, w, h, row, features, col);
		int cost = matrix_diff(v, wh, row, col);

		// 差が完全にゼロになったらループを抜ける.
		if (cost == 0){
			break;
		}

		//-------------------------------------------------------------------------------------------
		// 特徴の行列を更新する.
		//-------------------------------------------------------------------------------------------
		// hn　転置した重みの行列にデータ行列を掛け合わせたもの
		matrix_tx_left(hn, w, v, features, row, col);

		// hd  転置した重みの行列に重みの行列を掛け合わせたものに特徴の行列を掛け合わせたもの
		matrix_tx_left(hd_t, w, w, features, row, features);
		matrix_x(hd, hd_t, h, features, features, col);


		matrix_x(h, h, hn, features, col);
		matrix_divide(h, h, hd, features, col);



		//-------------------------------------------------------------------------------------------
		// 重みの行列を更新する.
		//-------------------------------------------------------------------------------------------

		// wn  データ行列に転置した特長の行列を掛け合わせたもの.
		matrix_tx_right(wn, v, h, row, col, features);

		// wd  重みの行列に、特徴の行列を掛け合わせたものに転置した特長の行列を掛け合わせたもの.
		matrix_x(wd_t, w, h, row, features, col);
		matrix_tx_right(wd, wd_t, h, row, col, features);

		matrix_x(w, w, wn, row, features);
		matrix_divide(w, w, wd, row, features);

	}

	mf->h = h;
	mf->w = w;


	//matrix_x(wh, w, h, ROW, FEATURES, COL);
	//matrix_print(wh, ROW, COL);
	//matrix_print(w, ROW, FEATURES);
	//matrix_print(h, FEATURES, COL);

	delete[] hn;
	delete[] hd;
	delete[] hd_t;
	delete[] wn;
	delete[] wd;
	delete[] wd_t;
	delete[] wh;
}


int main(int argc, char* argv[])
{
	int FEATURES = 2;
	int ROW = 5;
	int COL = 4;
	int COUNT = 5000;

	const double nARRAY[5][4]
		= { { 50.0, 50.0, 50.0, 50.0 },
		{ 100.0, 20.0, 30.0, 90.0 },
		{ 100.0, 30.0, 40.0, 100.0 },
		{ 90.0, 10.0, 10.0, 80.0 },
		{ 10.0, 100.0, 90.0, 20.0 }
	};

	MatrixFactor mf;
	factorize((double*)nARRAY, ROW, COL, FEATURES, COUNT, &mf);

	matrix_print(mf.w, ROW, FEATURES);
	matrix_print(mf.h, FEATURES, COL, 1000);

	return 0;
}
