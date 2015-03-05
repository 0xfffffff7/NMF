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
void matrix_add(double *c, double *a, double *b, int n, int m)
{
	int	i = 0;
	int j = 0;
	
	for(i =0; i <n; i++){
		for(j = 0; j < m; j++){
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
int matrix_diff(double *a, double *b, int n, int m)
{
	int	i = 0;
	int j = 0;
	int diff = 0;
	
	for(i =0; i < n; i++){
		for(j = 0; j < m; j++){
			diff = (int)pow( *(a + i * m + j) - *(b + i * m + j) , 2);
		}
	}
	
	return diff;
}

//-------------------------------------------
// マトリックスの初期化.
// ランダムな値で初期化する.
//-------------------------------------------
void matrix_init(double *a, int n, int m)
{
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			*( a + (i * m + j) ) = (double) (rand() % 10 + 7);
			srand( (unsigned int)time(NULL) ^ i + j);
		}
	}
}

//-------------------------------------------
// マトリックスの初期化.
// ゼロ初期化.
//-------------------------------------------
void matrix_init_zero(double *a, int n, int m)
{
	memset(a, 0x00, n * m * sizeof(double) );
}

//-------------------------------------------
// マトリックスのコピー.
// a = コピーされる行列.
// b = コピーする行列.
//-------------------------------------------
void matrix_copy(double *a, const double *b, int n, int m)
{
	memcpy(a, b, n * m * sizeof(double) );
}

//-------------------------------------------
// マトリックスの表示.
//-------------------------------------------
void matrix_print(double *a, int n, int m)
{
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			printf("%f ", *( a + (i * m + j )) );
		}
		printf("\n");
	}
    
	printf("\n");
}

//-------------------------------------------
// マトリックスの表示.
// 転置版.
//-------------------------------------------
void matrix_print_t(double *a, int n, int m)
{
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			printf("%f ", *( a + (j * m + i )) );
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
void matrix_x(double *c, double *a, double *b, int l, int m, int n)
{
	double t = 0;
	for(int i = 0; i < l; i++){
		for(int j = 0; j < n; j++){
			for(int k = 0; k < m; k++){
				t += *( a + (i * m + k)) * *( b + (k * n + j) );
			}
			*( c + (i * n + j) ) = t;
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
void matrix_divide(double *c, double *a, double *b, int l, int m, int n)
{
	double t = 0;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < l; j++){
			for(int k = 0; k < m; k++){
				t += *( a + i * m + k) / *( b + k * n + j );
			}
            
			*( c + i * l + j ) = t;
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
void matrix_x(double *c, double *a, double *b, int m, int n)
{
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			*( c + (i * n + j) ) = *( a + (i * n + j)) * *( b + (i * n + j) );
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
void matrix_divide(double *c, double *a, double *b, int m, int n)
{
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			*( c + (i * n + j) ) = *( a + (i * n + j)) / *( b + (i * n + j) );
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
void matrix_transpose(double *a, double *b, int n, int m)
{
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
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
void matrix_left_transpose_x(double *c, double *a, double *b, int l, int m, int n)
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
void matrix_right_transpose_x(double *c, double *a, double *b, int l, int m, int n)
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
void matrix_tx_left(double *c, double *a, double *b, int l, int m, int n)
{
	double t = 0;
	for(int i = 0; i < l; i++){
		for(int j = 0; j < n; j++){
			for(int k = 0; k < m; k++){
				t += *( a + (k * l + i)) * *( b + (k * n + j) );
			}
			*( c + (i * n + j) ) = t;
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
void matrix_tx_right(double *c, double *a, double *b, int l, int m, int n)
{
	double t = 0;
	for(int i = 0; i < l; i++){
		for(int j = 0; j < n; j++){
			for(int k = 0; k < m; k++){
				t += *(a + (i * m + k) ) * *(b + (j * m + k) );
			}
			*( c + (i * n + j) ) = t;
			t = 0;
		}
	}
}



#define FEATURES	2
#define DOCUMENTS	5
#define DATAS		4

const double nARRAY[DOCUMENTS][DATAS]
= {	{ 50.0,50.0,50.0,50.0 },
    { 100.0,20.0,30.0,90.0 },
    { 100.0,30.0,40.0,100.0 },
    { 90.0,10.0,10.0,80.0 },
    { 10.0,100.0,90.0,20.0 }
};

void factorize(double *v, int raw, int col, int features, unsigned int count)
{
    
	// 重みの行列と特徴の行列をランダムな値で初期化.
	double h[FEATURES][DATAS];
	double hn[FEATURES][DATAS];
	double hd[FEATURES][DATAS];
	double hd_t[FEATURES][FEATURES];
	double w[DOCUMENTS][FEATURES];
	double wn[DOCUMENTS][FEATURES];
	double wd[DOCUMENTS][FEATURES];
	double wd_t[DOCUMENTS][DATAS];
    
	// 重みと特徴の乗算結果用行列.
	double wh[DOCUMENTS][DATAS];
	matrix_init_zero((double*)wh, DOCUMENTS, DATAS);
    
    
	matrix_init((double*)h, FEATURES, DATAS);
	matrix_init_zero((double*)hn, FEATURES, DATAS);
	matrix_init_zero((double*)hd, FEATURES, DATAS);
	matrix_init_zero((double*)hd_t, FEATURES, FEATURES);
	matrix_init((double*)w, DOCUMENTS, FEATURES);
	matrix_init_zero((double*)wn, DOCUMENTS, FEATURES);
	matrix_init_zero((double*)wd, DOCUMENTS, FEATURES);
	matrix_init_zero((double*)wd_t, DOCUMENTS, DATAS);
    
    
	matrix_print((double*)h, FEATURES, DATAS);
	matrix_print((double*)w, DOCUMENTS, FEATURES);
    
    
	unsigned int i = 0;
	for(i = 0; i < count; i++){
        
		v = (double*)nARRAY;
        
		// 特徴の重みの行列の積から、元データとの差を計算する.
		matrix_x((double*)wh, (double*)w, (double*)h, DOCUMENTS, FEATURES, DATAS);
		int cost = matrix_diff(v, (double*)wh, DOCUMENTS, DATAS);
		
		// 差が完全にゼロになったらループを抜ける.
		if(cost == 0){
			break;
		}
        
		//-------------------------------------------------------------------------------------------
		// 特徴の行列を更新する.
		//-------------------------------------------------------------------------------------------
		// hn　転置した重みの行列にデータ行列を掛け合わせたもの
		matrix_tx_left((double*)hn, (double*)w, v, FEATURES, DOCUMENTS, DATAS);
        
		// hd  転置した重みの行列に重みの行列を掛け合わせたものに特徴の行列を掛け合わせたもの
		matrix_tx_left((double*)hd_t, (double*)w, (double*)w, FEATURES, DOCUMENTS, FEATURES);
		matrix_x((double*)hd, (double*)hd_t, (double*)h, FEATURES, FEATURES, DATAS);
        
        
		matrix_x((double*)h, (double*)h, (double*)hn, FEATURES, DATAS);
		matrix_divide((double*)h, (double*)h, (double*)hd, FEATURES, DATAS);
        
        
        
		//-------------------------------------------------------------------------------------------
		// 重みの行列を更新する.
		//-------------------------------------------------------------------------------------------
        
		// wn  データ行列に転置した特長の行列を掛け合わせたもの.
		matrix_tx_right((double*)wn, v, (double*)h, DOCUMENTS, DATAS, FEATURES);
        
		// wd  重みの行列に、特徴の行列を掛け合わせたものに転置した特長の行列を掛け合わせたもの.
		matrix_x((double*)wd_t, (double*)w, (double*)h, DOCUMENTS, FEATURES, DATAS);
		matrix_tx_right((double*)wd, (double*)wd_t, (double*)h, DOCUMENTS, DATAS, FEATURES);
		
		matrix_x((double*)w, (double*)w, (double*)wn, DOCUMENTS, FEATURES);
		matrix_divide((double*)w, (double*)w, (double*)wd, DOCUMENTS, FEATURES);
        
	}
	
	matrix_print((double*)h, FEATURES, DATAS);
	matrix_print((double*)w, DOCUMENTS, FEATURES);
    
	matrix_x((double*)wh, (double*)w, (double*)h, DOCUMENTS, FEATURES, DATAS);
	matrix_print((double*)wh, DOCUMENTS, DATAS);
}

int main(int argc, char* argv[])
{
	factorize( (double*)nARRAY, 1,1,FEATURES, 5000);
    
	return 0;
}
