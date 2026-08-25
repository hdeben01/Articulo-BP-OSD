#include <ap_fixed.h>
#ifndef OSD_H
#define OSD_H

int const N = 429;
int const M = 882;

typedef ap_uint<1> bit_t;
typedef double value_t;
// Declaraciones de funciones OSD
void decode(bit_t synd_in[N], value_t prob_ini_in[M], bit_t sol[M]);
int sort_columns(value_t prob[M], int sol[M]);
int make_H(bit_t H[N][M], int sorted_cols[M], bit_t H_sol[N][N], int used_cols[N]);
void eliminacion_gaussiana(bit_t col1[N], bit_t col2[N], int elem_estud);
int busca_posicion_columna(bit_t columna_colocar[N], bit_t H_lu[N][N], int numDetec);
void eliminacion_salida(bit_t columna_eliminar[N], bit_t H_lu[N][N], int index_col, int numDetec);
int solve_ecuationSys(bit_t synd[N], bit_t H_s[N][N], int used_cols[N], bit_t sol[M]);

// Álgebra lineal
void multiplicarMatrizVector(const bit_t synd[N], const bit_t H_s[N][N], bit_t result[N]);
void transponer_matriz_rect(bit_t entrada[N][M], bit_t salida[M][N]);
void transponer_matriz(bit_t entrada[N][N], bit_t salida[N][N]);

#endif // OSD_H