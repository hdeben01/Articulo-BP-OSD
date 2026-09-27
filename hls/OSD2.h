#ifndef OSD_H
#define OSD_H

#include <ap_fixed.h>
#include <ap_int.h>

int const N = 429;
int const M = 882;
// Es necesario añadir padding al valor de M para el bitonic sort, que necesita
// valores que sean potencia de 2
int const M_PADDED = 1024;
int const INF = 2e9;

typedef ap_uint<1> bit_t;
typedef double value_t;
typedef ap_uint<N> col_vec_t;

// Declaraciones de funciones OSD
void decode(bit_t synd_in[N], value_t prob_ini_in[M], bit_t sol[M]);
int sort_columns_bitonic(value_t prob[M_PADDED], int sol[M_PADDED]);

int make_H_vec(const col_vec_t H_T_cod[M], const int sorted_cols[M_PADDED], col_vec_t H_T_lu[N],int used_cols[N]);
int busca_posicion_columna_vec(col_vec_t &columna_colocar, const col_vec_t H_lu[N], int numDetec);
void eliminacion_salida_vec(col_vec_t &columna_eliminar, col_vec_t H_lu[N], int index_col, int numDetec);
int solve_ecuationSys_vec(col_vec_t synd, const col_vec_t H_s[N], const int used_cols[N], bit_t sol[M]);

#endif // OSD_H