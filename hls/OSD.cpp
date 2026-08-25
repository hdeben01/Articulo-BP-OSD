#include "OSD.h"

void read_synd(bit_t synd_in[N], bit_t synd[N]){
    for(int i = 0; i < N; i++){
        synd[i] = synd_in[i];
    }
}

void read_prob_ini(value_t prob_ini_in[M], value_t prob_ini[M]){
    for(int i = 0; i < M; i++){
        prob_ini[i] = prob_ini_in[i];
    }
}

void decode(bit_t synd_in[N], value_t prob_ini_in[M], bit_t sol[M]) {

    bit_t synd[N];
    value_t prob_ini[M];
    bit_t H_sol[N][N];
    bit_t H[N][M] = {
        #include "H1.h"
    };
/*#pragma HLS ARRAY_PARTITION variable=H_sol complete dim=2
#pragma HLS ARRAY_PARTITION variable=synd complete dim=0
#pragma HLS ARRAY_PARTITION variable=prob_ini complete dim=0
#pragma HLS ARRAY_PARTITION variable=H complete dim=2*/

    read_synd(synd_in, synd);
    read_prob_ini(prob_ini_in, prob_ini);

    int sorted_cols[M];
    sort_columns(prob_ini, sorted_cols);

    int used_cols[N];
    make_H(H, sorted_cols, H_sol, used_cols);

    solve_ecuationSys(synd, H_sol, used_cols, sol);

    //return sol;
}

int sort_columns(value_t prob[M], int sol[M]) {

    // Indices
    INIT_INDICES: for (int i = 0; i < M; ++i) {
        sol[i] = i;
    }

    int temp_indices[M];

    SORT_OUTER_LOOP: for (int width = 1; width < M; width *= 2) {
        SORT_INNER_LOOP: for (int i = 0; i < M; i += 2 * width) {

            int left = i;
            int right = (i + width < M) ? (i + width) : M;
            int end = (i + 2 * width < M) ? (i + 2 * width) : M;

            int l = left;
            int r = right;
            int t = left;

            SORT_WHILE_1: while (l < right && r < end) {
                /* ORDEN ASCENDENTE (LLR) */
                if (prob[sol[l]] <= prob[sol[r]]) {
                    temp_indices[t++] = sol[l++];
                }
                else {
                    temp_indices[t++] = sol[r++];
                }
                //*/

                /* ORDEN DESCENDENTE (probabilidades) 
                if (prob[sol[l]] >= prob[sol[r]]) {
                    temp_indices[t++] = sol[l++];
                }
                else {
                    temp_indices[t++] = sol[r++];
                }
                //*/
            }

            SORT_WHILE_2: while (l < right) {
                temp_indices[t++] = sol[l++];
            }

            SORT_WHILE_3: while (r < end) {
                temp_indices[t++] = sol[r++];
            }

            UPDATE_SOL: for (int k = left; k < end; k++) {
                sol[k] = temp_indices[k];
            }
        }
    }

    return 0;
}

int make_H(bit_t H[N][M], int sorted_cols[M], bit_t H_sol[N][N], int used_cols[N]) {

    bit_t H_T_rect[M][N];
    bit_t H_lu[N][N];
    bit_t H_lu_T[N][N];

/*#pragma HLS ARRAY_PARTITION variable=H_T_rect complete dim=2
#pragma HLS ARRAY_PARTITION variable=H_lu complete dim=2
#pragma HLS ARRAY_PARTITION variable=H_lu_T complete dim=2*/


    INIT_H_OUTER: for (int i = 0; i < N; i++) {
        used_cols[i] = -1;
//#pragma HLS UNROLL 
        INIT_H_INNER: for (int j = 0; j < N; j++) {
            H_lu[i][j] = 0;
            H_sol[i][j] = 0;
        }
    }

    transponer_matriz_rect(H, H_T_rect);

    int asigned_cols = 0;
    int i = 0;
    FIND_COLS_WHILE: while (i < M && asigned_cols < N) {

        int pos_insert = busca_posicion_columna(H_T_rect[sorted_cols[i]], H_lu, N);

        if (pos_insert < N) {
            // En caso de que sea una posicion valida, se coloca la columna
            used_cols[pos_insert] = sorted_cols[i];
            asigned_cols += 1;
            COPY_COL_LOOP: for (int j = 0; j < N; j++) {
                H_lu[j][pos_insert] = H_T_rect[sorted_cols[i]][j];
            }
        }
        i++;
    }

    transponer_matriz(H_lu, H_lu_T);

    ELIM_OUTER_LOOP: for (int col = 0; col < N; col++) {
        eliminacion_salida(H_lu_T[col], H_lu, col, N);
//#pragma HLS UNROLL 
        COPY_H_SOL_LOOP: for (int j = 0; j < N; j++) {
            H_sol[j][col] = H_lu_T[col][j];
        }
    }

    return 0;
}

void eliminacion_gaussiana(bit_t col1[N], bit_t col2[N], int elem_estud) {

    GAUSS_ELIM_LOOP: for (int i = 0; i < N; i++) {
        col2[i] = col1[i] ^ col2[i];
    }

    col2[elem_estud] = 1;
}

int busca_posicion_columna(bit_t columna_colocar[N], bit_t H_lu[N][N], int numDetec) {
    int j = 0;
    bit_t H_t[N][N];

    SEARCH_POS_WHILE: while (j < numDetec && (columna_colocar[j] == 0 || H_lu[j][j] == 1)) {
        if (columna_colocar[j] == 1 && H_lu[j][j] == 1) {
            transponer_matriz(H_lu, H_t);
            eliminacion_gaussiana(H_t[j], columna_colocar, j);
        }
        j += 1;
    }
    return j;
}

void eliminacion_salida(bit_t columna_eliminar[N], bit_t H_lu[N][N], int index_col, int numDetec) {

    bit_t H_t[N][N];
//#pragma HLS ARRAY_PARTITION variable=H_t complete dim=2
    OUTPUT_ELIM_LOOP: for (int j = index_col + 1; j < numDetec; j++) {

        if (columna_eliminar[j] == 1 && H_lu[j][j] == 1) {
            transponer_matriz(H_lu, H_t);
            eliminacion_gaussiana(H_t[j], columna_eliminar, j);
        }
    }
}

int solve_ecuationSys(bit_t synd[N], bit_t H_s[N][N], int used_cols[N], bit_t sol[M]) {
    bit_t result[N];

    multiplicarMatrizVector(synd, H_s, result);

    INIT_SOL_LOOP: for (int i = 0; i < M; i++) {
        sol[i] = 0;
    }

    MAP_SOL_LOOP: for (int i = 0; i < N; i++) {
        if (used_cols[i] >= 0 && used_cols[i] < M) {
            sol[used_cols[i]] = result[i];
        }
    }

    return 0;
}

// ÁLGEBRA LINEAL **************************************************************************************************/
void multiplicarMatrizVector(const bit_t synd[N], const bit_t H_s[N][N], bit_t result[N]) {
    MULT_MAT_VEC_OUTER: for (int i = 0; i < N; i++) {
        int sum = 0;
        MULT_MAT_VEC_INNER: for (int j = 0; j < N; j++) {
            sum ^= H_s[i][j] & synd[j];
        }
        result[i] = sum;
    }
}

void transponer_matriz_rect(bit_t entrada[N][M], bit_t salida[M][N]) {
    //#pragma HLS ARRAY_PARTITION variable=entrada complete dim=0
    //#pragma HLS ARRAY_PARTITION variable=salida complete dim=0

    // El pipeline le dice a HLS que solape las operaciones de los bucles
    //#pragma HLS PIPELINE II=1
    TRANS_RECT_ROW: for (int i = 0; i < N; i++) {
//#pragma HLS PIPELINE II=1
        TRANS_RECT_COL: for (int j = 0; j < M; j++) {
            salida[j][i] = entrada[i][j];
        }
    }
}

void transponer_matriz(bit_t entrada[N][N], bit_t salida[N][N]) {
    //#pragma HLS ARRAY_PARTITION variable=entrada complete dim=0
    //#pragma HLS ARRAY_PARTITION variable=salida complete dim=0

    // El pipeline le dice a HLS que solape las operaciones de los bucles
//#pragma HLS PIPELINE II=1
    TRANS_MAT_ROW: for (int i = 0; i < N; i++) {
        TRANS_MAT_COL: for (int j = 0; j < N; j++) {
            salida[j][i] = entrada[i][j];
        }
    }
}