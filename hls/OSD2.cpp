#include "OSD.h"
#include <cmath>

inline col_vec_t read_synd_vec(const bit_t synd_in[N]){
    #pragma HLS INLINE

    col_vec_t s = 0;

    READ_SYND: for(int i = 0; i < N;i++){
        #pragma HLS PIPELINE II=1
        s[i] = synd_in[i];
    }

    return s;
}

void read_prob_ini(value_t prob_ini_in[M], value_t prob_ini[M_PADDED]){
    for(int i = 0; i < M; i++){
        #pragma HLS PIPELINE II=1
        prob_ini[i] = prob_ini_in[i];
    }

    for(int i = M; i < M_PADDED; i++){
        #pragma HLS PIPELINE II=1
        prob_ini[i] = INF;
    }
}

void write_sol(bit_t sol_aux[N], bit_t sol[N]){
    for(int i = 0; i < N;i++){
        #pragma HLS PIPELINE II=1
        sol[i] = sol_aux[i];
    }
}

void decode(bit_t synd_in[N], value_t prob_ini_in[M], bit_t sol[M]) {

#pragma HLS INTERFACE mode=m_axi bundle=gmem0 port=synd_in     depth=N
#pragma HLS INTERFACE mode=m_axi bundle=gmem1 port=prob_ini_in depth=M
#pragma HLS INTERFACE mode=m_axi bundle=gmem2 port=sol         depth=M

    value_t prob_ini[M_PADDED];
    static const bit_t H[N][M] = {
        #include "H1.h"
    };
    col_vec_t H_T_cod[M];

    int sorted_cols[M_PADDED];
    int used_cols[N];
    bit_t sol_aux[M];
    col_vec_t H_T_lu[N];

#pragma HLS BIND_STORAGE variable=H_T_lu type=RAM_T2P impl=BRAM
#pragma HLS ARRAY_PARTITION variable=H_T_lu      complete dim=1

#pragma HLS BIND_STORAGE variable=H_T_cod type=RAM_1P impl=BRAM
// Para forzar a que se guarden en BRAM y no gaste todas las LUTS

#pragma HLS ARRAY_PARTITION variable=H           complete dim=1
#pragma HLS ARRAY_PARTITION variable=prob_ini    cyclic factor=16
#pragma HLS ARRAY_PARTITION variable=sorted_cols cyclic factor=16
#pragma HLS ARRAY_PARTITION variable=used_cols   complete dim=1
#pragma HLS ARRAY_PARTITION variable=sol_aux     complete dim=1

    TRANSP_H: for (int j = 0; j < M; j++) {
        #pragma HLS PIPELINE II=1
        col_vec_t col_words = 0;
        for (int i = 0; i < N; i++) {
            #pragma HLS UNROLL
            col_words[i] = H[i][j];
        }
        H_T_cod[j] = col_words;
    }

    read_prob_ini(prob_ini_in, prob_ini);
    col_vec_t synd = read_synd_vec(synd_in);

    sort_columns_bitonic(prob_ini, sorted_cols);
    make_H_vec(H_T_cod, sorted_cols, H_T_lu, used_cols);
    solve_ecuationSys_vec(synd, H_T_lu, used_cols, sol);
}

int sort_columns_bitonic(value_t prob[M_PADDED], int sol[M_PADDED]) {
    #pragma HLS INLINE

    INIT_INDICES: for (int i = 0; i < M_PADDED; ++i) {
        #pragma HLS UNROLL
        sol[i] = i;
    }

    BITONIC_STAGE: for (int k = 2; k <= M_PADDED; k <<= 1) {
        #pragma HLS LOOP_FLATTEN off
        BITONIC_STEP: for (int j = k >> 1; j > 0; j >>= 1) {
            #pragma HLS LOOP_FLATTEN off
            BITONIC_CAS: for (int i = 0; i < M_PADDED; i++) {
                #pragma HLS UNROLL factor=2

                #pragma HLS DEPENDENCE variable=sol inter false
                int ixj = i ^ j; 

                if (ixj > i) {
                    bool dir = ((i & k) == 0);

                    int idx1 = sol[i];
                    int idx2 = sol[ixj];

                    value_t val1 = prob[idx1];
                    value_t val2 = prob[idx2];

                    bool swap_needed = false;
                    if (val1 == val2) {
                        swap_needed = dir ? (idx1 > idx2) : (idx1 < idx2);
                    } else {
                        swap_needed = dir ? (val1 > val2) : (val1 < val2);
                    }

                    if (swap_needed) {
                        sol[i] = idx2;
                        sol[ixj] = idx1;
                    }
                }
            }
        }
    }

    return 0;
}

int make_H_vec(const col_vec_t H_T_cod[M], const int sorted_cols[M_PADDED], col_vec_t H_T_lu[N], int used_cols[N]) {

    INIT_H: for(int i = 0; i < N;i++){
        #pragma HLS UNROLL
        used_cols[i] = -1;
        H_T_lu[i] = 0;
    }

    int asigned_cols = 0;
    int i = 0;

    FIND_COLS: while(i < M && asigned_cols < N) {
        #pragma HLS LOOP_FLATTEN off
        col_vec_t col_curr = H_T_cod[sorted_cols[i]];
        int pos_insert = busca_posicion_columna_vec(col_curr, H_T_lu, N);

        if (pos_insert < N){
            used_cols[pos_insert] = sorted_cols[i];
            asigned_cols++;
            H_T_lu[pos_insert] = col_curr;
        }
        i++;
    }


    ELIM_OUTER_LOOP: for (int col = 0; col < N; col++) {
        #pragma HLS PIPELINE II=1

        #pragma HLS DEPENDENCE variable=H_T_lu type=inter false

        col_vec_t col_elim = H_T_lu[col];
        eliminacion_salida_vec(col_elim, H_T_lu, col, N);
        H_T_lu[col] = col_elim;
    }

    return 0;
}


inline void eliminacion_gaussiana_vec(col_vec_t col1, col_vec_t &col2, int elem_estud){
    #pragma HLS INLINE
    col2 ^= col1;
    col2[elem_estud] = 1;
}

int busca_posicion_columna_vec(col_vec_t &columna_colocar, const col_vec_t H_lu[N], int numDetec){

    int j = 0;
    SEARCH_POS_WHILE: while(j < numDetec && (columna_colocar[j] == 0 || H_lu[j][j] == 1)){
        #pragma HLS PIPELINE II=1
        if(columna_colocar[j] == 1 && H_lu[j][j] == 1){
            eliminacion_gaussiana_vec(H_lu[j], columna_colocar, j);
        }
        j++;
    }
    return j;
}

void eliminacion_salida_vec(col_vec_t &columna_eliminar, col_vec_t H_lu[N], int index_col, int numDetec){
    OUTPUT_ELIM: for(int j = 0; j < N; j++){
        #pragma HLS UNROLL
        if (j >= index_col + 1 && j < numDetec){
            if (columna_eliminar[j] == 1 && H_lu[j][j] == 1){
                eliminacion_gaussiana_vec(H_lu[j], columna_eliminar, j);
            }
        }
    }
}

inline void multiplicarMatrizVector_vec(const col_vec_t synd, const col_vec_t H_s[N], bit_t result[N]){
    #pragma HLS INLINE

    col_vec_t accum = 0;

    MULT_MAT_VEC_ACCUM: for(int j = 0; j < N; j++){
        #pragma HLS PIPELINE II=1
        if (synd[j] == 1) {
            accum ^= H_s[j];
        }
    }

    WRITE_RESULT: for(int i = 0; i < N; i++){
        #pragma HLS UNROLL
        result[i] = accum[i];
    }
}

int solve_ecuationSys_vec(const col_vec_t synd, const col_vec_t H_s[N], const int used_cols[N], bit_t sol_aux[M]) {
    bit_t result[N];

    #pragma HLS ARRAY_PARTITION variable=result complete dim=1

    multiplicarMatrizVector_vec(synd, H_s, result);

    INIT_SOL: for (int i = 0; i < M; i++) {
        #pragma HLS PIPELINE II=1
        sol_aux[i] = 0;
    }

    MAP_SOL: for (int i = 0; i < N; i++) {
        #pragma HLS PIPELINE II=1
        int idx = used_cols[i]; // Para no acceder dos veces
        if (idx >= 0 && idx < M) {
            sol_aux[idx] = result[i];
        }
    }

    return 0;
}