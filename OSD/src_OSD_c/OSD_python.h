#ifndef OSD_python_h
#define OSD_pthon_h

// CAMBIAR EN FUNCIÓN DEL CÓDIGO UTILIZADO
// ****
int const N = 429;
int const M = 882;
// ****

class osd_decoder {


public:

    int H[N][M];

    int H_T_rect[M][N];

    int H_t[N][N];

    int H_lu[N][N];

    int H_lu_T[N][N];

    int H_sol[N][N];

    osd_decoder() = default;

    osd_decoder(int H_cod[N][M]) {

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                H[i][j] = H_cod[i][j];
            }
        }
    };


    int* decode(int synd[N], double prob_ini[M], int sol[M]) {

        int sorted_cols[M];
        sort_columns(prob_ini, sorted_cols);

        int used_cols[N];
        make_H(H, sorted_cols, H_sol, used_cols);

        solve_ecuationSys(synd, H_sol, used_cols, sol);

        return sol;
    }


private:
    int sort_columns(double prob[M], int sol[M]) {

        // Indices
        for (int i = 0; i < M; ++i) {
            sol[i] = i;
        }

        int temp_indices[M];

        for (int width = 1; width < M; width *= 2) {
            for (int i = 0; i < M; i += 2 * width) {

                int left = i;
                int right = (i + width < M) ? (i + width) : M;
                int end = (i + 2 * width < M) ? (i + 2 * width) : M;

                int l = left;
                int r = right;
                int t = left;

                while (l < right && r < end) {
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

                while (l < right) {
                    temp_indices[t++] = sol[l++];
                }

                while (r < end) {
                    temp_indices[t++] = sol[r++];
                }

                for (int k = left; k < end; k++) {
                    sol[k] = temp_indices[k];
                }
            }
        }

        return 0;
    }


    int make_H(int H[N][M], int sort_columns[M], int H_sol[N][N], int used_cols[N]) {

        for (int i = 0; i < N; i++) {
            used_cols[i] = -1;
            for (int j = 0; j < N; j++) {
                H_lu[i][j] = 0;
                H_sol[i][j] = 0;
            }
        }

        transponer_matriz_rect(H, H_T_rect);

        int asigned_cols = 0;
        int i = 0;
        while (i < M && asigned_cols < N) {

            int pos_insert = busca_posicion_columna(H_T_rect[sort_columns[i]], H_lu, N);

            if (pos_insert < N) {
                // En caso de que sea una posicion valida, se coloca la columna
                used_cols[pos_insert] = sort_columns[i];
                asigned_cols += 1;
                for (int j = 0; j < N; j++) {
                    H_lu[j][pos_insert] = H_T_rect[sort_columns[i]][j];
                }
            }
            i++;
        }

        transponer_matriz(H_lu, H_lu_T);

        for (int col = 0; col < N; col++) {
            eliminacion_salida(H_lu_T[col], H_lu, col, N);
            for (int j = 0; j < N; j++) {
                H_sol[j][col] = H_lu_T[col][j];
            }
        }

        return 0;
    }


    void eliminacion_gaussiana(int col1[N], int col2[N], int elem_estud) {

        for (int i = 0; i < N; i++) {
            col2[i] = col1[i] ^ col2[i];
        }

        col2[elem_estud] = 1;
    }


    int busca_posicion_columna(int columna_colocar[N], int H_lu[N][N], int numDetec) {
        int j = 0;

        while (j < numDetec && (columna_colocar[j] == 0 || H_lu[j][j] == 1)) {
            if (columna_colocar[j] == 1 && H_lu[j][j] == 1) {
                transponer_matriz(H_lu, H_t);
                eliminacion_gaussiana(H_t[j], columna_colocar, j);
            }
            j += 1;
        }
        return j;
    }


    void eliminacion_salida(int columna_eliminar[N], int H_lu[N][N], int index_col, int numDetec) {

        for (int j = index_col + 1; j < numDetec; j++) {

            if (columna_eliminar[j] == 1 && H_lu[j][j] == 1) {
                transponer_matriz(H_lu, H_t);
                eliminacion_gaussiana(H_t[j], columna_eliminar, j);
            }
        }
    }


    int solve_ecuationSys(int synd[N], int H_s[N][N], int used_cols[N], int sol[M]) {
        int result[N];

        multiplicarMatrizVector(synd, H_s, result);

        for (int i = 0; i < M; i++) {
            sol[i] = 0;
        }

        for (int i = 0; i < N; i++) {
            if (used_cols[i] >= 0 && used_cols[i] < M) {
                sol[used_cols[i]] = result[i];
            }
        }

        return 0;
    }

    // ÁLGEBRA LINEAL **************************************************************************************************/
    void multiplicarMatrizVector(const int synd[N], const int H_s[N][N], int result[N]) {
        for (int i = 0; i < N; i++) {
            int sum = 0;
            for (int j = 0; j < N; j++) {
                sum ^= H_s[i][j] & synd[j];
            }
            result[i] = sum;
        }
    }

    void transponer_matriz_rect(int entrada[N][M], int salida[M][N]) {
        //#pragma HLS ARRAY_PARTITION variable=entrada complete dim=0
        //#pragma HLS ARRAY_PARTITION variable=salida complete dim=0

        // El pipeline le dice a HLS que solape las operaciones de los bucles
        //#pragma HLS PIPELINE II=1
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                salida[j][i] = entrada[i][j];
            }
        }
    }


    void transponer_matriz(int entrada[N][N], int salida[N][N]) {
        //#pragma HLS ARRAY_PARTITION variable=entrada complete dim=0
        //#pragma HLS ARRAY_PARTITION variable=salida complete dim=0

        // El pipeline le dice a HLS que solape las operaciones de los bucles
        //#pragma HLS PIPELINE II=1
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                salida[j][i] = entrada[i][j];
            }
        }
    }

};


class osd_decoder_opt {


public:

    int H_T_cod[M][N];

    int H_T_lu[N][N];

    int H_sol[N][N];

    osd_decoder_opt() = default;

    osd_decoder_opt(int H_T_cod[M][N]) {

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                this->H_T_cod[i][j] = H_T_cod[i][j];
            }
        }
    };


    int* decode(int synd[N], double prob_ini[M], int sol[M]) {

        int sorted_cols[M];
        sort_columns(prob_ini, sorted_cols);

        int used_cols[N];
        make_H(H_T_cod, sorted_cols, H_sol, used_cols);

        solve_ecuationSys(synd, H_sol, used_cols, sol);

        return sol;
    }


private:
    int sort_columns(double prob[M], int sol[M]) {

        // Indices
        for (int i = 0; i < M; ++i) {
            sol[i] = i;
        }

        int temp_indices[M];

        for (int width = 1; width < M; width *= 2) {
            for (int i = 0; i < M; i += 2 * width) {

                int left = i;
                int right = (i + width < M) ? (i + width) : M;
                int end = (i + 2 * width < M) ? (i + 2 * width) : M;

                int l = left;
                int r = right;
                int t = left;

                while (l < right && r < end) {
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

                while (l < right) {
                    temp_indices[t++] = sol[l++];
                }

                while (r < end) {
                    temp_indices[t++] = sol[r++];
                }

                for (int k = left; k < end; k++) {
                    sol[k] = temp_indices[k];
                }
            }
        }

        return 0;
    }


    int make_H(int H[M][N], int sort_columns[M], int H_sol[N][N], int used_cols[N]) {

        for (int i = 0; i < N; i++) {
            used_cols[i] = -1;
            for (int j = 0; j < N; j++) {
                H_T_lu[i][j] = 0;
                H_sol[i][j] = 0;
            }
        }

        int asigned_cols = 0;
        int i = 0;
        while (i < M && asigned_cols < N) {

            // Para no modificar H_T_cod
            int col_temp[N];
            for (int k = 0; k < N; k++) {
                col_temp[k] = H_T_cod[sort_columns[i]][k];
            }

            int pos_insert = busca_posicion_columna(H_T_cod[sort_columns[i]], H_T_lu, N);

            if (pos_insert < N) {
                // En caso de que sea una posicion valida, se coloca la columna
                used_cols[pos_insert] = sort_columns[i];
                asigned_cols += 1;
                for (int j = 0; j < N; j++) {
                    H_T_lu[pos_insert][j] = H_T_cod[sort_columns[i]][j];
                }
            }

            // Devolvemos al estado original H_T_cod
            for (int k = 0; k < N; k++) {
                H_T_cod[sort_columns[i]][k] = col_temp[k];
            }

            i++;
        }

        for (int col = 0; col < N; col++) {
            eliminacion_salida(H_T_lu[col], H_T_lu, col, N);
            for (int j = 0; j < N; j++) {
                H_sol[j][col] = H_T_lu[col][j];
            }
        }

        return 0;
    }


    void eliminacion_gaussiana(int col1[N], int col2[N], int elem_estud) {

        for (int i = 0; i < N; i++) {
            col2[i] = col1[i] ^ col2[i];
        }

        col2[elem_estud] = 1;
    }


    int busca_posicion_columna(int columna_colocar[N], int H_T_lu[N][N], int numDetec) {
        int j = 0;

        while (j < numDetec && (columna_colocar[j] == 0 || H_T_lu[j][j] == 1)) {
            if (columna_colocar[j] == 1 && H_T_lu[j][j] == 1) {
                eliminacion_gaussiana(H_T_lu[j], columna_colocar, j);
            }
            j += 1;
        }
        return j;
    }


    void eliminacion_salida(int columna_eliminar[N], int H_T_lu[N][N], int index_col, int numDetec) {

        for (int j = index_col + 1; j < numDetec; j++) {

            if (columna_eliminar[j] == 1 && H_T_lu[j][j] == 1) {
                eliminacion_gaussiana(H_T_lu[j], columna_eliminar, j);
            }
        }

        for (int i = 0; i < numDetec; i++) {
            H_T_lu[index_col][i] = columna_eliminar[i];
        }
    }


    int solve_ecuationSys(int synd[N], int H_s[N][N], int used_cols[N], int sol[M]) {
        int result[N];

        multiplicarMatrizVector(synd, H_s, result);

        for (int i = 0; i < M; i++) {
            sol[i] = 0;
        }

        for (int i = 0; i < N; i++) {
            if (used_cols[i] >= 0 && used_cols[i] < M) {
                sol[used_cols[i]] = result[i];
            }
        }

        return 0;
    }

    // ÁLGEBRA LINEAL **************************************************************************************************/
    void multiplicarMatrizVector(const int synd[N], const int H_s[N][N], int result[N]) {
        for (int i = 0; i < N; i++) {
            int sum = 0;
            for (int j = 0; j < N; j++) {
                sum ^= H_s[i][j] & synd[j];
            }
            result[i] = sum;
        }
    }

    void transponer_matriz_rect(int entrada[N][M], int salida[M][N]) {
        //#pragma HLS ARRAY_PARTITION variable=entrada complete dim=0
        //#pragma HLS ARRAY_PARTITION variable=salida complete dim=0

        // El pipeline le dice a HLS que solape las operaciones de los bucles
        //#pragma HLS PIPELINE II=1
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                salida[j][i] = entrada[i][j];
            }
        }
    }


    void transponer_matriz(int entrada[N][N], int salida[N][N]) {
        //#pragma HLS ARRAY_PARTITION variable=entrada complete dim=0
        //#pragma HLS ARRAY_PARTITION variable=salida complete dim=0

        // El pipeline le dice a HLS que solape las operaciones de los bucles
        //#pragma HLS PIPELINE II=1
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                salida[j][i] = entrada[i][j];
            }
        }
    }

};
#endif