#include "../include/min_sum_csc.h"



void min_sum(sparse_matrix_t *L, 
            int* syndrome, int size_checks, int size_vnode, 
            double Lj[VNODES], double alpha, int num_it, int *error_computed)
{

    for(int i = 0; i < num_it; i++){

        //color_printf(CYAN, "Iteration %d\n", i+1);

        double sum[VNODES];
        //int error[VNODES];

        compute_row_operations(L, syndrome, size_checks, size_vnode);
        
        /*
        printf("\tvalues after row ops:\n");
        for(int i = 0; i < L->nnz; i++){
            printf("%lf",L->values_csr[i]);
        }
        */

        csr_to_csc(L);

        compute_col_operations(L, syndrome, size_checks, size_vnode, alpha, Lj, sum);
        /*
        printf("\tvalues after col ops:\n");
        for(int i = 0; i < L->nnz; i++){
            printf("%lf",L->values_csc[i]);
        }
        */


        csc_to_csr(L);


        // Correct syndrome from the values of the codeword
        // if >= 0 then bit j = 0
        // if < 0  then bit j = 1 
        for(int j = 0; j < VNODES; j++){
            if (sum[j] >= 0) error_computed[j] = 0;
            else error_computed[j] = 1;
        }

        // ----------- DEBUG PRINT --------------
        /*
        printf("\tError computed: ");
        for(int j = 0; j < VNODES; j++){
            printf("%d ",error_computed[j]);
        }
        printf("\n");
        */

        /*
        printf("\tSum computed: ");
        for(int j = 0; j < VNODES; j++){
            printf("%lf ",sum[j]);
        }
        printf("\n");
        */

        /*
        printf("\tcsr values computed: ");
        for(int j = 0; j < L->nnz; j++){
            printf("%lf ",L->values_csr[j]);
        }
        printf("\n");
        */
        //------------------------------------------

        // Compute S = eH^T
        
        int resulting_syndrome[CHECK];
        for(int z = 0; z < CHECK; z++){
            int row_op = 0;
            int start_row = L->offset_rows[z];
            int end_row = L->offset_rows[z + 1];
            // To compute the error we only need the values where there is 1 in the H matrix
            for(int j = start_row; j < end_row; j++){
                int k = L->col_index[j];
                row_op ^= error_computed[k]; //we dont need the pcm because we are already iterating through the positions that are 1
            }
            resulting_syndrome[z] = row_op;
        }
        int error_found = 1;

        for(int i = 0; i < CHECK; i++){
            if(resulting_syndrome[i] != syndrome[i]) error_found = 0;
        }

        if(error_found) {
            
            break;
        }
        // ----------- DEBUG PRINT --------------
        /*
        for(int i = 0; i < CHECK; i++){
          
            printf("%d", syndrome[i]);
            
        }
        printf("\n");
        */

        /*
        printf("\tResulting syndrome: ");
        for(int i = 0; i < CHECK; i++){
            printf("%d ",resulting_syndrome[i]);
            if(resulting_syndrome[i] != syndrome[i]) error_found = 0;
        }
        printf("\n");

        if(error_found) {
            color_printf(GREEN, "\tERROR FOUND\n");
            break;
        }
        else if (i == num_it - 1) color_printf(RED, "\nUSED ALL ITERATIONS WITHOUT FINDING THE ERROR");

        printf("\n");
        */
        //------------------------------------------
    }
}

void compute_row_operations(sparse_matrix_t *L, 
                            int* syndrome, int size_checks, int size_vnode)
{

    for(int i = 0; i < CHECK; i++){
        if(i == size_checks) break;

        double min1 = DBL_MAX, min2 = DBL_MAX;
        int minpos = -1;
        int sign_minpos = 0;
        int row_sign = 0;

        // Search min1 and min2
        int start_row = L->offset_rows[i];
        int end_row = L->offset_rows[i + 1];
        for(int j = start_row; j < end_row; j++){
            //if(j == size_vnode) break;

            double val = L->values_csr[j];
            double abs_val = fabs(val);

            
            if(abs_val < min1){
                min2 = min1;
                min1 = abs_val;
                minpos = j;
                sign_minpos = (val >= 0 ? 0 : 1);
            }else if (abs_val < min2) {
                min2 = abs_val;
            }
            

            row_sign = row_sign ^ (val >= 0 ? 0 : 1);
        }

        for(int j = start_row; j < end_row; j++){
            //if(j == size_vnode) break;

            double val = L->values_csr[j];

            // Sign is negative (-1.0f) if the final signbit (operation in parethesis) is 0, 
            // and positive (1.0f) if its 1
            L->values_csr[j] = (1.0f - (2.0f * (row_sign ^ (val >= 0 ? 0 : 1) ^ syndrome[i]))) * min1;
            
        }

        // Assigning min2 to minpos only if we have atleast 2 elements non zero on that row
        if(end_row - start_row > 1)
            L->values_csr[minpos] = (1.0f - 2.0f * (row_sign ^ sign_minpos ^ syndrome[i])) * min2;
    }
}

void compute_col_operations(sparse_matrix_t *L,
                            int* syndrome, int size_checks, int size_vnode, double alpha, 
                            double Lj[CHECK], double sum[VNODES])
{
    
    for (int j = 0; j < VNODES; j++){
        if (j == size_vnode) break;

        // Possible optimization: Read entire column L[][j] to another variable beforehand and then add the values
        double sum_aux = 0.0f;
        int start_col = L->offset_cols[j];
        int end_col = L->offset_cols[j + 1];
        for(int i = start_col; i < end_col; i++){
            //if (i == size_checks) break;

            sum_aux += L->values_csc[i];

        
        }
        if(end_col - start_col > 0)
            sum[j] = Lj[j] + (alpha * sum_aux);
    }

    // Columnn iteration
    for (int j = 0; j < VNODES; j++){
        //if (j == size_vnode) break;

        // Possible optimization: Read entire column L[][j] to another variable beforehand and then add the values
        int start_col = L->offset_cols[j];
        int end_col = L->offset_cols[j + 1];
        for(int i = start_col; i < end_col; i++){
            //if (i == size_checks) break;

            L->values_csc[i] = sum[j] - (alpha * L->values_csc[i]);
        }
    
    }

}

int main() {
    double *L;//[CHECK][VNODES];
    L = (double*)malloc(CHECK*VNODES*sizeof(double));
    int *pcm_matrix;//[CHECK][VNODES];
    pcm_matrix = (int*)malloc(CHECK*VNODES*sizeof(int));
    double Lj[VNODES];
    FILE *file = fopen("input3.txt","r");
    if (file == NULL){
        perror("Error opening file");
        return 1;
    }

    // Read the probability p of the error model
    double p;
    if (fscanf(file, "%lf", &p) != 1) {
            fprintf(stderr, "Error reading probability p\n");
            fclose(file);
            return 1;
    }
    p = 1;//log((1.0f - p) /p);
    

    // Read rows and cols
    int rows,cols;
    if (fscanf(file, "%d %d", &rows, &cols) != 2) {
        fprintf(stderr, "Error reading dimensions\n");
        fclose(file);
        return 1;
    }

    // Read matrix L (initial beliefs)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int value;
            if (fscanf(file, "%d", &value) != 1) {
                fprintf(stderr, "Error reading valor of row %d col %d\n", i, j);
                fclose(file);
                return 1;
            }
            pcm_matrix[i * VNODES + j] = value;
            if(value == 1){
                L[i * VNODES + j] = 1.0f;
            }else{
                L[i * VNODES + j] = 0;
            }  
        }
    }

    // Read syndrome
    int syndrome[CHECK];
    for (int j = 0; j < rows; j++) {
        if (fscanf(file, "%d", &syndrome[j]) != 1) {
            fprintf(stderr, "Error reading syndrome[%d]\n", j);
            fclose(file);
            return 1;
        }
    }
    printf("Initial Syndrome:");
    for(int i = 0; i < CHECK; i++) printf(" %d", syndrome[i]);
    printf("\n");

    // Initialize Lj
    for (int j = 0; j < cols; j++) {
        Lj[j] = p;
    }
    printf("Lj: %.2f\n", p);

    sparse_matrix_t *L_sparse = malloc(sizeof(sparse_matrix_t));
    to_sparse_matrix_t(L,L_sparse, pcm_matrix);

    // Read alpha
    double alpha;
    if (fscanf(file, "%lf", &alpha) != 1) {
        fprintf(stderr, "Error reading alpha\n");
        fclose(file);
        return 1;
    }
    printf("Alpha: %.2f\n", alpha);

    int num_it = 10;
    if (fscanf(file, "%d", &num_it) != 1) {
        fprintf(stderr, "Error reading alpha\n");
        fclose(file);
        return 1;
    }
    //printf("Max Iterations: %d\n\n", num_it);

    //printf("Initial L Matrix:\n");
    //show_matrix(L, pcm_matrix, rows, cols);

    int error_computed[VNODES];

    min_sum(L_sparse, syndrome, rows, cols, Lj, alpha, num_it, &error_computed[0]);
    
    return 0;
}

void csr_to_csc(sparse_matrix_t *L){
    for(int i = 0; i < L->nnz; i++){
        L->values_csc[i] = L->values_csr[L->edges[i]];
    }
}

void csc_to_csr(sparse_matrix_t *L){
    for(int i = 0; i < L->nnz; i++){
        L->values_csr[L->edges[i]] = L->values_csc[i]; 
    }
}


// Recieves a flattened dense float matrix L (CHECK x VNODES) and fills out with CSR + edges for CSC
void to_sparse_matrix_t(double *L, sparse_matrix_t *out, int *pcm) {
    // Assumes out is already allocated and zeroed
    int nnz = 0;

    // Count non-zeros
    for (int i = 0; i < CHECK; i++) {
        for (int j = 0; j < VNODES; j++) {
           
            if (pcm[i * VNODES + j] != 0) nnz++;
        }
    }

    // Allocate arrays
    out->nnz = nnz;
    out->values_csr = (double*)malloc(nnz * sizeof(double));
    out->values_csc = (double*)malloc(nnz * sizeof(double));
    out->col_index = (int*)malloc(nnz * sizeof(int));
    out->row_index = (int*)malloc(nnz * sizeof(int));
    out->offset_rows = (int*)malloc((CHECK + 1) * sizeof(int));
    out->offset_cols = (int*)malloc((VNODES + 1) * sizeof(int));
    out->edges = (int*)malloc(nnz * sizeof(int));

    int *col_counts = (int*)calloc(VNODES, sizeof(int));

    // Fill CSR (row-wise)
    int idx = 0;
    for (int i = 0; i < CHECK; i++) {
        out->offset_rows[i] = idx;
        for (int j = 0; j < VNODES; j++) {
            double val = pcm[i * VNODES + j];
            if (val != 0.0f) {
                out->values_csr[idx] = 0;
                out->col_index[idx] = j;
                idx++;
            }
        }
    }
    out->offset_rows[CHECK] = idx;

    // Fill CSC (col-wise)
    idx = 0;
    for (int i = 0; i < VNODES; i++) {
        out->offset_cols[i] = idx;
        for (int j = 0; j < CHECK; j++) {
            double val = pcm[j * VNODES + i];
            if (val != 0.0f) {
                out->values_csc[idx] = 0;
                out->row_index[idx] = j;
                idx++;
            }
        }
    }
    out->offset_cols[VNODES] = idx;

    for(int i = 0; i < CHECK; i++){
        for(int k = out->offset_rows[i]; k < out->offset_rows[i+1]; k++){
            int col = out->col_index[k];
            int pos = out->offset_cols[col] + col_counts[col];
            out->edges[pos] = k;
            col_counts[col]++;
        }
    }

    /*
    for(int i = 0; i < nnz; i++){
        printf("%d ", out->edges[i]);
    }
    printf("\n");
    */
    free(col_counts);
}

/*
void show_matrix( double *matrix, int *non_zero,
                  int rows,  int cols)
{
    for(int i = 0; i < rows; i++){
        printf("\t");
        for(int j = 0; j < cols; j++){
            if(signbit(matrix[i * VNODES + j])) {
                if (non_zero[i * VNODES + j]) color_printf(YELLOW, " %.6f", matrix[i * VNODES + j]);
                else printf(" %.6f", matrix[i * VNODES + j]);
            }
            else {
                if (non_zero[i * VNODES + j]) color_printf(YELLOW, "  %.6f", matrix[i * VNODES + j]);
                else printf("  %.6f", matrix[i * VNODES + j]);
            }
        }
        printf("\n");
    }
    printf("\n");
}
*/
