#include "../include/min_sum.h"



void min_sum(double *L,  int *pcm_matrix, 
                            int* syndrome, int size_checks, int size_vnode, 
                            double Lj[VNODES], double alpha, int num_it, int *error_computed)
{

    for(int i = 0; i < num_it; i++){

        //color_printf(CYAN, "Iteration %d\n", i+1);


        double sum[VNODES];
        //int error[VNODES];

        compute_row_operations(L, pcm_matrix, syndrome, size_checks, size_vnode);
        //printf("\tL matrix after row ops:\n");
        //show_matrix(L, pcm_matrix, size_checks, size_vnode);

        compute_col_operations(L, pcm_matrix, syndrome, size_checks, size_vnode, alpha, Lj, sum);
        //printf("\tL matrix after col ops:\n");
        //show_matrix(L, pcm_matrix, size_checks, size_vnode);

        // Correct syndrome from the values of the codeword
        // if >= 0 then bit j = 0
        // if < 0  then bit j = 1 
        for(int j = 0; j < VNODES; j++){
            if (sum[j] >= 0) error_computed[j] = 0;
            else error_computed[j] = 1;
        }

        // ----------- DEBUG PRINT --------------
        //printf("\tError computed: ");
        for(int j = 0; j < VNODES; j++){
            //printf("%d ",error_computed[j]);
        }
        //printf("\n");
        //------------------------------------------

        // Compute S = eH^T
        int resulting_syndrome[CHECK];
        for(int i = 0; i < CHECK; i++){
            int row_op = 0;
            for(int j = 0; j < VNODES; j++){
                row_op ^= (error_computed[j] & pcm_matrix[i * VNODES + j]);
            }
            resulting_syndrome[i] = row_op;
        }
        int error_found = 1;

        // ----------- DEBUG PRINT --------------
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
        //else if (i == num_it - 1) color_printf(RED, "\nUSED ALL ITERATIONS WITHOUT FINDING THE ERROR");

        //printf("\n");
        //------------------------------------------
    }
}

void compute_row_operations(double *L,  int *non_zero, 
                            int* syndrome, int size_checks, int size_vnode)
{

    for(int i = 0; i < CHECK; i++){
        if(i == size_checks) break;

        double min1 = DBL_MAX, min2 = DBL_MAX;
        int minpos = -1;
        int sign_minpos = 0;
        int row_sign = 0;
        double product = 1.0;
        int nnz_row = 0;

        // Search min1 and min2
        for(int j = 0; j < VNODES; j++){
            if(j == size_vnode) break;

            double val = L[i * VNODES + j];
            double abs_val = fabs(val);

            if(non_zero[i * VNODES + j]){
                if(abs_val < min1){
                    min2 = min1;
                    min1 = abs_val;
                    minpos = j;
                    sign_minpos = (val >= 0 ? 0 : 1);
                }else if (abs_val < min2) {
                    min2 = abs_val;
                }
                nnz_row+=1;
            }

            row_sign = row_sign ^ (val >= 0 ? 0 : 1);
        }

        // Apply the corresponding value and sign to out[][]
        for(int j = 0; j < VNODES; j++){
            if(j == size_vnode) break;

            double val = L[i * VNODES + j];

            if(non_zero[i * VNODES + j]){
                // sign is negative (-1.0f) if the final signbit (operation in parethesis) is 0, 
                // and positive (1.0f) if its 1
                double sign =  1.0f - (2.0f * (row_sign ^ (val >= 0 ? 0 : 1) ^ syndrome[i]));

                // Assign min2 to minpos when loop ends to save if statements
                L[i * VNODES + j] = sign * min1;
            }
        }

        // Assigning min2 to minpos
        if(nnz_row > 1)
            L[i * VNODES + minpos] = (1.0f - 2.0f * (row_sign ^ sign_minpos ^ syndrome[i])) * min2;
    }
}

void compute_col_operations(double *L,  int *non_zero,
                            int* syndrome, int size_checks, int size_vnode, double alpha, 
                            double Lj[CHECK], double sum[VNODES])
{
    
    for (int j = 0; j < VNODES; j++){
        if (j == size_vnode) break;


        double sum_aux = 0.0f;
        int sum_count = 0;
        for(int i = 0; i < CHECK; i++){
            if (i == size_checks) break;

            sum_aux += L[i * VNODES + j];
            sum_count++;
        }
        if(sum_count > 0)
            sum[j] = Lj[j] + (alpha * sum_aux);
    }

    for (int i = 0; i < CHECK; i++){
        if(i == size_checks) break;

        for (int j = 0; j < VNODES; j++){
            if(j == size_vnode) break;

            if(non_zero[i * VNODES + j]) L[i * VNODES + j] = sum[j] - (alpha * L[i * VNODES + j]);
        }
    }

}

int main() {
    double *L;//[CHECK][VNODES];
    L = (double*)malloc(CHECK*VNODES*sizeof(double));
    int *pcm_matrix;//[CHECK][VNODES];
    pcm_matrix = (int*)malloc(CHECK*VNODES*sizeof(int));
    double Lj[VNODES];
    FILE *file = fopen("input_sliding.txt","r");
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
    p = 1.0;


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
                L[i * VNODES + j] = 0;
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
    int error_computed[VNODES];
    printf("Max Iterations: %d\n\n", num_it);

    printf("Initial L Matrix:\n");
    //show_matrix(L, pcm_matrix, rows, cols);

    min_sum(L, pcm_matrix, syndrome, rows, cols, Lj, alpha, num_it,&error_computed[0]);
    
    return 0;
}

void show_matrix( double *matrix, int *non_zero,
                  int rows,  int cols)
{
    for(int i = 0; i < rows; i++){
        printf("\t");
        for(int j = 0; j < cols; j++){
            if(signbit(matrix[i * VNODES + j])) {
                if (non_zero[i * VNODES + j]) color_printf(YELLOW," %.6f", matrix[i * VNODES + j]);
                else printf(" %.6f", matrix[i * VNODES + j]);
            }
            else {
                if (non_zero[i * VNODES + j]) color_printf(YELLOW,"  %.6f", matrix[i * VNODES + j]);
                else printf("  %.6f", matrix[i * VNODES + j]);
            }
        }
        printf("\n");
    }
    printf("\n");
}
