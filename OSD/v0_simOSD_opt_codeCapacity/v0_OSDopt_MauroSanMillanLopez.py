import numpy as np
import random
import time
from multiprocessing import Pool, cpu_count
from ldpc.bp_decoder import BpDecoder
from OSD_lib import OSD_decoder

np.random.seed(123)  # genera semilla de generador de números aleatorios

# FUNCIONES *******************************************************************************************************************************
def add_mod2(c1, c2):

    if len(c1) != len(c2):
        return -1

    for i in range(len(c1)):
        c1[i] = (c1[i] + c2[i]) % 2

    return c1


def to_mod2(mat):
    v_mod2 = []

    for a in mat:
        v_mod2.append(a % 2)

    return np.array(v_mod2)


def eliminacion_gaussiana(col1, col2, elem_estud):
    col = add_mod2(col1, col2)
    col[elem_estud] = 1  # Marcamos para saber que ahí se produjo una "colisión" con eliminación

    return col


def busca_posicion_columna(columna_colocar, H, numDetec):
    j = 0
    while j < numDetec and (columna_colocar[j] == 0 or H[j][j] == 1): # Mientras sea una columna válida y no encuentre su posición, por no tener capacidad o estar ya ocupada
        if columna_colocar[j] == 1 and H[j][j] == 1:  # Solo si está ya ocupado ese elemento de la diagonal
            columna_colocar = eliminacion_gaussiana(columna_colocar, H.T[j], j)

        j += 1

    return j


def eliminacion_salida(columna_eliminar, H, index_col, numDetec):   
    col = columna_eliminar

    for j in range(index_col+1, numDetec): # para toda columna hasta la salida

        if columna_eliminar[j] == 1 and H[j][j] == 1:  # Solo si está ya ocupado ese elemento de la diagonal
            col = eliminacion_gaussiana(columna_eliminar, H.T[j], j)

    return col


def fwrite(palabra, fichero):
    try: 
        with open(fichero): #player is the varible storing the username input
            with open(fichero, 'a') as f:
                f.write(str(palabra))

    except IOError:
        with open(fichero, 'w') as f:
            f.write(str(palabra))

# OSD *************************************************************************************************************************************
class OSD_module:

    def __init__(self, code_mtx):

        self.numDetec = code_mtx.shape[0]
        self.numQbits = code_mtx.shape[1]

        self.H = np.array((self.numDetec, self.numQbits))
        self.H = code_mtx


    def sort_probabilities(self, llr):  # O(n*log(n)), n =: número de bits del código. No merece la pena realizar el algoritmo descrito salvo si existe la posibilidad de paralelizar

        sort_columns = np.argsort(llr)  # Ordena índices de columna en función del valor asociado (menor a mayor)
        # ponderar llr con las probabilidades iniciales u otra información anterior

        return sort_columns


    def make_H(self, sort_columns):

        H = self.H.copy() # No queremos modificar la original

        H_lu = np.full((self.numDetec, self.numDetec), 0) # Creo la matriz cuadrada inicial

        used_cols = [-1] * self.numDetec

        for i in sort_columns:  # O(m^2), se puede mejorar si podemos hacer las búsquedas "segmentadas" -> O(m); m =: nº columnas matriz (qubits)
            # Busca su posición en la matriz Hs
            pos_insert = busca_posicion_columna(H.T[i], H_lu, self.numDetec)

            if  pos_insert < self.numDetec: # En caso de que sea una posición válida, se coloca la columna

                used_cols[pos_insert] = i
                H_lu.T[pos_insert] = H.T[i] # introduzco la columna

        H_s = np.full((self.numDetec, self.numDetec), 0) # esta será la matriz cuadrada inversa de la seleccionada

        i = 0
        for col in H_lu.T:  # O(n); m =: nº de detectores
            # introduzco columna definitiva al reducirse con todas las demás
            H_s.T[i] = eliminacion_salida(col, H_lu, i, self.numDetec)

            i += 1

        # No necesario porque realizamos permutación de columnas y se corresponde con la solución al sistema
        #H_s = H_s[used_cols] # aplicamos permutación realizada por la ordenación (por filas)

        return (H_s, used_cols)


    def solve_ecuationSys(self, synd, H, used_cols):

        sis_sol = to_mod2(H @ synd)

        est_error = [0 for i in range(self.numQbits)]

        for i in range(len(sis_sol)): # aplico los errores sacados del resultado del sistema resuelto

            if(sis_sol[i] == 1):
                est_error[used_cols[i]] = 1

        return np.array(est_error)


    def decode(self, synd, llr):

        sort_columns = self.sort_probabilities(llr)
        
        H_min = self.make_H(sort_columns)

        error_corr = self.solve_ecuationSys(synd, H_min[0], H_min[1])

        return error_corr

# LECTURA DE DATOS ************************************************************************************************************************
# ------- MATRIZ DEL CÓDIGO (control de paridad) -------
filenameIn = "../h_rom_syn.txt"
H = []

with open(filenameIn, 'r') as f:

    nRows = int(f.readline())
    nCols = int(f.readline())
    for j in range(nCols):
        col = f.readline()
        H.append([])
        for i in range(1, nRows + 1):
            H[j].append(int(col[i]))

    H = np.array(H, dtype=int)
    H = H.T  # se traspone porque la entrada es por cada fila una columna

# ------- MATRIZ LÓGICA --------------------------------
filename_log = "../B1_logic.txt"
H_log = []

with open(filename_log, 'r') as f:

    nRows = int(f.readline())
    nCols = int(f.readline())
    for j in range(nCols):
        col = f.readline().split()
        H_log.append(col)

    H_log = np.array(H_log, dtype=int)
    H_log = H_log.T  # se traspone porque la entrada es por cada fila una columna

# MULTIPROCESSING *************************************************************************************************************************

NUM_WORKERS = 20  # IMPORTANTE!! Este número debe dividir al número de episodios, ajustando al número de núcleos que se tengan
VEC_SIZE    = 882 # IMPORTANTE!! Debe coincidir con el número de columnas en la matriz H del código simulado


def generate_job(n_epochs, error_fis):

    error_fis_x = error_fis * (2/3) # (2/3 porque solo tengo que meter ruido en X)

    job_worker = []
    for _ in range(n_epochs//NUM_WORKERS):  # distribuimos la carga de trabajo de manera homogénea
        rand_error = [random.random() for _ in range(VEC_SIZE)]

        error_vector = []
        for p in rand_error:
            if(error_fis_x > p): 
                error_vector.append(1)
            else:
                error_vector.append(0)

        job_worker.append(np.array(error_vector))

    return job_worker


def worker(worker_id, n_epochs, error_fis):

    data = generate_job(n_epochs, error_fis)
    
    # BP counters
    global_nBP          = 0
    global_nBP_success  = 0
    global_nBP_fail     = 0
    # OSD counters
    global_nOSD         = 0
    global_nOSD_success = 0
    global_nOSD_fail    = 0
    # Decoder counter (global, "logic")
    global_logic_success = 0
    global_logic_fail    = 0

    error_fis_x = [ error_fis * (2/3) for _ in range(VEC_SIZE)] # (2/3 porque solo tengo que meter ruido en X)
    error_channel = np.ones(nCols) * error_fis_x # Inicialmente introducimos el nivel de ruido para los bits

    bp = BpDecoder(pcm=H, error_channel=error_channel, max_iter=100, bp_method='minimum_sum', ms_scaling_factor=0.625, schedule="parallel") # por defecto usa "minimun_sum"
    osd = OSD_decoder.OSD_decoder_opt(H.T) # módulo OSD no varía en función del error físico

    # BP local marks 
    nBP_marks = 0
    # OSD local marks
    nOSD_marks       = 0
    nOSD_fail_marks  = 0
    # Logical fail mark
    logic_fail_marks = 0
    # Time group
    t_total = 0
    t_osd   = 0

    ep = 0
    for error_vector in data:

        t_start = time.time()

        synd = to_mod2(H @ error_vector.T)

        sol = np.array([0 for _ in range(len(error_vector))])

        if synd.any():  # Se activa cuando no todos los elementos son 0. (Realiza Or(synd))
            # MÓDULO BP
            sol = bp.decode(synd)
            
            if bp.converge == False: # BP NO CONVERGE, debe entrar OSD para intentar corregir
                
                global_nBP_fail += 1

                # MÓDULO OSD
                t_osd_start = time.time()

                sol = osd.decode(synd, bp.log_prob_ratios)

                t_osd_end = time.time()
                t_osd += t_osd_end - t_osd_start

                # Comprobación de corrección
                if(not(np.equal(sol, error_vector).all())):
                    global_nOSD_fail += 1
                    nOSD_fail_marks  += 1
                else:
                    global_nOSD_success += 1

                global_nOSD += 1
                nOSD_marks  += 1
            else:                   # BP CONVERGE, comprobamos que haya corregido correctamente
                global_nBP_success += 1

            global_nBP += 1
            nBP_marks  += 1

        # ERROR LÓGICO
        error_logico = to_mod2(H_log @ (to_mod2(sol + error_vector)).T)

        if error_logico.any(): # S_l != 0
            global_logic_fail += 1
            logic_fail_marks  += 1
        else:
            global_logic_success += 1

        t_end = time.time()
        t_exe = t_end - t_start
        t_total += t_exe

        ep += 1

        if(ep % (len(data)//5) == 0): # relación 5:1 para agrupar simulaciones y mandarlas al final (50 -> 10; 500 -> 100; ...)
            fwrite("\nERROR " + str(error_fis) + ", EPISODE " + str(ep), "worker" + str(worker_id) + ".txt")
            fwrite("\tE_fis = "  + str(error_fis) + " ; E_log = " + str(global_logic_fail/ep) + " ( " + str(global_logic_fail) + " times fail )\n", "worker" + str(worker_id) + ".txt")
            fwrite("\tuses BP "  + str(global_nBP) + " times\n", "worker" + str(worker_id) + ".txt")
            fwrite("\tuses OSD " + str(global_nOSD) + " times\n", "worker" + str(worker_id) + ".txt")
            fwrite("\tOSD fails "+ str(global_nOSD_fail) + " times\n", "worker" + str(worker_id) + ".txt")
            fwrite("\ttime avrg " + str(t_total / (len(data)/5)) + " s\n", "worker" + str(worker_id) + ".txt")
            fwrite("\tosd time  " + str(t_osd / nOSD_marks) + " s\n", "worker" + str(worker_id) + ".txt")

    fwrite("_______________________________________________________________________________\n" , "worker" + str(worker_id) + ".txt")


    return [ep, global_logic_fail, global_nBP, global_nOSD, global_nOSD_fail, t_total, t_osd]

# MAIN ************************************************************************

def main():

    error_fisico_values = [0.1, 0.09, 0.08, 0.07, 0.06, 0.05]
    error_logico_values = []

    # para imitar la gráfica de https://arxiv.org/pdf/2205.06125 F-NMS-OSD
    epochs_values = [3, 4, 4, 5, 6, 7]

    start_exe = time.time()

    for i in range(len(error_fisico_values)):

        args = [(j, 10**epochs_values[i], error_fisico_values[i]) for j in range(NUM_WORKERS)]

        with Pool(NUM_WORKERS) as pool:
            results = pool.starmap(worker, args)

        # Recogemos resultados
        total_epochs       = 0
        total_error_logico = 0
        total_BP           = 0
        total_OSD          = 0
        total_OSD_fail     = 0
        total_time_exe     = 0
        total_time_osd     = 0
        total_epochs       = 0

        for res in results:
        
            epochs       = res[0] 
            error_logico = res[1]
            nBP          = res[2]
            nOSD         = res[3]
            nOSD_fails   = res[4]
            time_exe     = res[5]
            time_osd     = res[6]

            total_epochs       += epochs
            total_error_logico += error_logico
            total_BP           += nBP
            total_OSD          += nOSD
            total_OSD_fail     += nOSD_fails
            total_time_exe     += time_exe
            total_time_osd     += time_osd

        error_logico_values.append(total_error_logico/total_epochs)

        # Resultados finales para un valor de error físico
        print("\nE_fis = "  + str(error_fisico_values[i]) + " ; E_log = " + str(total_error_logico/total_epochs) + " ( " + str(total_error_logico) + " times fail )")
        print("\tnum epochs = " + str(total_epochs))
        print("\tuses BP "  + str(total_BP) + " times")
        print("\tuses OSD "  + str(total_OSD) + " times")
        print("\tOSD fails "+ str(total_OSD_fail) + " times")
        print("\taverage execution time " + str(total_time_exe/total_epochs) + " s")
        print("\taverage osd time "       + str(total_time_osd/total_OSD) + " s")

    end_exe = time.time()

    total_time_exe = end_exe - start_exe

    fwrite("RESUMEN:\n", "out.txt")
    fwrite("\nError fisico: " + str(error_fisico_values) + '\n', "out.txt")
    fwrite("Error logico: "   + str(error_logico_values) + '\n', "out.txt")
    fwrite("\nTotal execution time: " + str(total_time_exe) + " s\n", "out.txt")
    fwrite("*******\n", "out.txt")


if __name__ == "__main__":
    main()
