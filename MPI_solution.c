#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

typedef enum Location {
    INSIDE,
    OUTSIDE,
    ON_EDGE
} Location;

//создаем матрицу
float** create_matrix(int n_rows, int n_cols) {
    float** matrix = malloc(n_rows * sizeof(float*));
    for (int i = 0; i < n_rows; i++) {
        matrix[i] = malloc(n_cols * sizeof(float));
    }
    return matrix;
}

//удаляем матрицу
void free_matrix(float** matrix, int n_rows) {
    for (int i = 0; i < n_rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void fill_matrix(float** matrix, int n_rows, int n_cols, float val) {
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            matrix[i][j] = val;
        }
    }
}

void print_matrix(float** matrix, int n_rows, int n_cols) {
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            printf("%.7f\t", matrix[i][j]);
        }
        printf("\n");
    }
}

// Функция для проверки, является ли число степенью двойки
int is_power_of_two(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

// Функция для нахождения степени числа 2
int get_power_of_two(int n) {
    int power = 0;
    while (n > 1) {
        n /= 2;
        power++;
    }
    return power;
}

Location get_location(float x, float y, float eps) {

    if (((fabs(x) < eps) && (y >= -eps)) || ((fabs(y) < eps) && (x >= -eps))){
        return ON_EDGE;
    }
    if ((x > 0) && (y > 0)) {
        return OUTSIDE;
    }
    return INSIDE;
}

float F(float x, float y, Location loc, float h1, float h2, float eps){
    switch (loc){
        case INSIDE:
            return 1;
        case OUTSIDE:
            return 0;
        case ON_EDGE:
            if ((fabs(x) < eps) && (fabs(y) < eps)){
                return 0.75; // S/(h1*h2), S = h1*h2/4
            }
            return 0.5;
    }
}

float ax(float x, float y, Location loc, float eps){
    switch (loc){
        case INSIDE:
            return 1;
        case OUTSIDE:
            return 1/eps;
        case ON_EDGE:
            if ((fabs(y) < eps) && (fabs(x) < eps)){
                return 1;
            }
            if (fabs(y) < eps){
                return 0.5 + (0.5)/eps;
            }
            return 1;
    }
}

//записываем матрицу
void writeMatrixToFile(const char* filename, int rows, int cols, float** matrix) {
    // Открываем файл для записи
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Ошибка открытия файла!\n");
        return;
    }

    // Проходим по строкам и столбцам матрицы и записываем её в файл
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(file, "%.8f ", matrix[i][j]);
        }
        fprintf(file, "\n");  // Перевод строки после каждой строки матрицы
    }

    // Закрываем файл
    fclose(file);
}


float compute_diff(int i, int j, float** w_matrix, float** ax_matrix, float** by_matrix, float h1, float h2){
    float x_diff = -(ax_matrix[i+1][j] * (w_matrix[i+1][j] - w_matrix[i][j]) - ax_matrix[i][j] * (w_matrix[i][j] - w_matrix[i-1][j]))/(h1*h1);
    float y_diff = -(by_matrix[i][j+1] * (w_matrix[i][j+1] - w_matrix[i][j]) - by_matrix[i][j] * (w_matrix[i][j] - w_matrix[i][j-1]))/(h2*h2);
    return x_diff + y_diff;
}

void send_neighbour_column(float** w_matrix, int column_number, int n_rows, int neighbor, int tag){
    
    float send_column[n_rows];
    for (int i = 0; i < n_rows; i++) {
        send_column[i] = w_matrix[i][column_number];  // Отправляем предпоследний столбец
    }
    MPI_Send(send_column, n_rows, MPI_FLOAT, neighbor, tag, MPI_COMM_WORLD);
    
}

void recv_neighbour_column(float** w_matrix, int column_number, int n_rows, int neighbor, MPI_Status status, int tag){
    
    float recv_column[n_rows];
    MPI_Recv(recv_column, n_rows, MPI_FLOAT, neighbor, tag, MPI_COMM_WORLD, &status);
    for (int i = 0; i < n_rows; i++) {
        w_matrix[i][column_number] = recv_column[i];  // Получаем данные в первый столбец
    }
    
}

void send_neighbour_row(float** w_matrix, int row_number, int n_cols, int neighbor, int tag){
    
    MPI_Send(w_matrix[row_number], n_cols, MPI_FLOAT, neighbor, tag, MPI_COMM_WORLD);
    
}

void recv_neighbour_row(float** w_matrix, int row_number, int n_cols, int neighbor, MPI_Status status, int tag){
    
    MPI_Recv(w_matrix[row_number], n_cols, MPI_FLOAT, neighbor, tag, MPI_COMM_WORLD, &status);
    
}

void exchange_data(float** w_matrix, int world_rank, int n_cols, int n_rows, MPI_Status status, 
    int right_neighbor, int left_neighbor, int bot_neighbor, int top_neighbor,
    int right_flag, int left_flag, int bot_flag, int top_flag, int rows_flag, int cols_flag){
    // Горизонтальный обмен (с левым и правым соседями)

    if (world_rank % 2 == 0) {  // Процессы с чётным рангом
        // Отправляем данные правому соседу, получаем данные слева
        if (!right_flag) {
            send_neighbour_column(w_matrix, n_cols - 2 - cols_flag, n_rows, right_neighbor,0);
            recv_neighbour_column(w_matrix, n_cols - 1, n_rows, right_neighbor, status,1);
        }
        if (!left_flag) {
            recv_neighbour_column(w_matrix, 0, n_rows, left_neighbor, status,0);
            send_neighbour_column(w_matrix, 1 + cols_flag, n_rows, left_neighbor,1);
        }
    } else {  // Процессы с нечётным рангом
        if (!left_flag) {
            recv_neighbour_column(w_matrix, 0, n_rows, left_neighbor, status,0);
            send_neighbour_column(w_matrix, 1 + cols_flag, n_rows, left_neighbor,1);
        }
        if (!right_flag) {
            send_neighbour_column(w_matrix, n_cols - 2 - cols_flag, n_rows, right_neighbor,0);
            recv_neighbour_column(w_matrix, n_cols - 1, n_rows, right_neighbor, status,1);
        }
    }


    if ((world_rank / n_cols)%2 == 0) {  // Процессы с чётными строками
        // Отправляем данные вниз, получаем данные сверху
        if (!bot_flag) {
            send_neighbour_row(w_matrix, n_cols - 2 - rows_flag, n_cols, bot_neighbor,2);
            recv_neighbour_row(w_matrix, n_cols - 1, n_cols, bot_neighbor, status,3);
        }
        if (!top_flag) {
            recv_neighbour_row(w_matrix, 0, n_cols, top_neighbor, status,2);
            send_neighbour_row(w_matrix, 1 + rows_flag, n_cols, top_neighbor,3);
        }
    } else {  // Процессы с нечётным рангом
        if (!top_flag) {
            recv_neighbour_row(w_matrix, 0, n_cols, top_neighbor, status,2);
            send_neighbour_row(w_matrix, 1 + rows_flag, n_cols, top_neighbor,3);
        }
        if (!bot_flag) {
            send_neighbour_row(w_matrix, n_cols - 2 - rows_flag, n_cols, bot_neighbor,2);
            recv_neighbour_row(w_matrix, n_cols - 1, n_cols, bot_neighbor, status,3);
        }
    }
}

int main (int argc, char *argv[])
{

    MPI_Init(&argc, &argv);

    double start_time = MPI_Wtime();

    if (argc < 2) {
        printf("Usage: %s <var1>\n", argv[0]);
        return 1;
    }

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // Ранг процесса

    if (!is_power_of_two(world_size)) {
        if (world_rank == 0) {
            printf("Ошибка: количество процессов (%d) должно быть степенью двойки.\n", world_size);
            MPI_Abort(MPI_COMM_WORLD, 1);  // Завершаем программу с ошибкой
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Получаем степень числа 2
    int power = get_power_of_two(world_size);
    int global_n_rows = power / 2;
    int global_n_cols = power - global_n_rows;
    global_n_rows += 1;
    global_n_cols += 1;

    int local_row = world_rank / global_n_cols;
    int local_col = world_rank % global_n_cols;

    if (world_rank == 0){
        printf("global shape: %d, %d\n", global_n_rows, global_n_cols);
    }

    int M, N;
    M = N = atoi(argv[1]);

    int n_rows, n_cols;
    n_rows = (global_n_rows - 1)?M/global_n_rows+2:M+1;
    n_cols = (global_n_cols - 1)?N/global_n_cols+2:N+1;

    //понимание границ для перекидывания данных
    int left_flag = (local_col == 0)?1:0;
    int right_flag = (local_col == global_n_cols-1)?1:0;
    int top_flag = (local_row == 0)?1:0;
    int bot_flag = (local_row == global_n_rows-1)?1:0;

    //fprintf(stderr, "proc:%i, lf:%i,rf:%i,tf:%i,bf:%i\n", world_rank, left_flag, right_flag, top_flag, bot_flag);

    MPI_Status status;
    
    int left_neighbor = world_rank - 1;
    int right_neighbor = world_rank + 1;

    int top_neighbor = world_rank - global_n_cols;
    int bot_neighbor = world_rank + global_n_cols;


    //начиная отсюда код касающийся вычислений

    float A1, B1, A2, B2;
    A1 = A2 = -1;
    B1 = B2 = 1;

    float h1 = (B1 - A1)/M;
    float h2 = (B2 - A2)/N;

    float h = (h1 > h2) ? h1 : h2; //max(h1, h2)
    float eps = h*h;

    float** global_w_matrix = create_matrix(M+1, N+1);
    float** global_Aw_matrix = create_matrix(M+1, N+1);
    float** global_F_matrix = create_matrix(M+1, N+1);
    float** global_ax_matrix = create_matrix(M+1, N+1);
    float** global_by_matrix = create_matrix(M+1, N+1);
    float** global_r_matrix = create_matrix(M+1, N+1);
    float** global_Ar_matrix = create_matrix(M+1, N+1);

    fill_matrix(global_w_matrix, M+1, N+1, 0);
    fill_matrix(global_Aw_matrix, M+1, N+1, 0);
    fill_matrix(global_F_matrix, M+1, N+1, 0);
    fill_matrix(global_ax_matrix, M+1, N+1, 0);
    fill_matrix(global_by_matrix, M+1, N+1, 0);
    fill_matrix(global_r_matrix, M+1, N+1, 0);
    fill_matrix(global_Ar_matrix, M+1, N+1, 0);

    float x, y;
    Location loc;
    for (int i=1; i<M+1; i++){
        x = A1 + i * h1;
        for (int j=1; j<N+1; j++){
            y = A2 + j * h2;
            //printf("x:%.2f,y:%.2f\n", x, y);
            loc = get_location(x, y, eps);
            global_F_matrix[i][j] = F(x,y,loc,h1,h2,eps);
            global_ax_matrix[i][j] = ax(x, y, loc, eps);
            global_by_matrix[i][j] = ax(y, x, loc, eps);
        }
    }

    //if (world_rank==0){
    //    printf("F matrix\n");
    //    print_matrix(global_F_matrix, M+1, N+1);
    //}
    
    MPI_Barrier(MPI_COMM_WORLD);

    //считаем размеры матриц, которые мы будем использовать

    float** w_matrix = create_matrix(n_rows, n_cols);
    float** Aw_matrix = create_matrix(n_rows, n_cols);
    float** F_matrix = create_matrix(n_rows, n_cols);
    float** ax_matrix = create_matrix(n_rows, n_cols);
    float** by_matrix = create_matrix(n_rows, n_cols);
    float** r_matrix = create_matrix(n_rows, n_cols);
    float** Ar_matrix = create_matrix(n_rows, n_cols);

    fill_matrix(w_matrix, n_rows, n_cols, 0);
    fill_matrix(Aw_matrix, n_rows, n_cols, 0);
    fill_matrix(F_matrix, n_rows, n_cols, 0);
    fill_matrix(ax_matrix, n_rows, n_cols, 0);
    fill_matrix(by_matrix, n_rows, n_cols, 0);
    fill_matrix(r_matrix, n_rows, n_cols, 0);
    fill_matrix(Ar_matrix, n_rows, n_cols, 0);

    if (world_rank==0){
        fprintf(stderr,"w_matrix shape:%i\n", n_rows);
    }

    int offset_row = ((M-1)/4)*2;//n_rows - (n_rows%2+2);
    int offset_col = ((N-1)/4)*2;//n_cols - (n_cols%2+2);
    if (M%4==0){
        offset_row +=1;
    }
    if (N%4==0){
        offset_col +=1;
    }

    if (world_rank==0){
        fprintf(stderr,"offset_row shape:%i\n", offset_row);
    }

    for (int i=0; i<n_rows; i++){
        for (int j=0; j<n_cols; j++){
            //w_matrix[i][j] = global_F_matrix[offset_row*local_row + i][offset_col*local_col + j];

            F_matrix[i][j] = global_F_matrix[offset_row*local_row + i][offset_col*local_col + j];
            ax_matrix[i][j] = global_ax_matrix[offset_row*local_row + i][offset_col*local_col + j];
            by_matrix[i][j] = global_by_matrix[offset_row*local_row + i][offset_col*local_col + j];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    float tau, rr, Arr, global_rr, global_Arr;
    global_rr = 1000;

    if (world_rank==0){
        fprintf(stderr, "exchange 1 started\n");
    }

    exchange_data(w_matrix, world_rank, n_cols, n_rows, status, 
                right_neighbor, left_neighbor, bot_neighbor, top_neighbor,
                right_flag, left_flag, bot_flag, top_flag, (M+1)%2, (N+1)%2);

    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank==0){
        fprintf(stderr, "exchange 1 complete\n");
    }

    float threshold = 0.005 > eps ? 0.005 : eps;
    int even_flag = (M+1) % 2;
    for (int k = 0; threshold < global_rr; k++) {
        rr = 0;
        Arr = 0;
        global_rr = 0;
        global_Arr = 0;

        {

            for (int i = 1; i < n_rows-1; i++) {
                float x = A1 + i * h1;
                for (int j = 1; j < n_cols-1; j++) {
                    float y = A2 + j * h2;
                    Aw_matrix[i][j] = compute_diff(i, j, w_matrix, ax_matrix, by_matrix, h1, h2);
                    r_matrix[i][j] = Aw_matrix[i][j] - F_matrix[i][j];
                }
            }

            exchange_data(r_matrix, world_rank, n_cols, n_rows, status, 
                right_neighbor, left_neighbor, bot_neighbor, top_neighbor,
                right_flag, left_flag, bot_flag, top_flag, (M+1)%2, (N+1)%2);

            for (int i = 1 + even_flag * 3 * local_row; i < n_rows-1; i++) {
                for (int j = 1 + even_flag * 3 * local_col; j < n_cols-1; j++) {
                    Ar_matrix[i][j] = compute_diff(i, j, r_matrix, ax_matrix, by_matrix, h1, h2);
                    
                    rr += r_matrix[i][j] * r_matrix[i][j] * h1*h2;
                    Arr += Ar_matrix[i][j] * r_matrix[i][j] * h1*h2;
                }
            }

            MPI_Reduce(&rr, &global_rr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_rr, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

            MPI_Reduce(&Arr, &global_Arr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Bcast(&global_Arr, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

            tau = global_rr / global_Arr;

            for (int i = 1; i < n_rows-1; i++) {
                for (int j = 1; j < n_cols-1; j++) {
                    w_matrix[i][j] = w_matrix[i][j] - tau * r_matrix[i][j];
                }
            }

            exchange_data(w_matrix, world_rank, n_cols, n_rows, status, 
                right_neighbor, left_neighbor, bot_neighbor, top_neighbor,
                right_flag, left_flag, bot_flag, top_flag, (M+1)%2, (N+1)%2);
        }

        // Условие остановки
        if ((global_rr < threshold) & (world_rank==0)) {
            printf("success! rr=%f, threshold=%f, k=%i\n", rr, threshold, k);
        }
    }

    //exchange_data(w_matrix, world_rank, n_cols, n_rows, status, 
    //    right_neighbor, left_neighbor, bot_neighbor, top_neighbor,
    //    right_flag, left_flag, bot_flag, top_flag, (M+1)%2, (N+1)%2);


    //printf("Выход из цикла\n");
    

    MPI_Barrier(MPI_COMM_WORLD);


    //объединение результатов
    if (world_rank == 0) {
        // Allocate global matrix
        float** global_matrix = create_matrix(global_n_rows * n_rows, global_n_cols * n_cols);

        printf("n_rows:%i, n_cols:%i\n", n_rows, n_cols);
        
        // Copy its own data into the global matrix
        for (int i = 0; i < n_rows; i++) {
            for (int j = 0; j < n_cols; j++) {
                int global_i = local_row * n_rows + i;
                int global_j = local_col * n_cols + j;
                global_matrix[global_i][global_j] = w_matrix[i][j];
            }
        }

        int displacement = 4 - ((M)%2)*2;

        // Receive data from other processes 
        for (int p = 1; p < world_size; p++) {
            int proc_local_row = p / global_n_cols;
            int proc_local_col = p % global_n_cols;

            printf("proc_local_row:%i,proc_local_col%i\n", proc_local_row, proc_local_col);

            float* recv_buf = malloc(n_rows * n_cols * sizeof(float));
            MPI_Recv(recv_buf, n_rows * n_cols, MPI_FLOAT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // Place data into the global matrix
            for (int i = 1; i < n_rows; i++) {
                for (int j = 1; j < n_cols; j++) {
                    int global_i = proc_local_row * (n_rows+even_flag) + i;
                    int global_j = proc_local_col * (n_cols+even_flag) + j;
                    global_matrix[global_i-displacement*proc_local_row][global_j-displacement*proc_local_col] = recv_buf[i * n_cols + j];
                }
            }
            free(recv_buf);

        }

        writeMatrixToFile("matrix.txt", global_n_rows * n_rows-displacement*(global_n_rows-1), global_n_cols * n_cols-displacement*(global_n_cols-1), global_matrix);
        
        free_matrix(global_matrix, global_n_rows * n_rows);
    } else {
        // Flatten the local matrix to send it
        float* w_matrix_flat = malloc(n_rows * n_cols * sizeof(float));
        for (int i = 0; i < n_rows; i++) {
            memcpy(&w_matrix_flat[i * n_cols], w_matrix[i], n_cols * sizeof(float));
        }
        
        // Send the local matrix to the root process
        MPI_Send(w_matrix_flat, n_rows * n_cols, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        free(w_matrix_flat);
    }

    free_matrix(w_matrix, n_rows);
    free_matrix(Aw_matrix, n_rows);
    free_matrix(F_matrix, n_rows);
    free_matrix(ax_matrix, n_rows);
    free_matrix(by_matrix, n_rows);
    free_matrix(r_matrix, n_rows);
    free_matrix(Ar_matrix, n_rows);

    free_matrix(global_w_matrix, n_rows);
    free_matrix(global_Aw_matrix, n_rows);
    free_matrix(global_F_matrix, n_rows);
    free_matrix(global_ax_matrix, n_rows);
    free_matrix(global_by_matrix, n_rows);
    free_matrix(global_r_matrix, n_rows);
    free_matrix(global_Ar_matrix, n_rows);

    // Конец замера времени
    double end_time = MPI_Wtime();

    // Вывод времени выполнения для каждого процесса
    printf("Process %d: Time taken = %f seconds\n", world_rank, end_time - start_time);


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

}