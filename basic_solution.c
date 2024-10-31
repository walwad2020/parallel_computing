#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

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

float compute_diff(int i, int j, float** w_matrix, float** ax_matrix, float** by_matrix, float h1, float h2){
    float x_diff = -(ax_matrix[i+1][j] * (w_matrix[i+1][j] - w_matrix[i][j]) - ax_matrix[i][j] * (w_matrix[i][j] - w_matrix[i-1][j]))/(h1*h1);
    float y_diff = -(by_matrix[i][j+1] * (w_matrix[i][j+1] - w_matrix[i][j]) - by_matrix[i][j] * (w_matrix[i][j] - w_matrix[i][j-1]))/(h2*h2);
    return x_diff + y_diff;
}

int main (int argc, char *argv[])
{
    if (argc < 2) {
        printf("Usage: %s <var1>\n", argv[0]);
        return 1;
    }

    int M, N;
    M = N = atoi(argv[1]);

    clock_t begin = clock();
    
    int n_rows = M+1;
    int n_cols = N+1;

    float A1, B1, A2, B2;
    A1 = A2 = -1;
    B1 = B2 = 1;

    float h1 = (B1 - A1)/M;
    float h2 = (B2 - A2)/N;

    float h = (h1 > h2) ? h1 : h2; //max(h1, h2)
    float eps = h*h;

    printf("eps = %f\n",eps);


    float** w_matrix = create_matrix(n_rows, n_cols);
    float** Aw_matrix = create_matrix(n_rows, n_cols);
    float** F_matrix = create_matrix(n_rows, n_cols);
    float** ax_matrix = create_matrix(n_rows, n_cols);
    float** by_matrix = create_matrix(n_rows, n_cols);
    float** r_matrix = create_matrix(n_rows, n_cols);
    float** Ar_matrix = create_matrix(n_rows, n_cols);

    printf("xd2\n");

    fill_matrix(w_matrix, n_rows, n_cols, 0);
    fill_matrix(Aw_matrix, n_rows, n_cols, 0);
    fill_matrix(F_matrix, n_rows, n_cols, 0);
    fill_matrix(ax_matrix, n_rows, n_cols, 0);
    fill_matrix(by_matrix, n_rows, n_cols, 0);
    fill_matrix(r_matrix, n_rows, n_cols, 0);
    fill_matrix(Ar_matrix, n_rows, n_cols, 0);

    float x, y;
    Location loc;
    for (int i=1; i<M+1; i++){
        x = A1 + i * h1;
        for (int j=1; j<N+1; j++){
            y = A2 + j * h2;
            //printf("x:%.2f,y:%.2f\n", x, y);
            loc = get_location(x, y, eps);
            F_matrix[i][j] = F(x,y,loc,h1,h2,eps);
            ax_matrix[i][j] = ax(x, y, loc, eps);
            by_matrix[i][j] = ax(y, x, loc, eps);

        }
    }

    float ww_metric;

    float tau, rr, Arr;
    rr = 100;

    float timer = 10;

    float threshold = 0.005 > eps ? 0.005 : eps;

    for (int k=0;threshold < rr;k++){

    
        rr = Arr = 0;

        for (int i=1; i<M; i++){
            x = A1 + i * h1;
            for (int j=1; j<N; j++){
                y = A2 + j * h2;
                Aw_matrix[i][j] = compute_diff(i, j, w_matrix, ax_matrix, by_matrix, h1, h2);
                r_matrix[i][j] = Aw_matrix[i][j] - F_matrix[i][j];
            }
        }

        for (int i=1; i<M; i++){
            for (int j=1; j<N; j++){
                Ar_matrix[i][j] = compute_diff(i, j, r_matrix, ax_matrix, by_matrix, h1, h2);
                rr = rr + r_matrix[i][j] * r_matrix[i][j] * h1 * h2;
                Arr = Arr + Ar_matrix[i][j] * r_matrix[i][j] * h1 * h2;
            }
        }

        tau = rr / Arr;

        for (int i=1; i<M; i++){
            for (int j=1; j<N; j++){
                w_matrix[i][j] = w_matrix[i][j] - tau * r_matrix[i][j];
                
            }
        }
        
        if (rr < threshold){
            clock_t end = clock();
            double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
            printf("success! rr=%f, threshold=%f, time:%f, k=%i\n", rr, threshold, time_spent, k);
        }
    }




    printf("tau=%f\n",tau);
    writeMatrixToFile("matrix.txt", n_rows, n_cols, w_matrix);
    
    free_matrix(w_matrix, n_rows);
    free_matrix(Aw_matrix, n_rows);
    free_matrix(F_matrix, n_rows);
    free_matrix(ax_matrix, n_rows);
    free_matrix(by_matrix, n_rows);
    free_matrix(r_matrix, n_rows);
    free_matrix(Ar_matrix, n_rows);

}