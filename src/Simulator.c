// Exemplo de simulador em C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void plotar_resultados(int N, double r[N], double resultados_P[5][N], double P_init, int dias_para_salvar[5]) {
    FILE *file = fopen("../reports/resultados.csv", "w");
    if (file == NULL) {
        printf("Erro ao criar o arquivo, veja se a pasta '../reports/' existe\n");
        return;
    }

    fprintf(file, "Raio(m)");
    for(int i = 0; i < 5; i++) {
        fprintf(file, ",Dia_%d(MPa)", dias_para_salvar[i]);
    }
    fprintf(file, "\n");

    for(int i = 0; i < N; i++) {
        fprintf(file, "%.4f", r[i]);
        for(int j = 0; j < 5; j++) {
            fprintf(file, ",%.4f", resultados_P[j][i] / 1e6);
        }
        fprintf(file, "\n");
    }

    fclose(file);
    printf("Resultados exportados para '../reports/resultados.csv'.\n");
}

void simular_injecao_co2() {
    double r_w = 0.1;
    double r_e = 2000.0;
    double h = 50.0;
    double phi = 0.25;
    double k = 100e-15;
    double mu = 0.06e-3;
    double ct = 5e-10;
    double q_inj = 0.02;
    double P_init = 15e6;

    FILE *fin = fopen("../data/parametros.txt", "r");
    if (fin != NULL) {
        fscanf(fin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &r_w, &r_e, &h, &phi, &k, &mu, &ct, &q_inj, &P_init);
        fclose(fin);
        printf("Parametros carregados de '../data/parametros.txt'\n");
    } else {
        printf("O arquivo '../data/parametros.txt' não foi encontrado, usando valores padrão.\n"); //Caso não ache valores no /data, usa valores padrão de exemplo
    }

    int N = 150;
    double dr = (r_e - r_w) / N;
    
    double r[N];
    for(int i = 0; i < N; i++) {
        r[i] = r_w + dr / 2.0 + i * dr;
    }
    
    double r_faces[N + 1];
    for(int i = 0; i <= N; i++) {
        r_faces[i] = r_w + i * dr;
    }

    int dias_simulacao = 365;
    double dt = 24.0 * 3600.0;
    int num_steps = (int)((dias_simulacao * 24.0 * 3600.0) / dt);

    double T[N + 1];
    T[0] = 0.0;
    T[N] = 0.0;
    for(int i = 1; i < N; i++) {
        T[i] = (2.0 * M_PI * k * h / mu) * r_faces[i] / dr;
    }

    double V[N];
    double C[N];
    for(int i = 0; i < N; i++) {
        double r_out = r[i] + dr / 2.0;
        double r_in = r[i] - dr / 2.0;
        V[i] = M_PI * (r_out * r_out - r_in * r_in) * h;
        C[i] = V[i] * phi * ct / dt;
    }

    double diag_inferior[N];
    double diag_principal[N];
    double diag_superior[N];

    for(int i = 0; i < N; i++) {
        diag_principal[i] = C[i] + T[i] + T[i+1];
        diag_inferior[i] = 0.0; 
        diag_superior[i] = 0.0;
    }
    for(int i = 1; i < N; i++) {
        diag_inferior[i] = -T[i];
    }
    for(int i = 0; i < N - 1; i++) {
        diag_superior[i] = -T[i+1];
    }

    double P[N];
    for(int i = 0; i < N; i++) {
        P[i] = P_init;
    }

    int dias_para_salvar[5] = {1, 30, 90, 180, 365};
    double resultados_P[5][N];

    double B[N];

    for(int step = 1; step <= num_steps; step++) {
        for(int i = 0; i < N; i++) {
            B[i] = C[i] * P[i];
        }
        B[0] += q_inj;
        
        spsolve(N, diag_inferior, diag_principal, diag_superior, B, P);
        
        for(int d = 0; d < 5; d++) {
            if(step == dias_para_salvar[d]) {
                for(int i = 0; i < N; i++) {
                    resultados_P[d][i] = P[i];
                }
            }
        }
    }

    plotar_resultados(N, r, resultados_P, P_init, dias_para_salvar);
}

int main() {
    simular_injecao_co2();
    return 0;
}
