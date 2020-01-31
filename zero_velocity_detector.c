#include "zero_velocity_detector.h"

int *zero_velocity_detector(double u[][WINDOWS_SIZE])
{

    int i, k;
    int idx = 0;
    // Allocate memory
    int W = simdata.Window_size;
    int N = WINDOWS_SIZE;
    static int zupt[WINDOWS_SIZE]; //WINDOWS_SIZE here also contains buffer size. I'm thinking more about the embedding part of the code than the local version of it...
    memset(zupt, 0, WINDOWS_SIZE * sizeof(int));
    double T[N - W + 1];
    memset(T, 0, W * sizeof(double));

    // Run the desired detector type. Each detector return a vector with their
    // calculated test statistics T.
    switch (simdata.detector_type) {

        case GLRT:
            GLRT_f(u, T);
            break;

        case MV:
            //MV_f(u, T);
            break;

        case MAG:
            //MAG_f(u, T);
            break;

        case ARE:
            ARE_f(u, T);
            break;

        default:
            printf("The choosen detector type not recognized. The GLRT detector is used");
            GLRT_f(u, T);
            break;
    }

    // Check if the test statistics T are below the detector threshold. If so,
    // chose the hypothesis that the system has zero velocity
    for(k = 0; k < LENGTH(T); k++) {
        if(T[k] < simdata.gamma) {
            for(i = 0; i < W; i++)
                zupt[k + i] = 1;
        }
    }

    // Fix the edges of the detector statistics
    //T=[max(T)*ones(1,floor(W/2)) T max(T)*ones(1,floor(W/2))];
    //TODO
    int size = LENGTH(T) + floor(W / 2) * 2;
    double T_n[size];
    memset(T_n, 0, W * sizeof(double));

    for(i = 0; i < floor(W / 2); i++) {
        T_n[i] = max(T, LENGTH(T));
        idx++;
    }

    memcpy(T_n + idx, T, sizeof(T));
    idx += LENGTH(T);

    for(i = 0; i < floor(W / 2); i++)
        T_n[i + idx] = max(T, LENGTH(T));

    for(i = 0; i < LENGTH(T_n); i++) {
        //printf("%f ", T_n[i]);
    }

    for(i = 0; i < LENGTH(zupt); i++) {
        //printf("%f ", zupt[i]);
        //if(zupt[i] == 1)
        //printf("oui ");
        //else
        //printf("Non ");
    }

    return zupt;

}

void GLRT_f(double u[][WINDOWS_SIZE], double *T)
{

    double g = simdata.g;
    double sigma2_a = simdata.sigma_a * simdata.sigma_a;
    double sigma2_g = simdata.sigma_g * simdata.sigma_g;
    int W = simdata.Window_size;
    double ya_m[3] = {0};
    double tmp[3] = {0};
    double a = 0, b = 0;
    int i, k, l;

    int N = WINDOWS_SIZE;
    //double T[N-W+1] = {0};

    for(k = 0; k < N - W + 1; k++) {

        for(i = k; i < k + W - 2; i++) {

            ya_m[0] = u[0][i] + u[0][i + 1] + u[0][i + 2];
            ya_m[1] = u[1][i] + u[1][i + 1] + u[1][i + 2];
            ya_m[2] = u[2][i] + u[2][i + 1] + u[2][i + 2];

            //printf("%f %f %f\n", u[0][i], u[0][i+1], u[0][i+2]);
            //printf("%f %f %f\n", u[1][i], u[1][i+1], u[1][i+2]);
            //printf("%f %f %f\n", u[2][i], u[2][i+1], u[2][i+2]);
            //printf("\n");

        }

        ya_m[0] /= W;
        ya_m[1] /= W;
        ya_m[2] /= W;
        //printf("%f %f %f\n", ya_m[0], ya_m[1], ya_m[2]);




        for(l = k; l <= k + W - 1; l++) {

            tmp[0] = u[0][l] - (g * ya_m[0] / sqrt(ya_m[0] * ya_m[0] + ya_m[1] * ya_m[1] + ya_m[2] * ya_m[2]));
            tmp[1] = u[1][l] - (g * ya_m[1] / sqrt(ya_m[0] * ya_m[0] + ya_m[1] * ya_m[1] + ya_m[2] * ya_m[2]));
            tmp[2] = u[2][l] - (g * ya_m[2] / sqrt(ya_m[0] * ya_m[0] + ya_m[1] * ya_m[1] + ya_m[2] * ya_m[2]));
            //printf("%d \n", l);
            a = (u[3][l] * u[3][l] + u[4][l] * u[4][l] + u[5][l] * u[5][l]) / sigma2_g; //a = u(4:6,l)'*u(4:6,l)/sigma2_g
            //printf("%f \n", a);
            b = (tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]) / sigma2_a; //b = tmp'*tmp/sigma2_a
            T[k] = T[k] + a + b; // T(k)+u(4:6,l)'*u(4:6,l)/sigma2_g+tmp'*tmp/sigma2_a;


        }
    }

    for(i = 0; i < N - W + 1; i++)
        T[i] /= W;



}

void ARE_f(double u[][WINDOWS_SIZE], double *T)
{

    double g = simdata.g;
    double sigma2_g = simdata.sigma_g * simdata.sigma_g;
    int W = simdata.Window_size;
    int i = 0, k = 0, l = 0;

    int N = WINDOWS_SIZE;

    for(k = 0; k < N - W + 1; k++) {
        for(l = k; l <= k + W - 1; l++) {
            T[k] = T[k] + sqrt(u[3][l] * u[3][l] + u[4][l] * u[4][l] + u[5][l] * u[5][l]) * sqrt(u[3][l] * u[3][l] + u[4][l] * u[4][l] + u[5][l] * u[5][l]); //T(k)+norm(u(4:6,l))^2;
        }
    }

    for(i = 0; i < N - W + 1; i++)
        T[i] /= (sigma2_g * W);

}


int max(double arr[], int n)
{
    int i;
    int max = arr[0];

    for (i = 1; i < n; i++)  {
        if (arr[i] > max)
            max = arr[i];
    }

    return max;
}
