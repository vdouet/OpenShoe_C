#include "smoothed_ZUPTaidedINS.h"

void smoothed_ZUPTaidedINS(double u[][WINDOWS_SIZE], int zupt[], double x[][WINDOWS_SIZE], double cov_smoothed[][WINDOWS_SIZE], double P_smooth[][9][WINDOWS_SIZE], double dx[][WINDOWS_SIZE], double dx_smooth[][WINDOWS_SIZE], double P[][9][WINDOWS_SIZE])
{
    // Length of data
    int N = LENGTH_2D(u);
    int i, l, n;

    // Allocate state vectors
    static double quat[4][WINDOWS_SIZE];
    memset(quat, 0, 4 * WINDOWS_SIZE * sizeof (double));

    //static double dx[9][WINDOWS_SIZE];
    //memset(dx, 0, 9 * WINDOWS_SIZE * sizeof (double));

    static double dx_timeupd[9][WINDOWS_SIZE];
    memset(dx_timeupd, 0, 9 * WINDOWS_SIZE * sizeof (double));

    //static double dx_smooth[9][WINDOWS_SIZE];
    //memset(dx_smooth, 0, 9 * WINDOWS_SIZE * sizeof (double));

    // Allocate covariance matrices
    static double cov[9][WINDOWS_SIZE];
    memset(cov, 0, 9 * WINDOWS_SIZE * sizeof (double));

    //static double P[9][9][WINDOWS_SIZE];
    //memset (P, 0, 9 * 9 * WINDOWS_SIZE * sizeof (double));

    static double P_timeupd[9][9][WINDOWS_SIZE];
    memset (P_timeupd, 0, 9 * 9 * WINDOWS_SIZE * sizeof (double));

    //static double P_smooth[9][9][WINDOWS_SIZE];
    //memset (P_smooth, 0, 9 * 9 * WINDOWS_SIZE * sizeof (double));

    static double Fn[9][9][WINDOWS_SIZE];
    memset (Fn, 0, 9 * 9 * WINDOWS_SIZE * sizeof (double));

    double G[9][6] = {0};
    double F[9][9] = {0};

    // Constant matrices
    double Q[6][6] = {0};
    double R[3][3] = {0};
    int H[3][9] = {0};
    int Id[9][9] = {0};
    init_filter (P, Q, R, H);
    init_vec (cov, Id, P);

    // Initialize state vectors
    init_Nav_eq (u, x, quat);

    // Segment start/stop place holders
    int seg_start = 1;
    int seg_end = N;

    // Segment counter
    int c = 0;

    while(1) {

        /********************************/
        /**   Forward Kalman Filter   **/
        /******************************/

        for(n = seg_start; n <= seg_end; n++) {

            /*******************************/
            /**        Time Update       **/
            /*****************************/

            // Update mechanization equations
            double u_h[6] = {u[0][n], u[1][n], u[2][n], u[3][n], u[4][n], u[5][n]};
            Navigation_equations (x, u_h, quat, n);

            // Calculate state matrices
            state_matrix (quat, u_h, F, G, n);

            // Propagate errors
            for (i = 0; i < 9; i++) {
                for (l = 0; l < 9; l++) {
                    dx[i][n] += F[i][l] * dx[l][n - 1];
                }
            }

            // Update covariance matrices
            update_state_cov_smoothed(F, P, G, Q, cov, n); // I can optimize loop in that function for better performance

            // Save updated values for the smoothing
            for(i = 0; i < 9; i++) {
                dx_timeupd[i][n] = dx[i][n];

                for(l = 0; l < 9; l++) {
                    P_timeupd[i][l][n] = P[i][l][n];
                    Fn[i][l][n] = F[i][l]; //Did that to be able to use ZUPTaidedINS function where F is not declared similarly
                }
            }

            /*******************************/
            /**   Zero-velocity update   **/
            /*****************************/

            if (zupt[n] == true) {

                // Calculate the Kalman filter gain
                double K[9][3] = {0};
                kalman_filter_gain (K, P, H, R, n); // I can optimize loop in that function for better performance

                // Update deviation estimate
                update_deviation_estimate(dx, K, x, n);

                // Update the filter state covariance matrix P.
                update_state_cov_zupt_smoothed(P, Id, K, H, cov, n);
            }

            // Make sure the filter state covariance matrix is symmetric.
            double P_t[9][9] = {0};

            // P(:, :, k)'
            for (i = 0; i < 9; i++) {
                for (l = 0; l < 9; l++)
                    P_t[i][l] = P[l][i][n];
            }

            // P(:, :, k)=(P(:, :, k)+P(:, :, k)')/2;
            for (i = 0; i < 9; i++) {
                for (l = 0; l < 9; l++)
                    P[i][l][n] = (P[i][l][n] + P_t[i][l]) / 2;
            }

            // Store the diagonal of the state covariance matrix P.
            cov[0][n] = P[0][0][n];
            cov[1][n] = P[1][1][n];
            cov[2][n] = P[2][2][n];
            cov[3][n] = P[3][3][n];
            cov[4][n] = P[4][4][n];
            cov[5][n] = P[5][5][n];
            cov[6][n] = P[6][6][n];
            cov[7][n] = P[7][7][n];
            cov[8][n] = P[8][8][n];

            /********************************/
            /**   Segmentation decision   **/
            /******************************/

            if(c > 0) {
                c++;
            }

            double sum = cov[3][n - 1] + cov[4][n - 1] + cov[5][n - 1];
            double sum_n = cov[3][n] + cov[4][n] + cov[5][n];

            if(sum > 0.0001 && sum_n < 0.0001 && c == 0) {
                c = 1;
            }

            if(c == 30) { //What's that magic number ?
                seg_end = n;
                c = 0;
                break;
            }
        }

        /************************/
        /**   RTS smoothing   **/
        /**********************/

        // Initialize smoothing variables
        dx_smooth[0][seg_end] = dx[0][seg_end];
        dx_smooth[1][seg_end] = dx[1][seg_end];
        dx_smooth[2][seg_end] = dx[2][seg_end];
        dx_smooth[3][seg_end] = dx[3][seg_end];
        dx_smooth[4][seg_end] = dx[4][seg_end];
        dx_smooth[5][seg_end] = dx[5][seg_end];
        dx_smooth[6][seg_end] = dx[6][seg_end];
        dx_smooth[7][seg_end] = dx[7][seg_end];
        dx_smooth[8][seg_end] = dx[8][seg_end];

        for(i = 0; i < 9; i++) {
            for(l = 0; l < 9; l++) {
                P_smooth[i][l][seg_end] = P[i][l][seg_end];
            }
        }

        cov_smoothed[0][seg_end] = P_smooth[0][0][seg_end];
        cov_smoothed[1][seg_end] = P_smooth[1][1][seg_end];
        cov_smoothed[2][seg_end] = P_smooth[2][2][seg_end];
        cov_smoothed[3][seg_end] = P_smooth[3][3][seg_end];
        cov_smoothed[4][seg_end] = P_smooth[4][4][seg_end];
        cov_smoothed[5][seg_end] = P_smooth[5][5][seg_end];
        cov_smoothed[6][seg_end] = P_smooth[6][6][seg_end];
        cov_smoothed[7][seg_end] = P_smooth[7][7][seg_end];
        cov_smoothed[8][seg_end] = P_smooth[8][8][seg_end];

        for(n = seg_end - 1; n >= seg_start; n--) {

            // Variables needed for update equations
            double A[9][9] = {0};
            compute_A(A, P, Fn, P_timeupd, n);

            // Update states
            state_update(dx_smooth, dx, A, dx_timeupd, n);

            // Update covariance
            update_cov_smoothed(P_smooth, P, A, P_timeupd, cov_smoothed, n);

          }

        /*************************************/
        /**   Internal state compensation   **/
        /*************************************/

        for(n = seg_start; n <= seg_end; n++) {

            double dx_smooth_copy[9] = {
                -dx_smooth[0][n], -dx_smooth[1][n], -dx_smooth[2][n],
                    -dx_smooth[3][n], -dx_smooth[4][n], -dx_smooth[5][n],
                    -dx_smooth[6][n], -dx_smooth[7][n], -dx_smooth[8][n],
                };

            comp_internal_states(x, dx_smooth_copy, quat, n);

        }

        /************************/
        /**   Miscellaneous   **/
        /**********************/

        // Zero out last dx-value which will be used in next iteration
        dx[0][seg_end] = 0;
        dx[1][seg_end] = 0;
        dx[2][seg_end] = 0;
        dx[3][seg_end] = 0;
        dx[4][seg_end] = 0;
        dx[5][seg_end] = 0;
        dx[6][seg_end] = 0;
        dx[7][seg_end] = 0;
        dx[8][seg_end] = 0;

        // Zero out these values to avoid linearization problems
        P[0][8][seg_end] = 0;
        P[1][8][seg_end] = 0;
        P[8][0][seg_end] = 0;
        P[8][1][seg_end] = 0;

        // Save segmentation points
        // seg = [seg seg_end];

        // Check if end of data or correct segment place holders
        if(seg_end != N-1) {
            seg_start = seg_end + 1;
            seg_end = N-1;
        }

        else {
            // End of data -> break loop
            break;
        }

        // Go back and process new segment
    }
}

void update_cov_smoothed(double P_smooth[][9][WINDOWS_SIZE], double P[][9][WINDOWS_SIZE], double A[][9], double P_timeupd[][9][WINDOWS_SIZE], double cov_smoothed[][WINDOWS_SIZE], int n)
{

    //P_smooth(:, :, n) = P(:, :, n) + A*(P_smooth(:, :, n+1) - P_timeupd(:, :, n + 1))*A';
    int i, l, j;

    // A'
    double A_t[9][9] = {0};

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            A_t[i][l] = A[l][i];
    }

    //A*(P_smooth(:, :, n+1) - P_timeupd(:, :, n + 1))
    double A_mult_P[9][9] = {0};

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++) {
            for(j = 0; j < 9; j++) {
                A_mult_P[i][l] += A[i][j] * (P_smooth[j][l][n + 1] - P_timeupd[j][l][n + 1]);
            }
        }
    }

    // P(:, :, n) + A*(P_smooth(:, :, n+1) - P_timeupd(:, :, n + 1))*A'
    double part_A[9][9] = {0};

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++) {
            for(j = 0; j < 9; j++) {
                part_A[i][l] += A_mult_P[i][j] * A_t[j][l];
            }

            P_smooth[i][l][n] = P[i][l][n] + part_A[i][l];
        }
    }



    // Make sure the filter state covariance matrix is symmetric.
    double P_smooth_t[9][9] = {0};

    // P_smooth(:, :, n)'
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            P_smooth_t[i][l] = P_smooth[l][i][n];
    }

    // P_smooth(:, :, n) = (P_smooth(:, :, n) + P_smooth(:, :, n)')/2;
    for (i =
             0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            P_smooth[i][l][n] = (P_smooth[i][l][n] + P_smooth_t[i][l]) / 2;
    }

    // Store the diagonal of the state covariance matrix P.
    cov_smoothed[0][n] = P_smooth[0][0][n];
    cov_smoothed[1][n] = P_smooth[1][1][n];
    cov_smoothed[2][n] = P_smooth[2][2][n];
    cov_smoothed[3][n] = P_smooth[3][3][n];
    cov_smoothed[4][n] = P_smooth[4][4][n];
    cov_smoothed[5][n] = P_smooth[5][5][n];
    cov_smoothed[6][n] = P_smooth[6][6][n];
    cov_smoothed[7][n] = P_smooth[7][7][n];
    cov_smoothed[8][n] = P_smooth[8][8][n];

}

void state_update(double dx_smooth[9][WINDOWS_SIZE], double dx[9][WINDOWS_SIZE], double A[][9], double dx_timeupd[9][WINDOWS_SIZE], int n)
{

    //dx_smooth(:, n) = dx(:, n) + A*(dx_smooth(:, n + 1) - dx_timeupd(:, n + 1));
    int i, l;

    double part_A[9] = {0};

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++) {
            part_A[i] += A[i][l] * (dx_smooth[l][n + 1] - dx_timeupd[l][n + 1]);
        }

        dx_smooth[i][n] = dx[i][n] + part_A[i];
    }

}

void compute_A(double A[][9], double P[][9][WINDOWS_SIZE], double Fn[][9][WINDOWS_SIZE], double P_timeupd[][9][WINDOWS_SIZE], int k)
{
    //A = P(:, :, n)*F(:, :, n)'/P_timeupd(:, :, n + 1)
    int i, l, j;

    // F(:, :, n)'
    double F_t[9][9] = {0};

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            F_t[i][l] = Fn[l][i][k];
    }

    double part_A[9][9] = {0};

    // P(:, :, n)*F(:, :, n)'
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++) {
            for (j = 0; j < 9; j++)
                part_A[i][l] += P[i][j][k] * F_t[j][l];

        }
    }

    /* In the matlab implementation they used A/B which is matlab solving
     systems of linear equations xA = B for x.
     The matlab documentation tells us : "If A is a square matrix, then B/A is
     roughly equal to B*inv(A), but MATLAB processes B/A differently and more robustly."
     See : https://www.mathworks.com/help/matlab/ref/mrdivide.html */

    // A = P(:, :, n)*F(:, :, n)' * inv(P_timeupd(:, :, n + 1))
    double P_timeupd_inv[9][9] = {0};
    inv_mat_9(P_timeupd_inv, P_timeupd, k + 1);

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++) {
            for (j = 0; j < 9; j++)
                A[i][l] += part_A[i][j] * P_timeupd_inv[j][l];
        }
    }

}

/* Credit to Tapan Kumar Mishra : mishra.tapankumar [at] gmail.com for the Gauss Jordan algorithm */
// http://www.sourcecodesworld.com/source/show.asp?ScriptID=1086
void inv_mat_9(double I[][9], double P_timeupd[][9][WINDOWS_SIZE], int n)
{

    /* Original code from Tapan Kumar Mishra at http://www.sourcecodesworld.com/source/show.asp?ScriptID=1086 */

    double temp;
    int i, j, k, l, matsize = 9;

    // We make a copy of the matrix to invert
    double A[9][9] = {0};

    for(i = 0; i < 9; i++) {
      for(l = 0; l < 9; l++) {
        A[i][l] = P_timeupd[i][l][n];
      }
    }


    // Transform to identity matrix
    I[0][0] = 1;
    I[1][1] = 1;
    I[2][2] = 1;
    I[3][3] = 1;
    I[4][4] = 1;
    I[5][5] = 1;
    I[6][6] = 1;
    I[7][7] = 1;
    I[8][8] = 1;

    // We do not check if the matrix can be inverted or not

    /*---------------LoGiC starts here------------------*/		//procedure to make the matrix A to unit matrix

    for(k = 0; k < matsize; k++) {								//by some row operations,and the same row operations of
        //Unit mat. I gives the inverse of matrix A
        temp = A[k][k];										//'temp' stores the A[k][k] value so that A[k][k] will not change

        for(j = 0; j < matsize; j++) {							//during the operation A[i][j]/=A[k][k] when i=j=k
            A[k][j] /= temp;									//it performs the following row operations tomake A to unit matrix
            I[k][j] /= temp;									//R0=R0/A[0][0],similarly for I also R0=R0/A[0][0]
        }													              //R1=R1-R0*A[1][0] similarly for I

        for(i = 0; i < matsize; i++) {							//R2=R2-R0*A[2][0]		,,
            temp = A[i][k];									//R1=R1/A[1][1]

            for(j = 0; j < matsize; j++) {						//R0=R0-R1*A[0][1]
                //R2=R2-R1*A[2][1]
                if(i == k)
                    break;									//R2=R2/A[2][2]

                A[i][j] -= A[k][j] * temp;						//R0=R0-R2*A[0][2]
                I[i][j] -= I[k][j] * temp;						//R1=R1-R2*A[1][2]
            }
        }
    }

    /*---------------LoGiC ends here--------------------*/


}

void update_deviation_estimate(double dx[][WINDOWS_SIZE], double K[][3], double x[][WINDOWS_SIZE], int n)
{

    //dx(:, n)= dx(:, n) - K*(dx(4:6, n) - x(4:6,n));


    //(dx(4:6, n) - x(4:6,n))

    double res1[3] = {dx[3][n] - x[3][n], dx[4][n] - x[4][n], dx[5][n] - x[5][n]};

    // K*(dx(4:6, n) - x(4:6,n))

    double part_B[9] = {0};

    for (int i = 0; i < 9; i++) {
        for (int l = 0; l < 3; l++) {
            part_B[i] += K[i][l] * res1[l];
        }
    }

    //dx(:, n) - K*(dx(4:6, n) - x(4: 6, n));
    dx[0][n] -= part_B[0];
    dx[1][n] -= part_B[1];
    dx[2][n] -= part_B[2];
    dx[3][n] -= part_B[3];
    dx[4][n] -= part_B[4];
    dx[5][n] -= part_B[5];
    dx[6][n] -= part_B[6];
    dx[7][n] -= part_B[7];
    dx[8][n] -= part_B[8];

}

void update_state_cov_smoothed (double F[][9], double P[][9][WINDOWS_SIZE], double G[][6], double Q[][6], double cov[][WINDOWS_SIZE], int k)
{

    int i, l, j;

    // Transpose G and F.
    // They are doing a complex conjugate transpose but I don't think it matters here.
    double G_t[6][9] = {0};
    double F_t[9][9] = {0};

    for (i = 0; i < 6; i++) {
        for (l = 0; l < 9; l++)
            G_t[i][l] = G[l][i];
    }

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            F_t[i][l] = F[l][i];
    }

    double F_mult_P[9][9] = {0};
    double part_A[9][9] = {0};

    // F*P(:, :, k - 1)
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++) {
            for (j = 0; j < 9; j++)
                F_mult_P[i][l] += F[i][j] * P[j][l][k - 1];
        }
    }

    // Optimization here ? Merge both loops ?

    // F*P(:, :, k - 1)*F'
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++) {
            for (j = 0; j < 9; j++)
                part_A[i][l] += F_mult_P[i][j] * F_t[j][l];
        }
    }

    double G_mult_Q[9][6] = {0};
    double part_B[9][9] = {0};

    // G*Q
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 6; l++) {
            for (j = 0; j < 6; j++)
                G_mult_Q[i][l] += G[i][j] * Q[j][l];
        }
    }

    // Optimization here ? Merge both loops ?

    // G*Q*G'
    // I might have a problem here. Some signs are not the same. (value are 0's but still).
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++) {

            for (j = 0; j < 6; j++)
                part_B[i][l] += G_mult_Q[i][j] * G_t[j][l];
        }
    }

    // P(:, :, k) = F*P(:, :, k - 1)*F' + G*Q*G';
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            P[i][l][k] = part_A[i][l] + part_B[i][l];
    }

}

void update_state_cov_zupt_smoothed (double P[][9][WINDOWS_SIZE], int Id[][9], double K[][3], int H[][9], double cov[][WINDOWS_SIZE], int k)
{

    int i, l, j;

    // K*H
    double K_mult_H[9][9] = {0};

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++) {
            for (j = 0; j < 3; j++)
                K_mult_H[i][l] += K[i][j] * H[j][l];
        }
    }

    // (Id - K*H)
    double part_A[9][9] = {0};

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            part_A[i][l] = Id[i][l] - K_mult_H[i][l];
    }

    // P(:, :, k) = (Id - K*H)*P(:, :, k)
    double res[9][9] = {0};

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++) {
            for (j = 0; j < 9; j++)
                res[i][l] += part_A[i][j] * P[j][l][k];
        }
    }

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            P[i][l][k] = res[i][l];
    }

}
