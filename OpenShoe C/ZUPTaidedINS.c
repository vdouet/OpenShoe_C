/*
Modified work Copyright (c) 2020 Victor Douet <victor.douet@gmail.com>, ISC License (open source)
Original work: Copyright (c) 2011 OpenShoe, ISC License (open source)

Change: Translated original Matlab implementation in C.
*/

#include "ZUPTaidedINS.h"

void ZUPTaidedINS (double u[][WINDOWS_SIZE], int zupt[], double x_h[][WINDOWS_SIZE], double cov[][WINDOWS_SIZE])
{

    // Get  the length of the IMU data vector
    // Should be equal to the windows_size (+ buff if running in local)
    int N = LENGTH_2D (u);
    //int N = WINDOWS_SIZE;

    // Initialize the filter state covariance matrix P, the processes noise
    // covariance matrix Q, the pseudo measurement noise covariance R, and the
    // observation matrix H.

    // Initial state covariance matrix
    static double P[9][9][WINDOWS_SIZE];
    memset (P, 0, 9 * 9 * WINDOWS_SIZE * sizeof (double)); //Re-initialize the static array

    // Process noise covariance matrix
    double Q[6][6] = {0};

    double R[3][3] = {0};

    // Observation matrix
    int H[3][9] = {0}; // 3 = lines ; 9 = columns

    init_filter (P, Q, R, H);

    // Allocate vectors
    int Id[9][9] = {0};
    init_vec (cov, Id, P);

    // Initialize the navigation state vector x_h, and the quaternion vector quat.
    //[x_h(1:9,1) quat(:, 1)]=init_Nav_eq(u);  % Subfunction located further down in the file.
    static double quat[4][WINDOWS_SIZE];
    memset (quat, 0, 4 * WINDOWS_SIZE * sizeof (double));
    init_Nav_eq (u, x_h, quat);

    /*******************************/
    /** Run the filter algorithm **/
    /*****************************/

    for (int k = 1; k < N; k++) {

        /*******************************/
        /**        Time Update       **/
        /*****************************/

        double u_h[6] = {u[0][k], u[1][k], u[2][k], u[3][k], u[4][k], u[5][k]};

        // Update the navigation equations.
        Navigation_equations (x_h, u_h, quat, k);

        // Update state transition matrix
        double F[9][9] = {0};
        double G[9][6] = {0};
        state_matrix (quat, u_h, F, G, k);

        // Update the filter state covariance matrix P
        update_state_cov (F, P, G, Q, cov, k); // I can optimize loop in that function for better performance

        /*******************************/
        /**   Zero-velocity update   **/
        /*****************************/

        // Check if a zero velocity update should be done. If so, do the
        // following
        if (zupt[k] == true) {

            // Calculate the Kalman filter gain
            double K[9][3] = {0};
            kalman_filter_gain (K, P, H, R, k); // I can optimize loop in that function for better performance

            // Calculate the prediction error. Since the detector hypothesis
            // is that the platform has zero velocity, the prediction error is
            // equal to zero minu the estimated velocity.
            double z[3] = {-x_h[3][k], -x_h[4][k], -x_h[5][k]};

            // Estimation of the perturbations in the estimated navigation
            // states
            double dx[9] = {0};

            for (int i = 0; i < 9; i++) {
                for (int l = 0; l < 3; l++)
                    dx[i] += K[i][l] * z[l];
            }

            // Correct the navigation state using the estimated perturbations.
            // (Subfunction located further down in the file.)
            comp_internal_states (x_h, dx, quat, k);


            // Update the filter state covariance matrix P.
            update_state_cov_zupt (P, Id, K, H, cov, k);

        }

    }

}

void init_filter (double P[][9][WINDOWS_SIZE], double Q[][6], double R[][3], int H[][9])
{

    // General values for the observation matrix H
    H[0][3] = 1;
    H[1][4] = 1;
    H[2][5] = 1;

    // General values for the initial covariance matrix P
    P[0][0][0] = simdata.sigma_initial_pos[0] * simdata.sigma_initial_pos[0];
    P[1][1][0] = simdata.sigma_initial_pos[1] * simdata.sigma_initial_pos[1];
    P[2][2][0] = simdata.sigma_initial_pos[2] * simdata.sigma_initial_pos[2];
    P[3][3][0] = simdata.sigma_initial_vel[0] * simdata.sigma_initial_vel[0];
    P[4][4][0] = simdata.sigma_initial_vel[1] * simdata.sigma_initial_vel[1];
    P[5][5][0] = simdata.sigma_initial_vel[2] * simdata.sigma_initial_vel[2];
    P[6][6][0] = simdata.sigma_initial_att[0] * simdata.sigma_initial_att[0];
    P[7][7][0] = simdata.sigma_initial_att[1] * simdata.sigma_initial_att[1];
    P[8][8][0] = simdata.sigma_initial_att[2] * simdata.sigma_initial_att[2];

    // General values for the process noise covariance matrix Q
    Q[0][0] = simdata.sigma_acc[0] * simdata.sigma_acc[0];
    Q[1][1] = simdata.sigma_acc[1] * simdata.sigma_acc[2];
    Q[2][2] = simdata.sigma_acc[2] * simdata.sigma_acc[2];
    Q[3][3] = simdata.sigma_gyro[0] * simdata.sigma_gyro[0];
    Q[4][4] = simdata.sigma_gyro[1] * simdata.sigma_gyro[1];
    Q[5][5] = simdata.sigma_gyro[2] * simdata.sigma_gyro[2];

    // General values for the measurement noise matrix R
    R[0][0] = simdata.sigma_vel[0] * simdata.sigma_vel[0];
    R[1][1] = simdata.sigma_vel[1] * simdata.sigma_vel[1];
    R[2][2] = simdata.sigma_vel[2] * simdata.sigma_vel[2];

}

void init_vec (double cov[][WINDOWS_SIZE], int Id[][9], double P[][9][WINDOWS_SIZE])
{

    for (int i = 0; i < 9; i++) {
        Id[i][i] = 1;
        cov[i][0] = P[i][i][0];
    }

}

void init_Nav_eq (double u[][WINDOWS_SIZE], double x_h[][WINDOWS_SIZE], double quat[][WINDOWS_SIZE])
{

    // Under the assumption that the system is stationary during the first 20
    // samples, the initial roll and pitch is calculate from the 20 first
    // accelerometer readings.
    /* This needs to be at the beginning of the program in the embedded system */
    double f_u = 0, f_v = 0, f_w = 0;

    for (int i = 0; i < 20; i++) {
        f_u += u[0][i];
        f_v += u[1][i];
        f_w += u[2][i];
    }

    f_u /= 20;
    f_v /= 20;
    f_w /= 20;

    double roll = atan2 (-f_v, -f_w);
    double pitch = atan2 (f_u, sqrt (f_v * f_v + f_w * f_w));

    // Set the attitude vector
    double attitude[3] = {roll, pitch, simdata.init_heading};

    // Calculate quaternion corresponing to the initial attitude
    double Rb2t[3][3] = {0};
    Rt2b (attitude, Rb2t); //The resulting Rb2t matrix is transposed
    dcm2q (Rb2t, quat, 0);

    x_h[0][0] = simdata.init_pos[0];
    x_h[1][0] = simdata.init_pos[1];
    x_h[2][0] = simdata.init_pos[2];
    x_h[6][0] = attitude[0];
    x_h[7][0] = attitude[1];
    x_h[8][0] = attitude[2];

}

void Navigation_equations (double x[][WINDOWS_SIZE], double u[], double q[][WINDOWS_SIZE], int k)
{

    double Ts = simdata.Ts;
    double identity_omega[4][4] = {0}; //new identity + new omega matrix
    double sum = 0;
    int i, l;

    /*************************************************************************/
    /* Update the quaternion vector "q"  given the angular rate measurements */
    /*************************************************************************/

    double w_tb[3] = {u[3], u[4], u[5]};

    double P = w_tb[0] * (Ts / 2);
    double Q = w_tb[1] * (Ts / 2);
    double R = w_tb[2] * (Ts / 2);

    double v = sqrt ((w_tb[0] * w_tb[0] + w_tb[1] * w_tb[1] + w_tb[2] * w_tb[2])) * Ts;

    if (v != 0) {

        double sin_part = (2 / v) * sin (v / 2);

        double new_omega[4][4] = { // angular rate matrix multiplied by Ts/2
            0, sin_part * R, sin_part *(-Q), sin_part * P,
            sin_part *(-R), 0, sin_part * P, sin_part * Q,
            sin_part * Q, sin_part *(-P), 0, sin_part * R,
            sin_part *(-P), sin_part *(-Q), sin_part *(-R), 0,
        };

        //We now need to multiply that cos part by an identity 4 matrix
        double new_I4[4][4] = {
            cos (v / 2), 0, 0, 0,
            0, cos (v / 2), 0, 0,
            0, 0, cos (v / 2), 0,
            0, 0, 0, cos (v / 2),
        };

        //We add the cos(.)*I4 (new I4) matrix and ksin(.)*Omega (new_Omega) matrix together
        for (i = 0; i < 4; i++) {
            for (l = 0; l < 4; l++)
                identity_omega[i][l] = (new_I4[i][l] + new_omega[i][l]);
        }

        //we multiply the result of the previous matrix by qn-1. [cos(.)I4 + ksin(.)Omega]*qn-1
        for (i = 0; i < 4; i++) {
            q[i][k] = 0;

            for (l = 0; l < 4; l++)
                q[i][k] += identity_omega[i][l] * q[l][k - 1];
        }

        // Calculate the norm of the quaternion matrix
        for (i = 0; i < 4; i++)
            sum += (q[i][k] * q[i][k]);

        double norm_q = sqrt ((double) sum);

        // Normalize quaternions
        for (i = 0; i < 4; i++)
            q[i][k] /= norm_q;

    }

    else {
        q[0][k] = q[0][k - 1];
        q[1][k] = q[1][k - 1];
        q[2][k] = q[2][k - 1];
        q[3][k] = q[3][k - 1];
    }

    /*************************************************************************/
    /* Use the update quaternion to get attitude of the navigation system in */
    /* terms of Euler angles.                                                */
    /*************************************************************************/

    // Get the roll, pitch and yaw
    double Rb2t[3][3] = {0};
    q2dcm (Rb2t, q, k);

    // roll
    x[6][k] = atan2 (Rb2t[2][1], Rb2t[2][2]);

    // pitch
    x[7][k] = atan (- Rb2t[2][0] / sqrt (Rb2t[2][1] * Rb2t[2][1] + Rb2t[2][2] * Rb2t[2][2]));

    //yaw
    x[8][k] = atan2 (Rb2t[1][0], Rb2t[0][0]);

    /*************************************************************************/
    /* Update position and velocity states using the measured specific force, */
    /* and the newly calculated attitude.                                     */
    /*************************************************************************/

    //Transform the specific force vector into navigation coordinate frame.
    double f_t[3] = {0};

    for (i = 0; i < 3; i++) {
        f_t[i] = 0;

        for (l = 0; l < 3; l++)
            f_t[i] += Rb2t[i][l] * u[l];
    }

    // Subtract (add) the gravity, to obtain accelerations in navigation
    // coordinat system.
    double acc_t[3] = {f_t[0], f_t[1], f_t[2] + simdata.g};

    // State space model matrices
    double A[6][6] = {
        1, 0, 0, Ts, 0, 0,
        0, 1, 0, 0, Ts, 0,
        0, 0, 1, 0, 0, Ts,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1,
    };

    double B[6][3] = {
        (Ts * Ts) / 2, 0, 0,
        0, (Ts * Ts) / 2, 0,
        0, 0, (Ts * Ts) / 2,
        Ts, 0, 0,
        0, Ts, 0,
        0, 0, Ts,
    };


    double res1[6] = {0};
    double res2[6] = {0};

    //A*x(1:6)
    for (i = 0; i < 6; i++) {
        res1[i] = 0;

        for (l = 0; l < 6; l++)
            res1[i] += A[i][l] * x[l][k - 1];
    }

    //B*acc_t
    for (i = 0; i < 6; i++) {
        res2[i] = 0;

        for (l = 0; l < 3; l++)
            res2[i] += B[i][l] * acc_t[l];
    }

    //y(1:6)=A*x(1:6)+B*acc_t;
    for (i = 0; i < 6; i++)
        x[i][k] = res1[i] + res2[i];

}

void state_matrix (double q[][WINDOWS_SIZE], double u[], double F[][9], double G[][6], int k)
{

    double Ts = simdata.Ts;
    // Convert quaternion to a rotation matrix
    double Rb2t[3][3] = {0};
    q2dcm (Rb2t, q, k);

    // Transform measured force to force in
    // the navigation coordinate system.
    double f_t[3] = {0};

    for (int i = 0; i < 3; i++) {
        f_t[i] = 0;

        for (int l = 0; l < 3; l++)
            f_t[i] += Rb2t[i][l] * u[l];
    }

    // Create a ske symmetric matrix of the specific fore vector
    // This is a vectorial product matrix
    double St[3][3] = {
        0, -f_t[2], f_t[1],
        f_t[2], 0, -f_t[0],
        -f_t[1], f_t[0], 0,
    };

    // Approximation of the discret time state transition matrices
    // (The transition matrix and noise gain matrix are directly included in
    // the declaration of Fc and Gc. More efficient this way.).
    double Fc[9][9] = {
        1, 0, 0, 1 * Ts, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 1 * Ts, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 1 * Ts, 0, 0, 0,
        0, 0, 0, 1, 0, 0, St[0][0] *Ts, St[0][1] *Ts, St[0][2] *Ts,
        0, 0, 0, 0, 1, 0, St[1][0] *Ts, St[1][1] *Ts, St[1][2] *Ts,
        0, 0, 0, 0, 0, 1, St[2][0] *Ts, St[2][1] *Ts, St[2][2] *Ts,
        0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1,
    };

    double Gc[9][6] = {
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        Rb2t[0][0] *Ts, Rb2t[0][1] *Ts, Rb2t[0][2] *Ts, 0, 0, 0,
        Rb2t[1][0] *Ts, Rb2t[1][1] *Ts, Rb2t[1][2] *Ts, 0, 0, 0,
        Rb2t[2][0] *Ts, Rb2t[2][1] *Ts, Rb2t[2][2] *Ts, 0, 0, 0,
        0, 0, 0, -Rb2t[0][0] *Ts, -Rb2t[0][1] *Ts, -Rb2t[0][2] *Ts,
        0, 0, 0, -Rb2t[1][0] *Ts, -Rb2t[1][1] *Ts, -Rb2t[1][2] *Ts,
        0, 0, 0, -Rb2t[2][0] *Ts, -Rb2t[2][1] *Ts, -Rb2t[2][2] *Ts,
    };

    memcpy (F, Fc, sizeof (Fc));
    memcpy (G, Gc, sizeof (Gc));

}

void update_state_cov (double F[][9], double P[][9][WINDOWS_SIZE], double G[][6], double Q[][6], double cov[][WINDOWS_SIZE], int k)
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

    // Make sure the filter state covariance matrix is symmetric.

    double P_t[9][9] = {0};

    // P(:, :, k)'
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            P_t[i][l] = P[l][i][k];
    }

    // P(:, :, k)=(P(:, :, k)+P(:, :, k)')/2;
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            P[i][l][k] = (P[i][l][k] + P_t[i][l]) / 2;
    }

    // Store the diagonal of the state covariance matrix P.
    cov[0][k] = P[0][0][k];
    cov[1][k] = P[1][1][k];
    cov[2][k] = P[2][2][k];
    cov[3][k] = P[3][3][k];
    cov[4][k] = P[4][4][k];
    cov[5][k] = P[5][5][k];
    cov[6][k] = P[6][6][k];
    cov[7][k] = P[7][7][k];
    cov[8][k] = P[8][8][k];

}

void kalman_filter_gain (double K[][3], double P[][9][WINDOWS_SIZE], int H[][9], double R[][3], int k)
{

    //K=(P(:, :, k)*H')/(H*P(:, :, k)*H'+R);
    int i, l, j;

    // Transpose G and F.
    // They are doing a complex conjugate transpose but I don't think it matters here.
    double H_t[9][3] = {0};

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 3; l++)
            H_t[i][l] = H[l][i];
    }

    double part_A[9][3] = {0};

    // P(:, :, k)*H'
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 3; l++) {

            for (j = 0; j < 9; j++)
                part_A[i][l] += P[i][j][k] * H_t[j][l];
        }
    }

    double H_mult_P[3][9] = {0};

    // H*P(:, :, k)
    for (i = 0; i < 3; i++) {
        for (l = 0; l < 9; l++) {
            for (j = 0; j < 9; j++)
         H_mult_P[i][l] += H[i][j] * P[j][l][k];
        }
    }

    double part_B[3][3] = {0};

    // H*P(:, :, k)*H'
    for (i = 0; i < 3; i++) {

        for (l = 0; l < 3; l++) {
            for (j = 0; j < 9; j++)
                part_B[i][l] += H_mult_P[i][j] * H_t[j][l];

            // H*P(:, :, k)*H'+R
            part_B[i][l] += R[i][l];
        }
    }

    /* In the matlab implementation they used A/B which is matlab solving
     systems of linear equations xA = B for x.
     The matlab documentation tells us : "If A is a square matrix, then B/A is
     roughly equal to B*inv(A), but MATLAB processes B/A differently and more robustly."
     See : https://www.mathworks.com/help/matlab/ref/mrdivide.html */

    double part_Binv[3][3] = {0};
    inv_mat (part_Binv, part_B);

    for (i = 0; i < 9; i++) {
        for (l = 0; l < 3; l++) {
            for (j = 0; j < 3; j++)
                K[i][l] += part_A[i][j] * part_Binv[j][l];
        }
    }

}

void comp_internal_states (double x_h[][WINDOWS_SIZE], double dx[], double q[][WINDOWS_SIZE], int k)
{

    // [x_h(:,k) quat(:, k)]=comp_internal_states(x_h(:,k),dx,quat(:, k));

    // Convert quaternion to a rotation matrix
    double R[3][3] = {0};
    q2dcm (R, q, k);

    // Correct the state vector
    x_h[0][k] += dx[0];
    x_h[1][k] += dx[1];
    x_h[2][k] += dx[2];
    x_h[3][k] += dx[3];
    x_h[4][k] += dx[4];
    x_h[5][k] += dx[5];
    x_h[6][k] += dx[6];
    x_h[7][k] += dx[7];
    x_h[8][k] += dx[8];

    // Correct the rotation matrics
    double epsilon[3] = {dx[6], dx[7], dx[8]};
    double OMEGA[3][3] = {
        1, dx[8], -dx[7],
        -dx[8], 1, dx[6],
        dx[7], -dx[6], 1,
    };


    // Just trying an optimized loop here by unwrapping two loops for matrix multiplication
    double R1[3][3] = {0};

    for (int j = 0; j < 3; j++) {

        R1[0][j] =  OMEGA[0][0] * R[0][j] +
                    OMEGA[0][1] * R[1][j] +
                    OMEGA[0][2] * R[2][j];

        R1[1][j] =  OMEGA[1][0] * R[0][j] +
                    OMEGA[1][1] * R[1][j] +
                    OMEGA[1][2] * R[2][j];

        R1[2][j] =  OMEGA[2][0] * R[0][j] +
                    OMEGA[2][1] * R[1][j] +
                    OMEGA[2][2] * R[2][j];
    }

    // Get the corrected roll, pitch and heading from the corrected rotation
    // matrix
    x_h[6][k] = atan2 (R1[2][1], R1[2][2]);
    x_h[7][k] = atan (-R1[2][0] / sqrt (R1[2][1] * R1[2][1] + R1[2][2] * R1[2][2]));
    x_h[8][k] = atan2 (R1[1][0], R1[0][0]);

    // Calculate the corrected quaternions
    dcm2q (R1, q, k);

}

void update_state_cov_zupt (double P[][9][WINDOWS_SIZE], int Id[][9], double K[][3], int H[][9], double cov[][WINDOWS_SIZE], int k)
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

    // Make sure the filter state covariance matrix is symmetric.
    double P_t[9][9] = {0};

    // P(:, :, k)'
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            P_t[i][l] = P[l][i][k];
    }

    // P(:, :, k)=(P(:, :, k)+P(:, :, k)')/2;
    for (i = 0; i < 9; i++) {
        for (l = 0; l < 9; l++)
            P[i][l][k] = (P[i][l][k] + P_t[i][l]) / 2;
    }

    // Store the diagonal of the state covariance matrix P.
    cov[0][k] = P[0][0][k];
    cov[1][k] = P[1][1][k];
    cov[2][k] = P[2][2][k];
    cov[3][k] = P[3][3][k];
    cov[4][k] = P[4][4][k];
    cov[5][k] = P[5][5][k];
    cov[6][k] = P[6][6][k];
    cov[7][k] = P[7][7][k];
    cov[8][k] = P[8][8][k];

}

void Rt2b (double ang[], double R[][3])  //It also transpose the resulting matrix
{

    double cr = cos (ang[0]);
    double sr = sin (ang[0]);

    double cp = cos (ang[1]);
    double sp = sin (ang[1]);

    double cy = cos (ang[2]);
    double sy = sin (ang[2]);

    //This also transpose the R matrix
    R[0][0] = cy * cp;
    R[1][0] = sy * cp;
    R[2][0] = -sp;

    R[0][1] = -sy * cr + cy * sp * sr;
    R[1][1] = cy * cr + sy * sp * sr;
    R[2][1] = cp * sr;

    R[0][2] = sy * sr + cy * sp * cr;
    R[1][2] = -cy * sr + sy * sp * cr;
    R[2][2] = cp * cr;

}

void dcm2q (double R[][3], double q[][WINDOWS_SIZE], int k)
{

    // qw = q[3], qx = q[0], qy = q[1], qz = q[2]

    //T = 1 + R(1,1) + R(2,2) + R(3,3);
    double T = 1 + R[0][0] + R[1][1] + R[2][2];

    if (T > 0.00000001) {

        // S = 0.5 / sqrt(T);
        // qw = 0.25 / S;
        // qx = ( R(3,2) - R(2,3) ) * S;
        // qy = ( R(1,3) - R(3,1) ) * S;
        // qz = ( R(2,1) - R(1,2) ) * S;

        double S = 0.5 / sqrt (T);
        q[3][k] = 0.25 / S;
        q[0][k] = ( R[2][1] - R[1][2] ) * S;
        q[1][k] = ( R[0][2] - R[2][0] ) * S;
        q[2][k] = ( R[1][0] - R[0][1] ) * S;

    }

    else if (R[0][0] > R[1][1] && R[0][0] > R[2][2]) {   //if (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))

        // S = sqrt( 1 + R(1,1) - R(2,2) - R(3,3)) * 2;
        // qw = (R(3,2) - R(2,3)) / S;
        // qx = 0.25 * S;
        // qy = (R(1,2) + R(2,1)) / S;
        // qz = (R(1,3) + R(3,1)) / S;

        double S = sqrt (1 + R[0][0] - R[1][1] - R[2][2]) * 2;
        q[3][k] = (R[2][1] - R[1][2]) / S;
        q[0][k] = 0.25 * S;
        q[1][k] = (R[0][1] + R[1][0]) / S;
        q[2][k] = (R[0][2] + R[2][0]) / S;

    }

    else if (R[1][1] > R[2][2]) {   // elseif (R(2,2) > R(3,3))

        // S = sqrt( 1 + R(2,2) - R(1,1) - R(3,3) ) * 2;
        // qw = (R(1,3) - R(3,1)) / S;
        // qx = (R(1,2) + R(2,1)) / S;
        // qy = 0.25 * S;
        // qz = (R(2,3) + R(3,2)) / S;

        double S = sqrt (1 + R[1][1] - R[0][0] - R[2][2]) * 2;
        q[3][k] = (R[0][2] - R[2][0]) / S;
        q[0][k] = (R[0][1] + R[1][0]) / S;
        q[1][k] = 0.25 * S;
        q[2][k] = (R[1][2] + R[2][1]) / S;

    }

    else {

        // S = sqrt( 1 + R(3,3) - R(1,1) - R(2,2) ) * 2;
        // qw = (R(2,1) - R(1,2)) / S;
        // qx = (R(1,3) + R(3,1)) / S;
        // qy = (R(2,3) + R(3,2)) / S;
        // qz = 0.25 * S;

        double S = sqrt (1 + R[2][2] - R[0][0] - R[1][1]) * 2;
        q[3][k] = (R[1][0] - R[0][1]) / S;
        q[0][k] = (R[0][2] + R[2][0]) / S;
        q[1][k] = (R[1][2] + R[2][1]) / S;
        q[2][k] = 0.25 * S;

    }

    // Do we need to normalize the quaternions ? They do that in their embedded firmware but not in matlab...
}

void q2dcm (double R[][3], double q[][WINDOWS_SIZE], int k)
{

    double p[6] = { q[0][k] *q[0][k], q[1][k] *q[1][k], q[2][k] *q[2][k], q[3][k] *q[3][k], 0, 0 };

    p[4] = p[1] + p[2];

    if (p[0] + p[3] + p[4] != 0)
        p[5] = 2 / (p[0] + p[3] + p[4]);

    else
        p[5] = 0;

    R[0][0] = 1 - p[5] * p[4];
    R[1][1] = 1 - p[5] * (p[0] + p[2]);
    R[2][2] = 1 - p[5] * (p[0] + p[1]);

    p[0] = p[5] * q[0][k];
    p[1] = p[5] * q[1][k];
    p[4] = p[5] * q[2][k] * q[3][k];
    p[5] = p[0] * q[1][k];

    R[0][1] = p[5] - p[4];
    R[1][0] = p[5] + p[4];

    p[4] = p[1] * q[3][k];
    p[5] = p[0] * q[2][k];

    R[0][2] = p[5] + p[4];
    R[2][0] = p[5] - p[4];

    p[4] = p[0] * q[3][k];
    p[5] = p[1] * q[2][k];

    R[1][2] = p[5] - p[4];
    R[2][1] = p[5] + p[4];

}

void inv_mat (double mat_inv[][3], double mat[][3])
{

    // Test with Gauss Jordan method ?

    //finding determinant
    double determinant = 0;

    for (int i = 0; i < 3; i++)
        determinant = determinant + (mat[0][i] * (mat[1][ (i + 1) % 3] * mat[2][ (i + 2) % 3] - mat[1][ (i + 2) % 3] * mat[2][ (i + 1) % 3]));

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            mat_inv[i][j] = ((mat[ (j + 1) % 3][ (i + 1) % 3] * mat[ (j + 2) % 3][ (i + 2) % 3]) - (mat[ (j + 1) % 3][(i + 2) % 3] * mat[ (j + 2) % 3][ (i + 1) % 3])) / determinant;
    }

}
