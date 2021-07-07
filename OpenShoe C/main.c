/*
Modified work Copyright (c) 2020 Victor Douet <victor.douet@gmail.com>, ISC License (open source)
Original work: Copyright (c) 2011 OpenShoe, ISC License (open source)

Change: Translated original Matlab implementation in C.
*/

#include "main.h"

/* To compile and link in command line do the following :
 *
 * gcc main.c settings.c zero_velocity_detector.c smoothed_ZUPTaidedINS.c ZUPTaidedINS.c -o OpenShoeC
 *
 */

int main(void)
{

    load_file_and_process();
    printf("Finished\n");
    return 0;
}

void load_file_and_process()
{

    /* OpenShoe algorithm variables */
    int *zupt;
    static double cov_closed[9][WINDOWS_SIZE] = {0};
    static double x_closed[9][WINDOWS_SIZE] = {0};
    static double u[6][WINDOWS_SIZE] = {0};
    static double x_smoothed[9][WINDOWS_SIZE] = {0};
    static double cov_smoothed[9][WINDOWS_SIZE] = {0};
    static double P_smooth[9][9][WINDOWS_SIZE] = {0};
    static double dx[9][WINDOWS_SIZE] = {0};
    static double dx_smooth[9][WINDOWS_SIZE] = {0};
    static double P[9][9][WINDOWS_SIZE] = {0};
    // memset(u, 0, 6 * WINDOWS_SIZE * sizeof(double));
    // memset(cov, 0, 9*WINDOWS_SIZE * sizeof(double)); //Re-initialize the static
    // array memset(x_h, 0, 9*WINDOWS_SIZE * sizeof(double)); //Re-initialize the
    // static array

    // Loads the algorithm settings and the IMU data
    printf("Loads the algorithm settings and the IMU data\n");
    settings(u);

    // Run the zero-velocity detector
    printf("Runs the zero velocity detector\n");
    zupt = zero_velocity_detector(u);

    // Run the costumary Kalman filter
    printf("Runs the customary Kalman filter\n");
    ZUPTaidedINS(u, zupt, x_closed, cov_closed);

    // Run the smoothing filter
    printf("Run the smoothing filter\n");
    smoothed_ZUPTaidedINS(u, zupt, x_smoothed, cov_smoothed, P_smooth, dx, dx_smooth, P);
    printf("done\n");

    // Write the ouput
    printf("Writing ouput\n");
    write_Output(x_closed);
    //write_Output(x_smoothed);
    //write_Output_zupt(zupt);

}
