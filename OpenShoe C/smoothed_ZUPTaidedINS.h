#ifndef SMOOTHEDZUPTAIDEDINS_H_   /* Include guard */
#define SMOOTHEDZUPTAIDEDINS_H_

#include "settings.h"
#include "ZUPTaidedINS.h"

#define N_MAT 9

extern struct simdata_struct simdata;



void smoothed_ZUPTaidedINS(double u[][WINDOWS_SIZE], int zupt[], double x[][WINDOWS_SIZE], double cov_smoothed[][WINDOWS_SIZE], double P_smooth[][9][WINDOWS_SIZE], double dx[][WINDOWS_SIZE], double dx_smooth[][WINDOWS_SIZE], double P[][9][WINDOWS_SIZE]);
void update_cov_smoothed(double P_smooth[][9][WINDOWS_SIZE], double P[][9][WINDOWS_SIZE], double A[][9], double P_timeupd[][9][WINDOWS_SIZE], double cov_smoothed[][WINDOWS_SIZE] , int n);
void state_update(double dx_smooth[9][WINDOWS_SIZE], double dx[9][WINDOWS_SIZE], double A[][9], double dx_timeupd[9][WINDOWS_SIZE], int n);
void update_deviation_estimate(double dx[][WINDOWS_SIZE], double K[][3], double x[][WINDOWS_SIZE], int n);
void update_state_cov_zupt_smoothed (double P[][9][WINDOWS_SIZE], int Id[][9], double K[][3], int H[][9], double cov[][WINDOWS_SIZE], int k);
void update_state_cov_smoothed (double F[][9], double P[][9][WINDOWS_SIZE], double G[][6], double Q[][6], double cov[][WINDOWS_SIZE], int k);
void compute_A(double A[][9], double P[][9][WINDOWS_SIZE], double Fn[][9][WINDOWS_SIZE], double P_timeupd[][9][WINDOWS_SIZE], int k);
void inv_mat_9(double I[][9], double mat[][9][WINDOWS_SIZE], int k);




#endif // FOO_H_
