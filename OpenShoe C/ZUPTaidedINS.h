#ifndef ZUPTAIDEDINS_H_   /* Include guard */
#define ZUPTAIDEDINS_H_

#include <stdbool.h>
#include "settings.h"


#define LENGTH_2D(x)  (sizeof(x)[0] / sizeof((x)[0][0]))
extern struct simdata_struct simdata;

void ZUPTaidedINS(double u[][WINDOWS_SIZE], int zupt[], double x_h[][WINDOWS_SIZE], double cov[][WINDOWS_SIZE]);
void init_filter(double P[][9][WINDOWS_SIZE], double Q[][6], double R[][3], int H[][9]);
void init_vec(double cov[][WINDOWS_SIZE], int Id[][9], double P[][9][WINDOWS_SIZE]);
void init_Nav_eq(double u[][WINDOWS_SIZE], double x_h[][WINDOWS_SIZE], double quat[][WINDOWS_SIZE]);
void Navigation_equations(double x[][WINDOWS_SIZE], double u[], double q[][WINDOWS_SIZE], int k);
void state_matrix(double q[][WINDOWS_SIZE], double u[], double F[][9], double G[][6], int k);
void update_state_cov(double F[][9], double P[][9][WINDOWS_SIZE], double G[][6], double Q[][6], double cov[][WINDOWS_SIZE], int k);
void kalman_filter_gain(double K[][3], double P[][9][WINDOWS_SIZE], int H[][9], double R[][3], int k);
void comp_internal_states(double x_h[][WINDOWS_SIZE], double dx[], double q[][WINDOWS_SIZE], int k);
void update_state_cov_zupt(double P[][9][WINDOWS_SIZE], int Id[][9], double K[][3], int H[][9], double cov[][WINDOWS_SIZE], int k);
void Rt2b(double ang[], double R[][3]);
void dcm2q(double R[][3], double quat[][WINDOWS_SIZE], int k);
void q2dcm(double R[][3], double q[][WINDOWS_SIZE], int k);
void inv_mat(double mat_inv[][3], double mat[][3]);


#endif // FOO_H_
