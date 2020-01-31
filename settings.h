#ifndef SETTINGS_H_   /* Include guard */
#define SETTINGS_H_

#define pi    M_PI
#define GLRT  0
#define MV    1
#define MAG   2
#define ARE   3
#define WINDOWS_SIZE 3+51052 //Need to change the hardcoded file size

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void settings(double u[][WINDOWS_SIZE]);
double gravity(double lambda, double h);
void load_dataset(double u[][WINDOWS_SIZE]);
void write_Output(double x[9][WINDOWS_SIZE]);
void write_Output_zupt(int zupt[WINDOWS_SIZE]);

struct simdata_struct {
    double   altitude;
    double   latitude;
    double   init_heading;
    double   init_pos[3];
    double   sigma_acc[3];
    double   sigma_gyro[3];
    double   sigma_vel[3];
    double   sigma_a;
    double   sigma_g;
    int      Window_size;
    int      gamma;
    double   g;
    double   Ts;
    int      detector_type;
    double   sigma_initial_pos[3];
    double   sigma_initial_vel[3];
    double   sigma_initial_att[3];
};




#endif // FOO_H_
