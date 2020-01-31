#ifndef ZERO_VELOCITY_DETECTOR_H_   /* Include guard */
#define ZERO_VELOCITY_DETECTOR_H_

#include "settings.h"

#define LENGTH(x)  (sizeof(x) / sizeof((x)[0]))

extern struct simdata_struct simdata;
#define LEN simdata.Window_size

int *zero_velocity_detector(double u[][WINDOWS_SIZE]);
void GLRT_f(double u[][WINDOWS_SIZE], double *T);
void ARE_f(double u[][WINDOWS_SIZE], double *T);
int max(double arr[], int n);

#endif // FOO_H_
