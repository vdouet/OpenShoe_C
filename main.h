#ifndef MAIN_H_   /* Include guard */
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "settings.h"
#include "zero_velocity_detector.h"
#include "ZUPTaidedINS.h"
#include "smoothed_ZUPTaidedINS.h"



void navigation_equations(double ax, double ay, double az, double gx, double gy, double gz);
void step_estimation(double ax, double ay, double az, double gx, double gy, double gz);
void load_file_and_process();
void read_file(double u[][WINDOWS_SIZE]);

#endif // FOO_H_
