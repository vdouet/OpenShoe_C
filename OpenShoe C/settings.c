/*
Modified work Copyright (c) 2020 Victor Douet <victor.douet@gmail.com>, ISC License (open source)
Original work: Copyright (c) 2011 OpenShoe, ISC License (open source)

Change: Translated original Matlab implementation in C.
*/

#include "settings.h"

/*************** settings.c ****************/

// Initialize simdata
struct simdata_struct simdata;

void settings(double u[][WINDOWS_SIZE])
{

    /*************** GENERAL PARAMETERS ****************/

    // Rough altitude [m]
    simdata.altitude = 100.0;

    // Roug latitude [degrees]
    simdata.latitude = 58.0;

    // Magnitude of the local gravity vector [m/s^2]
    simdata.g = gravity(simdata.latitude, simdata.altitude);

    // Sampling period [s]
    simdata.Ts = (double) 1 / 820;

    // Load the data
    load_dataset(u);

    // Initial heading [rad]
    simdata.init_heading = 0 * pi / 180;

    // Initial position (x,y,z)-axis [m]
    //memset(simdata.init_pos, 0, sizeof(simdata.init_pos));
    simdata.init_pos[0] = 0;
    simdata.init_pos[1] = 0;
    simdata.init_pos[2] = 0;

    /*************** Detector Settings ****************/

    // Detector type to be used. You can chose between:
    // GLRT - Generalized likelihood ratio test
    // MV -  Accelerometer measurements variance test
    // MAG - Accelerometer measurements magnitude test
    // ARE - Angular rate measurement energy test
    simdata.detector_type = GLRT;

    // Standard deviation of the acceleromter noise [m/s^2]. This is used to
    // control the zero-velocity detectors trust in the accelerometer data.
    simdata.sigma_a = 0.01;

    // Standard deviation of the gyroscope noise [rad/s]. This is used to
    // control the zero-velocity detectors trust in the gyroscope data.
    simdata.sigma_g = 0.1 * pi / 180;

    // Window size of the zero-velocity detector [samples]
    simdata.Window_size = 3;

    // Threshold used in the zero-velocity detector. If the test statistics are
    // below this value the zero-velocity hypothesis is chosen.
    simdata.gamma = 5000;

    /*************** FILTER PARAMETERS ****************/

    // Settings for the process noise, measurement noise, and initial state
    // covariance matrices Q, R, and P. All three matices are assumed to be
    // diagonal matrices, and all settings are defined as standard deviations.

    // Process noise for modeling the accelerometer noise (x,y,z platform
    // coordinate axis) and other accelerometer errors [m/s^2].
    simdata.sigma_acc[0] = 0.5 * 1;
    simdata.sigma_acc[1] = 0.5 * 1;
    simdata.sigma_acc[2] = 0.5 * 1;

    // Process noise for modeling the gyroscope noise (x,y,z platform coordinate
    // axis) and other gyroscope errors [rad/s].
    simdata.sigma_gyro[0] = 0.5 * 1 * pi / 180; //[rad/s]
    simdata.sigma_gyro[1] = 0.5 * 1 * pi / 180; //[rad/s]
    simdata.sigma_gyro[2] = 0.5 * 1 * pi / 180; //[rad/s]


    // Pseudo zero-velocity update measurement noise covariance (R). The
    // covariance matrix is assumed diagonal.
    simdata.sigma_vel[0] = 0.01; //[m/s]
    simdata.sigma_vel[1] = 0.01; //[m/s]
    simdata.sigma_vel[2] = 0.01; //[m/s]

    // Diagonal elements of the initial state covariance matrix (P).
    simdata.sigma_initial_pos[0] = 1e-5;  // Position (x,y,z navigation coordinate axis) [m]
    simdata.sigma_initial_pos[1] = 1e-5;  // Position (x,y,z navigation coordinate axis) [m]
    simdata.sigma_initial_pos[2] = 1e-5;  // Position (x,y,z navigation coordinate axis) [m]
    simdata.sigma_initial_vel[0] = 1e-5;  // Velocity (x,y,z navigation coordinate axis) [m/s]
    simdata.sigma_initial_vel[1] = 1e-5;  // Velocity (x,y,z navigation coordinate axis) [m/s]
    simdata.sigma_initial_vel[2] = 1e-5;  // Velocity (x,y,z navigation coordinate axis) [m/s]
    simdata.sigma_initial_att[0] = (pi / 180 * 0.1); // Attitude (roll,pitch,heading) [rad]
    simdata.sigma_initial_att[1] = (pi / 180 * 0.1); // Attitude (roll,pitch,heading) [rad]
    simdata.sigma_initial_att[2] = (pi / 180 * 0.1); // Attitude (roll,pitch,heading) [rad]

}

double gravity(double lambda, double h)
{

    lambda = pi / 180 * lambda;
    double gamma = 9.780327 * (1 + 0.0053024 * (sin(lambda) * sin(lambda)) - 0.0000058 * (sin(2 * lambda) * sin(2 * lambda)));
    double g = gamma - ((3.0877e-6) - (0.004e-6) * (sin(lambda) * sin(lambda))) * h + (0.072e-12) * (h * h);
    return g;
}

void load_dataset(double u[][WINDOWS_SIZE])
{

    //OpenShoe dataset
    int index_max = 5;
    int index_gx = 3;
    int index_gy = 4;
    int index_gz = 5;
    int index_ax = 0;
    int index_ay = 1;
    int index_az = 2;

    /* Opening and reading file variables */
    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *tok;
    int i = 0, k = 0;
    int nb_line = 0;
    double ax = 0, ay = 0, az = 0, gx = 0, gy = 0, gz = 0;

    fp = fopen("OpenShoeDataSet.txt", "r");

    if (fp == NULL) {
        printf("failed");
        getchar();
        exit(EXIT_FAILURE);
    }

    printf("\nReading file...\n");

    while ((read = getline(&line, &len, fp)) != -1) {
        tok = strtok(line, ",");

        while(i <= index_max) {

            if(i == index_ax)
                u[0][k] = atof(tok);

            if(i == index_ay)
                u[1][k] = atof(tok);

            if(i == index_az)
                u[2][k] = atof(tok);

            if(i == index_gx)
                u[3][k] = atof(tok);

            if(i == index_gy)
                u[4][k] = atof(tok);

            if(i == index_gz)
                u[5][k] = atof(tok);

            i++;
            tok = strtok(0, ",");
        }

        i = 0;
        k++;
        nb_line++;

    }

    printf("nb line read : %d \n", nb_line);

    fclose(fp);

    if (line)
        free(line);

}

void write_Output(double x[9][WINDOWS_SIZE]) {

  FILE *fp_output = fopen("result.txt", "w+");

  if (fp_output == NULL) {
      printf("failed");
      getchar();
      exit(EXIT_FAILURE);
  }

  for(int i = 0; i < WINDOWS_SIZE; i++) {
      fprintf(fp_output, "%f, %f, %f, %f, %f, %f, %f, %f, %f\n", x[0][i], x[1][i], x[2][i], x[3][i], x[4][i], x[5][i], x[6][i], x[7][i], x[8][i]);
  }

  fclose(fp_output);
}

void write_Output_zupt(int zupt[WINDOWS_SIZE]) {

  FILE *fp_output = fopen("result_zupt.txt", "w+");

  if (fp_output == NULL) {
      printf("failed");
      getchar();
      exit(EXIT_FAILURE);
  }

  for(int i = 0; i < WINDOWS_SIZE; i++) {
      fprintf(fp_output, "%d\n", zupt[i]);
  }

  fclose(fp_output);
}
