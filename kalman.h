//
// Created by bllgg on 1/19/21.
//

#ifndef MADGWICK_TEST_KALMAN_H
#define MADGWICK_TEST_KALMAN_H

class Kalman{
public:
    /*==============================================================*/
    /*                       Public Variables                       */
    /*==============================================================*/
    // no public variables

    /*==============================================================*/
    /*                  Constructor and Destructor                  */
    /*==============================================================*/
    Kalman();
    ~Kalman();

    /*==============================================================*/
    /*                       Public Functions                       */
    /*==============================================================*/
    //ret val   function name               parameters
    void        KalmanFilter                (double gx, double gy, double gz, double ax, double ay, double az);

    /*==============================================================*/
    /*                     Getters and Setters                      */
    /*==============================================================*/
    double      get_roll                    (void);
    void        set_roll                    (double roll);

    double      get_pitch                   (void);
    void        set_pitch                   (double pitch);

private:

    /*==============================================================*/
    /*                      Private Variables                       */
    /*==============================================================*/
    double      roll = 0;
    double      pitch = 0;
    double      gyro_noise = 0;
    double      accel_noise = 0;

    double      theta_hat_roll = 0;
    double      theta_hat_pitch = 0;
    double      theta_t_minus_1_roll = 0;
    double      theta_t_minus_1_pitch = 0;
    double      delta_t = 0.01;
    double      p_hat_roll = 1;
    double      p_hat_pitch = 1;
    double      z_k_roll = 0;
    double      z_k_pitch = 0;
    double      K_roll = 0;
    double      K_pitch = 0;

    /*==============================================================*/
    /*                      Private Functions                       */
    /*==============================================================*/
    //ret val   function name               parameters
    double      radians_to_degrees          (double radian_val);
    double      degrees_to_radians          (double degree_val);
};

#endif //MADGWICK_TEST_KALMAN_H
