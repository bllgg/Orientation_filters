//
// Created by bllgg on 1/10/21.
// Directly Inherited from Madgwick filter
//

#ifndef MADGWICK_TEST_AHRS_H
#define MADGWICK_TEST_AHRS_H

class Madgwick {

public:

    /*==============================================================*/
    /*                       Public Variables                       */
    /*==============================================================*/
    // no public variables

    /*==============================================================*/
    /*                  Constructor and Destructor                  */
    /*==============================================================*/
    Madgwick(double beta_val);
    ~Madgwick();

    /*==============================================================*/
    /*                       Public Functions                       */
    /*==============================================================*/
    //ret val   function name               parameters
    void        MadgwickAHRSupdate          (double gx, double gy, double gz, double ax, double ay, double az, double mx, double my, double mz);
    void        MadgwickAHRSupdateIMU       (double gx, double gy, double gz, double ax, double ay, double az);
//    void        compute_orientation         (double *q);
//    double*      quaternionMul               (double *q_1, double *q_2);
//    double*      get_accel_jacobian          (double *q);
//    double*      get_accel_function          (double *q, double *a);
//    void        normalizeq                  (double *q);
//    void        update_roll_and_pitch       ();
//    void        get_mag_jacob               ();
//    void        get_mag_func                ();
//    void        get_rotation_mat            ();
//    void        update_roll_pitch_yaw       ();

    /*==============================================================*/
    /*                     Getters and Setters                      */
    /*==============================================================*/
    double       get_roll                    (void);
    void        set_roll                    (int roll);

    double       get_pitch                   (void);
    void        set_pitch                   (int pitch);

    double       get_yaw                     (void);
    void        set_yaw                     (int yaw);

    double       get_beta                    (void);
    void        set_beta                    (double beta);

    double*      get_q                       (void);
    void        set_q                       (double q_0, double q_1, double q_2, double q_3);

private:

    /*==============================================================*/
    /*                      Private Variables                       */
    /*==============================================================*/
    double      beta = 0.5;
    double      q0 = 1.0 , q1 = 0.0 , q2 = 0.0, q3 = 0.0;
    int         roll = 0, pitch = 0, yaw = 0;
    double      time_interval = 0.1;

    /*==============================================================*/
    /*                      Private Functions                       */
    /*==============================================================*/
    //ret val   function name               parameters
    double      radians_to_degrees          (double radian_val);
    double      degrees_to_radians          (double degree_val);
    double      invSqrt                     (double x);
};

#endif //MADGWICK_TEST_AHRS_H
