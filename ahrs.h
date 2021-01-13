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
    Madgwick(float beta_val);
    ~Madgwick();

    /*==============================================================*/
    /*                       Public Functions                       */
    /*==============================================================*/
    //ret val   function name               parameters
    void        MadgwickAHRSupdate          (float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz);
    void        MadgwickAHRSupdateIMU       (float gx, float gy, float gz, float ax, float ay, float az);
//    void        compute_orientation         (float *q);
//    float*      quaternionMul               (float *q_1, float *q_2);
//    float*      get_accel_jacobian          (float *q);
//    float*      get_accel_function          (float *q, float *a);
//    void        normalizeq                  (float *q);
//    void        update_roll_and_pitch       ();
//    void        get_mag_jacob               ();
//    void        get_mag_func                ();
//    void        get_rotation_mat            ();
//    void        update_roll_pitch_yaw       ();

    /*==============================================================*/
    /*                     Getters and Setters                      */
    /*==============================================================*/
    float       get_roll                    (void);
    void        set_roll                    (int roll);

    float       get_pitch                   (void);
    void        set_pitch                   (int pitch);

    float       get_yaw                     (void);
    void        set_yaw                     (int yaw);

    float       get_beta                    (void);
    void        set_beta                    (float beta);

    float*      get_q                       (void);
    void        set_q                       (float q_0, float q_1, float q_2, float q_3);

private:

    /*==============================================================*/
    /*                      Private Variables                       */
    /*==============================================================*/
    float      beta = 0.5;
    float      q0 = 1.0 , q1 = 0.0 , q2 = 0.0, q3 = 0.0;
    int        roll = 0, pitch = 0, yaw = 0;
    float      sampleFreq = 10.0;

    /*==============================================================*/
    /*                      Private Functions                       */
    /*==============================================================*/
    //ret val   function name               parameters
    float      radians_to_degrees          (float radian_val);
    float      degrees_to_radians          (float degree_val);
    float      invSqrt                     (float x);
};

#endif //MADGWICK_TEST_AHRS_H
