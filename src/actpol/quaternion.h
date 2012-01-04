
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef double Quaternion[4];

// q = a + b
static inline void
Quaternion_add(Quaternion q, const Quaternion a, const Quaternion b)
{
    q[0] = a[0] + b[0];
    q[1] = a[1] + b[1];
    q[2] = a[2] + b[2];
    q[3] = a[3] + b[3];
}

// q = q*
static inline void
Quaternion_conj(Quaternion q)
{
    q[1] = -q[1];
    q[2] = -q[2];
    q[3] = -q[3];
}

// q = a
static inline void
Quaternion_copy(Quaternion q, const Quaternion a)
{
    q[0] = a[0];
    q[1] = a[1];
    q[2] = a[2];
    q[3] = a[3];
}

// q = (1,0,0,0)
static inline void
Quaternion_identity(Quaternion q)
{
    q[0] = 1.;
    q[1] = 0.;
    q[2] = 0.;
    q[3] = 0.;
}

// q = 1/q = q*/|q|^2
void Quaternion_inv(Quaternion q);

// q = a * b
void Quaternion_mul(Quaternion q, const Quaternion a, const Quaternion b);

// q = {w,x,y,z}
static inline void
Quaternion_new(Quaternion q, double w, double x, double y, double z)
{
    q[0] = w;
    q[1] = x;
    q[2] = y;
    q[3] = z;
}

// |q|
double Quaternion_norm(const Quaternion q);

// q = (rotate by angle around arbitrary vector v)
void Quaternion_rot(Quaternion q, double angle, const double v[3]);

// q = R_i(angle)
void Quaternion_r1(Quaternion q, double angle);
void Quaternion_r2(Quaternion q, double angle);
void Quaternion_r3(Quaternion q, double angle);

// q = R_i(angle) * q
void Quaternion_r1_mul(double angle, Quaternion q);
void Quaternion_r2_mul(double angle, Quaternion q);
void Quaternion_r3_mul(double angle, Quaternion q);

// q = q * scale
static inline void
Quaternion_scale(Quaternion q, double scale)
{
    q[0] *= scale;
    q[1] *= scale;
    q[2] *= scale;
    q[3] *= scale;
}

// q = a - b
static inline void
Quaternion_sub(Quaternion q, const Quaternion a, const Quaternion b)
{
    q[0] = a[0] - b[0];
    q[1] = a[1] - b[1];
    q[2] = a[2] - b[2];
    q[3] = a[3] - b[3];
}

void Quaternion_to_matrix(const Quaternion q, double mat[3][3]);

//

typedef struct
{
    Quaternion q0, q1;
    double alpha, sin_alpha;
}
QuaternionSlerp;

void
QuaternionSlerp_init(QuaternionSlerp *slerp, const Quaternion a, const Quaternion b);

void
QuaternionSlerp_interpolate(const QuaternionSlerp *slerp, double t, Quaternion z);

#ifdef __cplusplus
}
#endif

