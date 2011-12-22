
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef double Quaternion[4];

void Quaternion_new(Quaternion q, double w, double x, double y, double z);

// q = a + b
void Quaternion_add(Quaternion q, Quaternion a, Quaternion b);

// q = q*
void Quaternion_conj(Quaternion q);

// q = a
void Quaternion_copy(Quaternion q, Quaternion a);

// q = (1,0,0,0)
void Quaternion_identity(Quaternion q);

// q = 1/q = q*/|q|^2
void Quaternion_inv(Quaternion q);

// q = a * b
void Quaternion_mul(Quaternion q, Quaternion a, Quaternion b);

// |q|
double Quaternion_norm(Quaternion q);

// q = (rotate by angle around arbitrary vector v)
void Quaternion_rot(Quaternion q, double angle, double v[3]);

// q = R_i(angle)
void Quaternion_r1(Quaternion q, double angle);
void Quaternion_r2(Quaternion q, double angle);
void Quaternion_r3(Quaternion q, double angle);

// q = R_i(angle) * q
void Quaternion_r1_mul(double angle, Quaternion q);
void Quaternion_r2_mul(double angle, Quaternion q);
void Quaternion_r3_mul(double angle, Quaternion q);

// q = q * scale
void Quaternion_scale(Quaternion q, double scale);

// q = a - b
void Quaternion_sub(Quaternion q, Quaternion a, Quaternion b);

void Quaternion_to_matrix(Quaternion q, double mat[3][3]);

#ifdef __cplusplus
}
#endif

