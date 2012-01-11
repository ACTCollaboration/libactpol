
#include <assert.h>

#include "actpol/math.h"
#include "actpol/quaternion.h"

void
Quaternion_inv(Quaternion q)
{
    double norm2 = Quaternion_norm2(q);
    q[0] = q[0]/norm2;
    q[1] = -q[1]/norm2;
    q[2] = -q[2]/norm2;
    q[3] = -q[3]/norm2;
}

double
Quaternion_norm(const Quaternion q)
{
    return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
}

void
Quaternion_rot(Quaternion q, double angle, const double v[3])
{
    double angle_2 = 0.5*angle;
    double s = sin(angle_2);
    double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    assert(norm > 0.);
    q[0] = cos(angle_2);
    q[1] = s*v[0]/norm;
    q[2] = s*v[1]/norm;
    q[3] = s*v[2]/norm;
}

void
Quaternion_r1(Quaternion q, double angle)
{
    double angle_2 = 0.5*angle;
    q[0] = cos(angle_2);
    q[1] = sin(angle_2);
    q[2] = 0.;
    q[3] = 0.;
}

void
Quaternion_r2(Quaternion q, double angle)
{
    double angle_2 = 0.5*angle;
    q[0] = cos(angle_2);
    q[1] = 0.;
    q[2] = sin(angle_2);
    q[3] = 0.;
}

void
Quaternion_r3(Quaternion q, double angle)
{
    double angle_2 = 0.5*angle;
    q[0] = cos(angle_2);
    q[1] = 0.;
    q[2] = 0.;
    q[3] = sin(angle_2);
}

void
Quaternion_r1_mul(double angle, Quaternion q)
{
    Quaternion a, b;
    Quaternion_r1(a, angle);
    Quaternion_copy(b, q);

    q[0] = a[0]*b[0] - a[1]*b[1]; // - a[2]*b[2] - a[3]*b[3];
    q[1] = a[0]*b[1] + a[1]*b[0]; // + a[2]*b[3] - a[3]*b[2];
    q[2] = a[0]*b[2] - a[1]*b[3]; // + a[2]*b[0] + a[3]*b[1];
    q[3] = a[0]*b[3] + a[1]*b[2]; // - a[2]*b[1] + a[3]*b[0];
}

void
Quaternion_r2_mul(double angle, Quaternion q)
{
    Quaternion b;
    double angle_2 = 0.5*angle;
    double c = cos(angle_2), s = sin(angle_2);
    Quaternion_copy(b, q);

    q[0] = c*b[0] - s*b[2];
    q[1] = c*b[1] + s*b[3];
    q[2] = c*b[2] + s*b[0];
    q[3] = c*b[3] - s*b[1];

    /*
    Quaternion_r2(angle, a);

    q[0] = a[0]*b[0] - a[2]*b[2]; //- a[3]*b[3];
    q[1] = a[0]*b[1] + a[2]*b[3]; //- a[3]*b[2];
    q[2] = a[0]*b[2] + a[2]*b[0]; //+ a[3]*b[1];
    q[3] = a[0]*b[3] - a[2]*b[1]; //+ a[3]*b[0];
    */
}

void
Quaternion_r3_mul(double angle, Quaternion q)
{
    Quaternion a, b;
    Quaternion_r3(a, angle);
    Quaternion_copy(b, q);

    q[0] = a[0]*b[0] /*- a[1]*b[1] - a[2]*b[2]*/ - a[3]*b[3];
    q[1] = a[0]*b[1] /*+ a[1]*b[0] + a[2]*b[3]*/ - a[3]*b[2];
    q[2] = a[0]*b[2] /*- a[1]*b[3] + a[2]*b[0]*/ + a[3]*b[1];
    q[3] = a[0]*b[3] /*+ a[1]*b[2] - a[2]*b[1]*/ + a[3]*b[0];
}

void
Quaternion_to_matrix(const Quaternion q, double mat[3][3])
{
    Quaternion u;
    Quaternion_copy(u, q);
    Quaternion_unit(u);

    double a2 = u[0]*u[0], b2 = u[1]*u[1], c2 = u[2]*u[2], d2 = u[3]*u[3];
    mat[0][0] = a2 + b2 - c2 - d2;
    mat[1][1] = a2 - b2 + c2 - d2;
    mat[2][2] = a2 - b2 - c2 + d2;
    mat[0][1] = 2.*(u[1]*u[2] - u[0]*u[3]);
    mat[0][2] = 2.*(u[1]*u[3] + u[0]*u[2]);
    mat[1][2] = 2.*(u[2]*u[3] - u[0]*u[1]);
    mat[1][0] = 2.*(u[1]*u[2] + u[0]*u[3]);
    mat[2][0] = 2.*(u[1]*u[3] - u[0]*u[2]);
    mat[2][1] = 2.*(u[2]*u[3] + u[0]*u[1]);
}

void
Quaternion_to_matrix_col1(const Quaternion u, double col1[3])
{
    // make sure quaternion is normalized before calling
    double a2 = u[0]*u[0], b2 = u[1]*u[1], c2 = u[2]*u[2], d2 = u[3]*u[3];
    col1[0] = a2 + b2 - c2 - d2;
    col1[1] = 2.*(u[1]*u[2] + u[0]*u[3]);
    col1[2] = 2.*(u[1]*u[3] - u[0]*u[2]);
}

void
Quaternion_to_matrix_col2(const Quaternion u, double col2[3])
{
    // make sure quaternion is normalized before calling
    double a2 = u[0]*u[0], b2 = u[1]*u[1], c2 = u[2]*u[2], d2 = u[3]*u[3];
    col2[0] = 2.*(u[1]*u[2] - u[0]*u[3]);
    col2[1] = a2 - b2 + c2 - d2;
    col2[2] = 2.*(u[2]*u[3] + u[0]*u[1]);
}

void
Quaternion_to_matrix_col3(const Quaternion u, double col3[3])
{
    // make sure quaternion is normalized before calling
    double a2 = u[0]*u[0], b2 = u[1]*u[1], c2 = u[2]*u[2], d2 = u[3]*u[3];
    col3[0] = 2.*(u[1]*u[3] + u[0]*u[2]);
    col3[1] = 2.*(u[2]*u[3] - u[0]*u[1]);
    col3[2] = a2 - b2 - c2 + d2;
}

void
Quaternion_unit(Quaternion q)
{
    double norm2 = Quaternion_norm2(q);
    double invnorm = invsqrt(norm2);
    Quaternion_scale(q, invnorm);
}

void
QuaternionSlerp_init(QuaternionSlerp *slerp, const Quaternion a, const Quaternion b)
{
    double cos_alpha = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
    slerp->alpha = acos(cos_alpha);
    slerp->sin_alpha = sqrt(1. - cos_alpha*cos_alpha);
    Quaternion_copy(slerp->q0, a);
    Quaternion_copy(slerp->q1, b);
}

void
QuaternionSlerp_interpolate(const QuaternionSlerp *slerp, double t, Quaternion q)
{
    double s0 = sin((1.-t)*slerp->alpha)/slerp->sin_alpha;
    double s1 = sin(t*slerp->alpha)/slerp->sin_alpha;
    for (int i = 0; i != 4; ++i)
        q[i] = s0*slerp->q0[i] + s1*slerp->q1[i];
}

