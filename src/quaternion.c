
#include <assert.h>
#include <math.h>

#include "actpol/quaternion.h"

void
Quaternion_add(Quaternion q, Quaternion a, Quaternion b)
{
    q[0] = a[0] + b[0];
    q[1] = a[1] + b[1];
    q[2] = a[2] + b[2];
    q[3] = a[3] + b[3];
}

void
Quaternion_conj(Quaternion q)
{
    q[1] = -q[1];
    q[2] = -q[2];
    q[3] = -q[3];
}

void
Quaternion_copy(Quaternion q, Quaternion a)
{
    q[0] = a[0];
    q[1] = a[1];
    q[2] = a[2];
    q[3] = a[3];
}

void
Quaternion_identity(Quaternion q)
{
    q[0] = 1.;
    q[1] = 0.;
    q[2] = 0.;
    q[3] = 0.;
}

void
Quaternion_inv(Quaternion q)
{
    double norm = Quaternion_norm(q);
    double norm2 = norm*norm;
    q[0] = q[0]/norm2;
    q[1] = -q[1]/norm2;
    q[2] = -q[2]/norm2;
    q[3] = -q[3]/norm2;
}

void
Quaternion_mul(Quaternion q, Quaternion a, Quaternion b)
{
    q[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
    q[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
    q[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
    q[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
}

void
Quaternion_new(Quaternion q, double w, double x, double y, double z)
{
    q[0] = w;
    q[1] = x;
    q[2] = y;
    q[3] = z;
}

double
Quaternion_norm(Quaternion q)
{
    return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
}

void
Quaternion_scale(Quaternion q, double scale)
{
    q[0] *= scale;
    q[1] *= scale;
    q[2] *= scale;
    q[3] *= scale;
}

void
Quaternion_sub(Quaternion q, Quaternion a, Quaternion b)
{
    q[0] = a[0] - b[0];
    q[1] = a[1] - b[1];
    q[2] = a[2] - b[2];
    q[3] = a[3] - b[3];
}

void
Quaternion_rot(Quaternion q, double angle, double v[3])
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
Quaternion_to_matrix(Quaternion q, double mat[3][3])
{
    Quaternion u;
    Quaternion_copy(u, q);
    double norm = Quaternion_norm(u);
    Quaternion_scale(u, 1./norm);

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

