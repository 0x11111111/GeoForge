#ifndef ROTATION_H
#define ROTATION_H

#include <algorithm>
#include <cmath>
#include <limits>

//////////////////////////////////////////////////////////////////
// 旋转所需的公式

// 点乘和叉乘和归一化 

template<typename T>
inline T DotProduct(const T x[3], const T y[3])
{
    return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

template<typename T>
inline void CrossProduct(const T x[3], const T y[3], T result[3])
{
    result[0] = x[1] * y[2] - x[2] * y[1];
    result[1] = x[2] * y[0] - x[0] * y[2];
    result[2] = x[0] * y[1] - x[1] * y[0];
}

template<typename T>
inline void Normalize(T* vec)
{
    T length = sqrt(DotProduct(vec, vec));
    vec[0] /= length;
    vec[1] /= length;
    vec[2] /= length;
}

template<typename T>
inline T GetLength(T* vec)
{
    return sqrt(DotProduct(vec, vec));
}
//////////////////////////////////////////////////////////////////


// 将角度转换成四元数
template<typename T>
inline void AngleAxisToQuaternion(const T* angle_axis, T* quaternion)
{
    const T& a0 = angle_axis[0];
    const T& a1 = angle_axis[1];
    const T& a2 = angle_axis[2];
    const T theta_squared = a0 * a0 + a1 * a1 + a2 * a2;

    if (theta_squared > T(std::numeric_limits<double>::epsilon()))
    {
        const T theta = sqrt(theta_squared);
        const T half_theta = theta * T(0.5);
        const T k = sin(half_theta) / theta;
        quaternion[0] = cos(half_theta);
        quaternion[1] = a0 * k;
        quaternion[2] = a1 * k;
        quaternion[3] = a2 * k;
    }
    else
    {
        // 如果theta的平方为0
        const T k(0.5);
        quaternion[0] = T(1.0);
        quaternion[1] = a0 * k;
        quaternion[2] = a1 * k;
        quaternion[3] = a2 * k;
    }
}

template<typename T>
inline void QuaternionToAngleAxis(const T* quaternion, T* angle_axis)
{
    const T& q1 = quaternion[1];
    const T& q2 = quaternion[2];
    const T& q3 = quaternion[3];
    const T sin_squared_theta = q1 * q1 + q2 * q2 + q3 * q3;

    // 对于表示非零旋转的四元数,转换在数值上是稳定的
    if (sin_squared_theta > T(std::numeric_limits<double>::epsilon()))
    {
        const T sin_theta = sqrt(sin_squared_theta);
        const T& cos_theta = quaternion[0];

        // 如果 cos_theta 为负,则 theta 大于 pi/2,这意味着 2 * theta 的角度轴矢量的角度将大于 pi...

        const T two_theta = T(2.0) * ((cos_theta < 0.0)
            ? atan2(-sin_theta, -cos_theta)
            : atan2(sin_theta, cos_theta));
        const T k = two_theta / sin_theta;

        angle_axis[0] = q1 * k;
        angle_axis[1] = q2 * k;
        angle_axis[2] = q3 * k;
    }
    else
    {
        // 对于近似于零的旋转, sqrt() 将出现NaN ,因为参数为零. 通过使用泰勒级数进行逼近,并在一阶处截断,当使用 Jets 时,值和一阶导数将被正确计算.
        const T k(2.0);
        angle_axis[0] = q1 * k;
        angle_axis[1] = q2 * k;
        angle_axis[2] = q3 * k;
    }
}

template<typename T>
inline void AngleAxisRotatePoint(const T angle_axis[3], const T pt[3], T result[3])
{
    const T theta2 = DotProduct(angle_axis, angle_axis);
    if (theta2 > T(std::numeric_limits<double>::epsilon()))
    {
        // 非零时,使用罗德里格斯公式
        //
        //   result = pt costheta +
        //            (w x pt) * sintheta +
        //            w (w . pt) (1 - costheta)
        //
        //如果 angle_axis 向量的范数大于零,我们要只计算平方根,否则会被零除.
        //
        const T theta = sqrt(theta2);
        const T costheta = cos(theta);
        const T sintheta = sin(theta);
        const T theta_inverse = 1.0 / theta;

        const T w[3] = { angle_axis[0] * theta_inverse,
                        angle_axis[1] * theta_inverse,
                        angle_axis[2] * theta_inverse };
        T w_cross_pt[3];
        CrossProduct(w, pt, w_cross_pt);

        const T tmp = DotProduct(w, pt) * (T(1.0) - costheta);
        //    (w[0] * pt[0] + w[1] * pt[1] + w[2] * pt[2]) * (T(1.0) - costheta);

        result[0] = pt[0] * costheta + w_cross_pt[0] * sintheta + w[0] * tmp;
        result[1] = pt[1] * costheta + w_cross_pt[1] * sintheta + w[1] * tmp;
        result[2] = pt[2] * costheta + w_cross_pt[2] * sintheta + w[2] * tmp;
    }
    else
    {
        // Near zero, the first order Taylor approximation of the rotation
        // matrix R corresponding to a vector w and angle w is
        //
        //   R = I + hat(w) * sin(theta)
        //
        // But sintheta ~ theta and theta * w = angle_axis, which gives us
        //
        //  R = I + hat(w)
        //
        // and actually performing multiplication with the point pt, gives us
        // R * pt = pt + w x pt.
        //
        // Switching to the Taylor expansion near zero provides meaningful
        // derivatives when evaluated using Jets.
        //
        // Explicitly inlined evaluation of the cross product for
        // performance reasons.
        /*const T w_cross_pt[3] = { angle_axis[1] * pt[2] - angle_axis[2] * pt[1],
                                  angle_axis[2] * pt[0] - angle_axis[0] * pt[2],
                                  angle_axis[0] * pt[1] - angle_axis[1] * pt[0] };*/
        T w_cross_pt[3];
        CrossProduct(angle_axis, pt, w_cross_pt);

        result[0] = pt[0] + w_cross_pt[0];
        result[1] = pt[1] + w_cross_pt[1];
        result[2] = pt[2] + w_cross_pt[2];
    }
}

#endif