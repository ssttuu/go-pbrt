package pbrt

import "github.com/stupschwartz/go-pbrt/pkg/math"

type Quaternion struct {
	v *Vector3f
	w float64
}

func (q *Quaternion) Add(other *Quaternion) *Quaternion {
	return &Quaternion{
		v: q.v.Add(other.v),
		w: q.w * other.w,
	}
}

func (q *Quaternion) Sub(other *Quaternion) *Quaternion {
	return &Quaternion{
		v: q.v.Sub(other.v),
		w: q.w - other.w,
	}
}

func (q *Quaternion) MulScalar(other float64) *Quaternion {
	return &Quaternion{
		v: q.v.MulScalar(other),
		w: q.w * other,
	}
}

func (q *Quaternion) DivScalar(other float64) *Quaternion {
	return &Quaternion{
		v: q.v.DivScalar(other),
		w: q.w / other,
	}
}

func (q *Quaternion) Dot(other *Quaternion) float64 {
	return q.v.Dot(other.v) + q.w*other.w
}

func (q *Quaternion) Normalized() *Quaternion {
	return q.DivScalar(math.Sqrt(q.Dot(q)))
}

func (q *Quaternion) ToTransform() *Transform {
	xx, yy, zz := q.v.X*q.v.X, q.v.Y*q.v.Y, q.v.Z*q.v.Z
	xy, xz, yz := q.v.X*q.v.Y, q.v.X*q.v.Z, q.v.Y*q.v.Z
	wx, wy, wz := q.v.X*q.w, q.v.Y*q.w, q.v.Z*q.w

	var m *Matrix4x4
	m[0][0] = 1 - 2*(yy+zz)
	m[0][1] = 2 * (xy + wz)
	m[0][2] = 2 * (xz - wy)
	m[1][0] = 2 * (xy - wz)
	m[1][1] = 1 - 2*(xx+zz)
	m[1][2] = 2 * (yz + wx)
	m[2][0] = 2 * (xz + wy)
	m[2][1] = 2 * (yz - wx)
	m[2][2] = 1 - 2*(xx+yy)

	// Transpose since we are left-handed
	return &Transform{m.Transpose(), m}
}

func Slerp(t float64, q1, q2 *Quaternion) *Quaternion {
	cosTheta := q1.Dot(q2)
	if cosTheta > 0.9995 {
		return q1.MulScalar(1 - t).Add(q2.MulScalar(t)).Normalized()
	}

	theta := math.Acos(math.Clamp(cosTheta, -1, 1))
	thetap := theta * t
	qperp := q2.Sub(q1.MulScalar(cosTheta))
	return q1.MulScalar(math.Cos(thetap)).Add(qperp.MulScalar(math.Sin(thetap)))
}
