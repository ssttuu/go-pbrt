package pbrt

import (
	"math"
)

type Point2i = XYInt64
type Point2f = XYFloat64

type Vector2i = XYInt64
type Vector2f = XYFloat64

type Point3i = XYZInt64
type Point3f = XYZFloat64

type Vector3i = XYZInt64
type Vector3f = XYZFloat64

type Normal3i = XYZInt64
type Normal3f = XYZFloat64

func CoordinateSystem(v1 *Vector3f) (v2, v3 *Vector3f) {
	if math.Abs(v1.X) > math.Abs(v1.Y) {
		v := v1.X*v1.X + v1.Z*v1.Z
		v2 = &Vector3f{-v1.Z, 0, v1.X}
		v2.DivAssign(&Vector3f{v, v, v})
	} else {
		v := v1.Y*v1.Y + v1.Z*v1.Z
		v2 = &Vector3f{0, v1.Z, -v1.Y}
		v2.DivAssign(&Vector3f{v, v, v})
	}
	v3 = v1.Cross(v2)

	return v2, v3
}

func SphericalDirection(sinTheta, cosTheta, phi float64) *Vector3f {
	return &Vector3f{sinTheta * math.Cos(phi), sinTheta * math.Sin(phi), cosTheta}
}

func SphericalDirectionXYZ(sinTheta, cosTheta, phi float64, x, y, z *Vector3f) *Vector3f {
	return x.MulScalar(sinTheta * math.Cos(phi)).
		Add(y.MulScalar(sinTheta * math.Sin(phi))).
		Add(z.MulScalar(cosTheta))
}

func MinPoint(p0, p1 *Point3f) *Point3f {
	return &Point3f{
		X: math.Min(p0.X, p1.X),
		Y: math.Min(p0.Y, p1.Y),
		Z: math.Min(p0.Z, p1.Z),
	}
}

func MaxPoint(p0, p1 *Point3f) *Point3f {
	return &Point3f{
		X: math.Max(p0.X, p1.X),
		Y: math.Max(p0.Y, p1.Y),
		Z: math.Max(p0.Z, p1.Z),
	}
}

func Quadratic(a, b, c float64) (hasRoot bool, t0, t1 float64) {
	discrim := b*b - 4.0*a*c
	if discrim < 0.0 {
		return false, 0, 0
	}

	rootDiscrim := math.Sqrt(discrim)

	var q float64
	if b < 0 {
		q = -0.5 * (b - rootDiscrim)
	} else {
		q = -0.5 * (b + rootDiscrim)
	}
	t0 = q / a
	t1 = c / q

	if t0 > t1 {
		t0, t1 = t1, t0
	}
	return true, t0, t1
}

func FaceForward(n1, n2 *Normal3f) *Normal3f {
	if n1.Dot(n2) < 0.0 {
		return n1.MulScalar(-1)
	}
	return n1
}
