package pbrt

import "math"

func CoordinateSystem(v1, v2, v3 *Vector3) {
	if math.Abs(v1.x) > math.Abs(v1.y) {
		v := v1.x * v1.x + v1.z * v1.z
		*v2 = Vector3{-v1.z, 0, v1.x}
		v2.DivAssign(&Vector3{v, v, v})
	} else {
		v := v1.y * v1.y + v1.z * v1.z
		*v2 = Vector3{0, v1.z, -v1.y}
		v2.DivAssign(&Vector3{v, v, v})
	}
	*v3 = *v1.Cross(v2)
}

func SphericalDirection(sinTheta, cosTheta, phi float64) *Vector3 {
	return &Vector3{sinTheta * math.Cos(phi), sinTheta * math.Sin(phi), cosTheta}
}

func SphericalDirectionXYZ(sinTheta, cosTheta, phi float64, x, y, z *Vector3) *Vector3 {
	return x.MulScalar(sinTheta * math.Cos(phi)).
		Add(y.MulScalar(sinTheta * math.Sin(phi))).
		Add(z.MulScalar(cosTheta))
}