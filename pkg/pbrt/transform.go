package pbrt

import (
	"math"
	"errors"
)

var ErrSingularMatrix = errors.New("singular matrix in matrix invert")
var ErrSameDirectionVectors = errors.New("singular matrix in matrix invert")

type Matrix4x4 [4][4]float64

func NewMatrix4x4() Matrix4x4 {
	return Matrix4x4{
		{1.0, 0.0, 0.0, 0.0},
		{0.0, 1.0, 0.0, 0.0},
		{0.0, 0.0, 1.0, 0.0},
		{0.0, 0.0, 0.0, 1.0},
	}
}

func (m *Matrix4x4) Equals(other *Matrix4x4) bool {
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if m[i][j] != other[i][j] {
				return false
			}
		}
	}
	return true
}

func (m *Matrix4x4) NotEquals(other *Matrix4x4) bool {
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if m[i][j] != other[i][j] {
				return true
			}
		}
	}
	return false
}

func (m *Matrix4x4) Transpose() *Matrix4x4 {
	return &Matrix4x4{
		{m[0][0], m[1][0], m[2][0], m[3][0]},
		{m[0][1], m[1][1], m[2][1], m[3][1]},
		{m[0][2], m[1][2], m[2][2], m[3][2]},
		{m[0][3], m[1][3], m[2][3], m[3][3]},
	}
}

func (m *Matrix4x4) Mul(other *Matrix4x4) *Matrix4x4 {
	var r *Matrix4x4
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			r[i][j] = m[i][0]*other[0][j] + m[i][1]*other[1][j] + m[i][2]*other[2][j] + m[i][3]*m[3][j]
		}
	}
	return r
}

func (m *Matrix4x4) Inverse() (*Matrix4x4, error) {
	var indxc, indxr, ipiv [4]int

	minv := *m

	for i := 0; i < 4; i++ {
		irow, icol := 0, 0
		big := 0.0

		// Choose pivot
		for j := 0; j < 4; j++ {
			if ipiv[j] != 1 {
				for k := 0; k < 4; k++ {
					if ipiv[k] == 0 {
						if math.Abs(minv[j][k]) >= big {
							big = float64(math.Abs(minv[j][k]))
							irow = j
							icol = k
						}
					} else if ipiv[k] > 1 {
						return nil, ErrSingularMatrix
					}
				}
			}
		}

		ipiv[icol]++

		// Swap rows _irow_ and _icol_ for pivot
		if irow != icol {
			for k := 0; k < 4; k++ {
				minv[irow][k], minv[icol][k] = minv[icol][k], minv[irow][k]
			}
		}

		indxr[i] = irow
		indxc[i] = icol
		if minv[icol][icol] == 0.0 {
			return nil, ErrSingularMatrix
		}

		// Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
		pivinv := 1.0 / minv[icol][icol]
		minv[icol][icol] = 1.0
		for j := 0; j < 4; j++ {
			minv[icol][j] *= pivinv
		}

		// Subtract this row from others to zero out their columns
		for j := 0; j < 4; j++ {
			if j != icol {
				save := minv[j][icol]
				minv[j][icol] = 0.0
				for k := 0; k < 4; k++ {
					minv[j][k] -= minv[icol][k] * save
				}
			}
		}
	}

	// Swap columns to reflect permutation
	for j := 3; j >= 0; j-- {
		if indxr[j] != indxc[j] {
			for k := 0; k < 4; k++ {
				minv[k][indxr[j]], minv[k][indxc[j]] = minv[k][indxc[j]], minv[k][indxr[j]]
			}
		}
	}

	return &minv, nil
}

type Transform struct {
	matrix, matrixInverse *Matrix4x4
}

func NewTransform(m *Matrix4x4) *Transform {
	mInv, err := m.Inverse()
	if err != nil {
		// TODO
		return nil
	}
	return &Transform{
		matrix:        m,
		matrixInverse: mInv,
	}
}

//func TransformVector(v *Vector3f) *Transform {
//	return
//}
//
//func NewTransformFromRay(r *Ray) *Transform {
//
//}

func (t *Transform) IsIdentity() bool {
	return t.matrix[0][0] == 1.0 && t.matrix[0][1] == 0.0 && t.matrix[0][2] == 0.0 && t.matrix[0][3] == 0.0 &&
		t.matrix[1][0] == 0.0 && t.matrix[1][1] == 1.0 && t.matrix[1][2] == 0.0 && t.matrix[1][3] == 0.0 &&
		t.matrix[2][0] == 0.0 && t.matrix[2][1] == 0.0 && t.matrix[2][2] == 1.0 && t.matrix[2][3] == 0.0 &&
		t.matrix[3][0] == 0.0 && t.matrix[3][1] == 0.0 && t.matrix[3][2] == 0.0 && t.matrix[3][3] == 1.0

}

func (t *Transform) Inverse() *Transform {
	return &Transform{t.matrixInverse, t.matrix}
}

func (t *Transform) Mul(other *Transform) *Transform {
	return &Transform{
		matrix:        t.matrix.Mul(other.matrix),
		matrixInverse: t.matrixInverse.Mul(other.matrixInverse),
	}
}

func (t *Transform) TransformPoint(p *Point3f) *Point3f {
	xp := t.matrix[0][0]*p.X + t.matrix[0][1]*p.Y + t.matrix[0][2]*p.Z + t.matrix[0][3]
	yp := t.matrix[1][0]*p.X + t.matrix[1][1]*p.Y + t.matrix[1][2]*p.Z + t.matrix[1][3]
	zp := t.matrix[2][0]*p.X + t.matrix[2][1]*p.Y + t.matrix[2][2]*p.Z + t.matrix[2][3]
	wp := t.matrix[3][0]*p.X + t.matrix[3][1]*p.Y + t.matrix[3][2]*p.Z + t.matrix[3][3]

	pTransformed := &Point3f{xp, yp, zp}
	if wp == 1.0 {
		return pTransformed
	}

	return pTransformed.DivScalar(wp)
}

func (t *Transform) TransformVector(v *Vector3f) *Vector3f {
	return &Vector3f{
		t.matrix[0][0]*v.X + t.matrix[0][1]*v.Y + t.matrix[0][2]*v.Z,
		t.matrix[1][0]*v.X + t.matrix[1][1]*v.Y + t.matrix[1][2]*v.Z,
		t.matrix[2][0]*v.X + t.matrix[2][1]*v.Y + t.matrix[2][2]*v.Z,
	}
}

func (t *Transform) TransformNormal(n *Vector3f) *Vector3f {
	return &Vector3f{
		t.matrixInverse[0][0]*n.X + t.matrixInverse[1][0]*n.Y + t.matrixInverse[2][0]*n.Z,
		t.matrixInverse[0][1]*n.X + t.matrixInverse[1][1]*n.Y + t.matrixInverse[2][1]*n.Z,
		t.matrixInverse[0][2]*n.X + t.matrixInverse[1][2]*n.Y + t.matrixInverse[2][2]*n.Z,
	}
}

func (t *Transform) TransformRay(r *Ray) *Ray {
	origin := t.TransformPoint(r.origin)
	direction := t.TransformVector(r.direction)
	return &Ray{origin, direction, r.tMax, r.time, r.medium}
}

func (t *Transform) TransformSurfaceInteraction(si *SurfaceInteraction) *SurfaceInteraction {
	ret := *si
	ret.point = t.TransformPoint(si.point)

	ret.normal = t.TransformNormal(si.normal).Normalized()
	ret.wo = t.TransformVector(si.wo).Normalized()
	ret.time = si.time
	ret.mediumAccessor = si.mediumAccessor
	ret.uv = si.uv
	ret.shape = si.shape
	ret.dpdu = t.TransformVector(ret.dpdu)
	ret.dpdv = t.TransformVector(ret.dpdv)
	ret.dndu = t.TransformNormal(ret.dndu)
	ret.dndv = t.TransformNormal(ret.dndv)
	ret.shading.normal = t.TransformNormal(si.shading.normal)
	ret.shading.dpdu = t.TransformVector(si.shading.dpdu)
	ret.shading.dpdv = t.TransformVector(si.shading.dpdv)
	ret.shading.dndu = t.TransformNormal(si.shading.dndu)
	ret.shading.dndv = t.TransformNormal(si.shading.dndv)
	ret.dudx = si.dudx
	ret.dvdx = si.dvdx
	ret.dudy = si.dudy
	ret.dvdy = si.dvdy
	//ret.dpdx = t.TransformVector(si.dpdx)
	//ret.dpdy = t.TransformVector(si.dpdy)
	ret.bsdf = si.bsdf
	ret.primitive = si.primitive
	ret.shading.normal = FaceForward(ret.shading.normal, ret.normal)
	ret.faceIndex = si.faceIndex
	return &ret
}

func (t *Transform) TransformBounds(b *Bounds3) *Bounds3 {
	ret := &Bounds3{min: t.TransformPoint(b.min), max: t.TransformPoint(b.max)}
	ret.Union(t.TransformPoint(&Point3f{b.max.X, b.min.Y, b.min.Z}))
	ret.Union(t.TransformPoint(&Point3f{b.min.X, b.max.Y, b.min.Z}))
	ret.Union(t.TransformPoint(&Point3f{b.min.X, b.min.Y, b.max.Z}))
	ret.Union(t.TransformPoint(&Point3f{b.min.X, b.max.Y, b.max.Z}))
	ret.Union(t.TransformPoint(&Point3f{b.max.X, b.max.Y, b.min.Z}))
	ret.Union(t.TransformPoint(&Point3f{b.max.X, b.min.Y, b.max.Z}))
	return ret
}

func Translate(delta *Vector3f) *Transform {
	return &Transform{
		matrix: &Matrix4x4{
			{1, 0, 0, delta.X},
			{0, 1, 0, delta.Y},
			{0, 0, 1, delta.Z},
			{0, 0, 0, 1},
		},
		matrixInverse: &Matrix4x4{
			{1, 0, 0, -delta.X},
			{0, 1, 0, -delta.Y},
			{0, 0, 1, -delta.Z},
			{0, 0, 0, 1},
		},
	}
}

func Scale(x, y, z float64) *Transform {
	return &Transform{
		matrix: &Matrix4x4{
			{x, 0, 0, 0},
			{0, y, 0, 0},
			{0, 0, z, 0},
			{0, 0, 0, 1},
		},
		matrixInverse: &Matrix4x4{
			{1.0 / x, 0, 0, 0},
			{0, 1.0 / y, 0, 0},
			{0, 0, 1.0 / z, 0},
			{0, 0, 0, 1},
		},
	}
}

func RotateX(theta float64) *Transform {
	sinTheta := math.Sin(Radians(theta))
	cosTheta := math.Cos(Radians(theta))
	m := &Matrix4x4{
		{1, 0, 0, 0},
		{0, cosTheta, -sinTheta, 0},
		{0, sinTheta, cosTheta, 0},
		{0, 0, 0, 1},
	}
	return &Transform{
		matrix:        m,
		matrixInverse: m.Transpose(),
	}
}

func RotateY(theta float64) *Transform {
	sinTheta := math.Sin(Radians(theta))
	cosTheta := math.Cos(Radians(theta))
	m := &Matrix4x4{
		{cosTheta, 0, sinTheta, 0},
		{0, 1, 0, 0},
		{-sinTheta, 0, cosTheta, 0},
		{0, 0, 0, 1},
	}
	return &Transform{
		matrix:        m,
		matrixInverse: m.Transpose(),
	}
}

func RotateZ(theta float64) *Transform {
	sinTheta := math.Sin(Radians(theta))
	cosTheta := math.Cos(Radians(theta))
	m := &Matrix4x4{
		{cosTheta, -sinTheta, 0, 0},
		{sinTheta, cosTheta, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
	}
	return &Transform{
		matrix:        m,
		matrixInverse: m.Transpose(),
	}
}

func Rotate(theta float64, axis Vector3f) *Transform {
	a := axis.Normalized()
	sinTheta := math.Sin(Radians(theta))
	cosTheta := math.Cos(Radians(theta))

	var m *Matrix4x4

	// Compute rotation of first basis vector
	m[0][0] = a.X*a.X + (1-a.X*a.X)*cosTheta
	m[0][1] = a.X*a.Y*(1-cosTheta) - a.Z*sinTheta
	m[0][2] = a.X*a.Z*(1-cosTheta) + a.Y*sinTheta
	m[0][3] = 0

	// Compute rotations of second and third basis vectors
	m[1][0] = a.X*a.Y*(1-cosTheta) + a.Z*sinTheta
	m[1][1] = a.Y*a.Y + (1-a.Y*a.Y)*cosTheta
	m[1][2] = a.Y*a.Z*(1-cosTheta) - a.X*sinTheta
	m[1][3] = 0

	m[2][0] = a.X*a.Z*(1-cosTheta) - a.Y*sinTheta
	m[2][1] = a.Y*a.Z*(1-cosTheta) + a.X*sinTheta
	m[2][2] = a.Z*a.Z + (1-a.Z*a.Z)*cosTheta
	m[2][3] = 0

	return &Transform{matrix: m, matrixInverse: m.Transpose()}
}

func LookAt(pos *Point3f, look *Point3f, up *Vector3f) (*Transform, error) {
	var cameraToWorld *Matrix4x4
	// Initialize fourth column of viewing matrix
	cameraToWorld[0][3] = pos.X
	cameraToWorld[1][3] = pos.Y
	cameraToWorld[2][3] = pos.Z
	cameraToWorld[3][3] = 1

	// Initialize first three columns of viewing matrix
	dir := look.Sub(pos).Normalized()
	if up.Normalized().Cross(dir).Length() == 0 {
		return nil, ErrSameDirectionVectors
	}
	left := up.Normalized().Cross(dir).Normalized()
	newUp := dir.Cross(left)
	cameraToWorld[0][0] = left.X
	cameraToWorld[1][0] = left.Y
	cameraToWorld[2][0] = left.Z
	cameraToWorld[3][0] = 0.
	cameraToWorld[0][1] = newUp.X
	cameraToWorld[1][1] = newUp.Y
	cameraToWorld[2][1] = newUp.Z
	cameraToWorld[3][1] = 0.
	cameraToWorld[0][2] = dir.X
	cameraToWorld[1][2] = dir.Y
	cameraToWorld[2][2] = dir.Z
	cameraToWorld[3][2] = 0.

	inv, err := cameraToWorld.Inverse()
	if err != nil {
		return nil, err
	}
	return &Transform{matrix: inv, matrixInverse: cameraToWorld}, nil
}

func Orthographic(zNear, zFar float64) *Transform {
	return Scale(1, 1, 1/(zFar-zNear)).Mul(Translate(&Vector3f{0, 0, -zNear}))
}

func Perspective(fov, n, f float64) *Transform {
	persp := Matrix4x4{
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, f / (f - n), -f * n / (f - n)},
		{0, 0, 1, 0},
	}

	invTanAng := 1.0 / math.Tan(Radians(fov)/2)
	return Scale(invTanAng, invTanAng, 1).Mul(NewTransform(&persp))
}

type DerivativeTerm struct {
	kc, kx, ky, kz float64
}

func (dt *DerivativeTerm) Eval(p *Point3f) float64 {
	return dt.kc + dt.kx*p.X + dt.ky*p.Y + dt.kz*p.Z
}

type AnimatedTransform struct {
	startTransform, endTransform *Transform
	startTime, endTime           float64
	actuallyAnimated             bool
	startT, endT                 *Vector3f
	startR, endR                 *Quaternion
	startS, endS                 *Matrix4x4
	hasRotation                  bool
	c1, c2, c3, c4, c5           [3]*DerivativeTerm
}

func (at *AnimatedTransform) Interpolate(time float64) *Transform {
	if !at.actuallyAnimated || time <= at.startTime {
		return at.startTransform
	}
	if time >= at.endTime {
		return at.endTransform
	}

	dt := (time - at.startTime) / (at.endTime - at.startTime)
	// interpolate translation
	trans := at.startT.MulScalar(1.0 - dt).Add(at.endT.MulScalar(dt))

	// interpolate rotation
	rotate := Slerp(dt, at.startR, at.endR)

	// interpolate scale
	var scale *Matrix4x4
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			scale[i][j] = Lerp(dt, at.startS[i][j], at.endS[i][j])
		}
	}

	// compute interpolated matrix as product of interpolated components
	return Translate(trans).Mul(rotate.ToTransform()).Mul(NewTransform(scale))

}
