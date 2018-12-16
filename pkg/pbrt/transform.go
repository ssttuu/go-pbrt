package pbrt

import (
	"log"

	"github.com/pkg/errors"
	"github.com/ssttuu/go-pbrt/pkg/math"
)

var ErrSingularMatrix = errors.New("singular Matrix in Matrix invert")
var ErrSameDirectionVectors = errors.New("same direction in Matrix invert")

func SolveLinearSystem2x2(A [2][2]float64, B [2]float64) (solvable bool, x0, x1 float64) {
	det := A[0][0]*A[1][1] - A[0][1]*A[1][0]
	if math.Abs(det) < 1e-10 {
		return false, 0, 0
	}
	x0 = (A[1][1]*B[0] - A[0][1]*B[1]) / det
	x1 = (A[0][0]*B[1] - A[1][0]*B[0]) / det
	if math.IsNaN(x0) || math.IsNaN(x1) {
		return false, 0, 0
	}
	return true, x0, x1

}

type Matrix4x4 [4][4]float64

func NewMatrix4x4() *Matrix4x4 {
	return &Matrix4x4{
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
	return !m.Equals(other)
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
	var r Matrix4x4
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			r[i][j] = m[i][0]*other[0][j] + m[i][1]*other[1][j] + m[i][2]*other[2][j] + m[i][3]*m[3][j]
		}
	}
	return &r
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
							big = math.Abs(minv[j][k])
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
	Matrix, MatrixInverse *Matrix4x4
}

func NewTransform(m *Matrix4x4) *Transform {
	mInv, err := m.Inverse()
	if err != nil {
		log.Fatal("TODO")
	}
	return &Transform{
		Matrix:        m,
		MatrixInverse: mInv,
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
	return t.Matrix[0][0] == 1.0 && t.Matrix[0][1] == 0.0 && t.Matrix[0][2] == 0.0 && t.Matrix[0][3] == 0.0 &&
		t.Matrix[1][0] == 0.0 && t.Matrix[1][1] == 1.0 && t.Matrix[1][2] == 0.0 && t.Matrix[1][3] == 0.0 &&
		t.Matrix[2][0] == 0.0 && t.Matrix[2][1] == 0.0 && t.Matrix[2][2] == 1.0 && t.Matrix[2][3] == 0.0 &&
		t.Matrix[3][0] == 0.0 && t.Matrix[3][1] == 0.0 && t.Matrix[3][2] == 0.0 && t.Matrix[3][3] == 1.0

}

func (t *Transform) Inverse() *Transform {
	return &Transform{t.MatrixInverse, t.Matrix}
}

func (t *Transform) Mul(other *Transform) *Transform {
	return &Transform{
		Matrix:        t.Matrix.Mul(other.Matrix),
		MatrixInverse: t.MatrixInverse.Mul(other.MatrixInverse),
	}
}

//func (t *Transform) TransformPoint(p *Point3f) *Point3f {
//	// compute transformed points from p
//	xp := t.Matrix[0][0]*p.X + t.Matrix[0][1]*p.Y + t.Matrix[0][2]*p.Z + t.Matrix[0][3]
//	yp := t.Matrix[1][0]*p.X + t.Matrix[1][1]*p.Y + t.Matrix[1][2]*p.Z + t.Matrix[1][3]
//	zp := t.Matrix[2][0]*p.X + t.Matrix[2][1]*p.Y + t.Matrix[2][2]*p.Z + t.Matrix[2][3]
//	wp := t.Matrix[3][0]*p.X + t.Matrix[3][1]*p.Y + t.Matrix[3][2]*p.Z + t.Matrix[3][3]
//
//	newPoint := &Point3f{xp, yp, zp}
//
//	if wp == 1.0 {
//		return newPoint
//	}
//
//	return newPoint.DivScalar(wp)
//}
//
//
//func (t *Transform) TransformPointAndError(p *Point3f) (newPoint *Point3f, newPointError *Vector3f) {
//	// compute transformed points from p
//	xp := t.Matrix[0][0]*p.X + t.Matrix[0][1]*p.Y + t.Matrix[0][2]*p.Z + t.Matrix[0][3]
//	yp := t.Matrix[1][0]*p.X + t.Matrix[1][1]*p.Y + t.Matrix[1][2]*p.Z + t.Matrix[1][3]
//	zp := t.Matrix[2][0]*p.X + t.Matrix[2][1]*p.Y + t.Matrix[2][2]*p.Z + t.Matrix[2][3]
//	wp := t.Matrix[3][0]*p.X + t.Matrix[3][1]*p.Y + t.Matrix[3][2]*p.Z + t.Matrix[3][3]
//
//	// compute absolute error for transformed point
//	*newPointError = Vector3f{
//		X: math.Abs(t.Matrix[0][0] * p.X) + math.Abs(t.Matrix[0][1] * p.Y) + math.Abs(t.Matrix[0][2] * p.Z) + math.Abs(t.Matrix[0][3]),
//		Y: math.Abs(t.Matrix[1][0] * p.X) + math.Abs(t.Matrix[1][1] * p.Y) + math.Abs(t.Matrix[1][2] * p.Z) + math.Abs(t.Matrix[0][3]),
//		Z: math.Abs(t.Matrix[2][0] * p.X) + math.Abs(t.Matrix[2][1] * p.Y) + math.Abs(t.Matrix[2][2] * p.Z) + math.Abs(t.Matrix[0][3]),
//	}
//	newPointError = newPointError.MulScalar(Gamma(3))
//
//	*newPoint = Point3f{xp, yp, zp}
//
//	if wp == 1.0 {
//		return newPoint, newPointError
//	}
//
//	return newPoint.DivScalar(wp), newPointError
//}

func (t *Transform) TransformPoint(p *Point3f, pError *Vector3f) (*Point3f, *Vector3f) {
	// compute transformed points from p
	xp := t.Matrix[0][0]*p.X + t.Matrix[0][1]*p.Y + t.Matrix[0][2]*p.Z + t.Matrix[0][3]
	yp := t.Matrix[1][0]*p.X + t.Matrix[1][1]*p.Y + t.Matrix[1][2]*p.Z + t.Matrix[1][3]
	zp := t.Matrix[2][0]*p.X + t.Matrix[2][1]*p.Y + t.Matrix[2][2]*p.Z + t.Matrix[2][3]
	wp := t.Matrix[3][0]*p.X + t.Matrix[3][1]*p.Y + t.Matrix[3][2]*p.Z + t.Matrix[3][3]

	absError := &Vector3f{
		X: (math.Gamma(3.0)+1.0)*(math.Abs(t.Matrix[0][0])*pError.X+math.Abs(t.Matrix[0][1])*pError.Y+math.Abs(t.Matrix[0][2])*pError.Z) + (math.Gamma(3) * (math.Abs(t.Matrix[0][0]*p.X) + math.Abs(t.Matrix[0][1])*p.Y + math.Abs(t.Matrix[0][2]*p.Z+math.Abs(t.Matrix[0][3])))),
		Y: (math.Gamma(3.0)+1.0)*(math.Abs(t.Matrix[1][0])*pError.X+math.Abs(t.Matrix[1][1])*pError.Y+math.Abs(t.Matrix[1][2])*pError.Z) + (math.Gamma(3) * (math.Abs(t.Matrix[1][0]*p.X) + math.Abs(t.Matrix[1][1])*p.Y + math.Abs(t.Matrix[1][2]*p.Z+math.Abs(t.Matrix[1][3])))),
		Z: (math.Gamma(3.0)+1.0)*(math.Abs(t.Matrix[2][0])*pError.X+math.Abs(t.Matrix[2][1])*pError.Y+math.Abs(t.Matrix[2][2])*pError.Z) + (math.Gamma(3) * (math.Abs(t.Matrix[2][0]*p.X) + math.Abs(t.Matrix[2][1])*p.Y + math.Abs(t.Matrix[2][2]*p.Z+math.Abs(t.Matrix[2][3])))),
	}

	newPoint := &Point3f{xp, yp, zp}

	if wp == 1.0 {
		return newPoint, absError
	}

	return newPoint.DivScalar(wp), absError
}

func (t *Transform) TransformVector(v *Vector3f) *Vector3f {
	return &Vector3f{
		t.Matrix[0][0]*v.X + t.Matrix[0][1]*v.Y + t.Matrix[0][2]*v.Z,
		t.Matrix[1][0]*v.X + t.Matrix[1][1]*v.Y + t.Matrix[1][2]*v.Z,
		t.Matrix[2][0]*v.X + t.Matrix[2][1]*v.Y + t.Matrix[2][2]*v.Z,
	}
}

func (t *Transform) TransformVectorWithAbsError(v *Vector3f) (newVector *Vector3f, absError *Vector3f) {
	newVector = &Vector3f{
		t.Matrix[0][0]*v.X + t.Matrix[0][1]*v.Y + t.Matrix[0][2]*v.Z,
		t.Matrix[1][0]*v.X + t.Matrix[1][1]*v.Y + t.Matrix[1][2]*v.Z,
		t.Matrix[2][0]*v.X + t.Matrix[2][1]*v.Y + t.Matrix[2][2]*v.Z,
	}
	absError = &Vector3f{
		X: math.Gamma(3) * (math.Abs(t.Matrix[0][0]*v.X) + math.Abs(t.Matrix[0][1]*v.Y) + math.Abs(t.Matrix[0][2]*v.Z)),
		Y: math.Gamma(3) * (math.Abs(t.Matrix[1][0]*v.X) + math.Abs(t.Matrix[1][1]*v.Y) + math.Abs(t.Matrix[1][2]*v.Z)),
		Z: math.Gamma(3) * (math.Abs(t.Matrix[2][0]*v.X) + math.Abs(t.Matrix[2][1]*v.Y) + math.Abs(t.Matrix[2][2]*v.Z)),
	}
	return newVector, absError
}

func (t *Transform) TransformNormal(n *Vector3f) *Vector3f {
	return &Vector3f{
		X: t.MatrixInverse[0][0]*n.X + t.MatrixInverse[1][0]*n.Y + t.MatrixInverse[2][0]*n.Z,
		Y: t.MatrixInverse[0][1]*n.X + t.MatrixInverse[1][1]*n.Y + t.MatrixInverse[2][1]*n.Z,
		Z: t.MatrixInverse[0][2]*n.X + t.MatrixInverse[1][2]*n.Y + t.MatrixInverse[2][2]*n.Z,
	}
}

func (t *Transform) TransformRay(ray *Ray) (r *Ray, originErr, directionErr *Vector3f) {
	origin, originErr := t.TransformPoint(ray.Origin, new(Vector3f))
	direction, directionErr := t.TransformVectorWithAbsError(ray.Direction)
	lengthSquared := direction.LengthSquared()
	if lengthSquared > 0 {
		dt := direction.Abs().Dot(originErr) / lengthSquared
		origin.AddAssign(direction.MulScalar(dt))
	}
	return &Ray{origin, direction, ray.TMax, ray.Time, ray.Medium}, originErr, directionErr
}

func (t *Transform) TransformSurfaceInteraction(si *SurfaceInteraction) *SurfaceInteraction {
	ret := *si

	ret.interaction = si.interaction
	ret.Point, ret.PointError = t.TransformPoint(si.Point, si.PointError)

	ret.Normal = t.TransformNormal(si.Normal).Normalized()
	ret.Wo = t.TransformVector(si.Wo).Normalized()
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
	//ret.dpdx = Type.TransformVector(si.dpdx)
	//ret.dpdy = Type.TransformVector(si.dpdy)
	ret.BSDF = si.BSDF
	ret.Primitive = si.Primitive
	ret.shading.normal = FaceForward(ret.shading.normal, ret.Normal)
	ret.faceIndex = si.faceIndex
	return &ret
}

func (t *Transform) TransformBounds(b *Bounds3) *Bounds3 {
	corner, _ := t.TransformPoint(b.Min, new(Vector3f))
	ret := &Bounds3{Min: corner, Max: corner}

	for i := 1; i < 8; i++ {
		corner, _ = t.TransformPoint(b.Corner(i), new(Vector3f))
		ret.UnionPoint(corner)
	}
	return ret
}

func Translate(delta *Vector3f) *Transform {
	return &Transform{
		Matrix: &Matrix4x4{
			{1, 0, 0, delta.X},
			{0, 1, 0, delta.Y},
			{0, 0, 1, delta.Z},
			{0, 0, 0, 1},
		},
		MatrixInverse: &Matrix4x4{
			{1, 0, 0, -delta.X},
			{0, 1, 0, -delta.Y},
			{0, 0, 1, -delta.Z},
			{0, 0, 0, 1},
		},
	}
}

func Scale(x, y, z float64) *Transform {
	return &Transform{
		Matrix: &Matrix4x4{
			{x, 0, 0, 0},
			{0, y, 0, 0},
			{0, 0, z, 0},
			{0, 0, 0, 1},
		},
		MatrixInverse: &Matrix4x4{
			{1.0 / x, 0, 0, 0},
			{0, 1.0 / y, 0, 0},
			{0, 0, 1.0 / z, 0},
			{0, 0, 0, 1},
		},
	}
}

func RotateX(degrees float64) *Transform {
	sinTheta := math.Sin(math.Radians(degrees))
	cosTheta := math.Cos(math.Radians(degrees))
	m := &Matrix4x4{
		{1, 0, 0, 0},
		{0, cosTheta, -sinTheta, 0},
		{0, sinTheta, cosTheta, 0},
		{0, 0, 0, 1},
	}
	return &Transform{
		Matrix:        m,
		MatrixInverse: m.Transpose(),
	}
}

func RotateY(degrees float64) *Transform {
	sinTheta := math.Sin(math.Radians(degrees))
	cosTheta := math.Cos(math.Radians(degrees))
	m := &Matrix4x4{
		{cosTheta, 0, sinTheta, 0},
		{0, 1, 0, 0},
		{-sinTheta, 0, cosTheta, 0},
		{0, 0, 0, 1},
	}
	return &Transform{
		Matrix:        m,
		MatrixInverse: m.Transpose(),
	}
}

func RotateZ(degrees float64) *Transform {
	sinTheta := math.Sin(math.Radians(degrees))
	cosTheta := math.Cos(math.Radians(degrees))
	m := &Matrix4x4{
		{cosTheta, -sinTheta, 0, 0},
		{sinTheta, cosTheta, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
	}
	return &Transform{
		Matrix:        m,
		MatrixInverse: m.Transpose(),
	}
}

func Rotate(theta float64, axis Vector3f) *Transform {
	a := axis.Normalized()
	sinTheta := math.Sin(math.Radians(theta))
	cosTheta := math.Cos(math.Radians(theta))

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

	return &Transform{Matrix: m, MatrixInverse: m.Transpose()}
}

func LookAt(pos *Point3f, look *Point3f, up *Vector3f) (*Transform, error) {
	cameraToWorld := new(Matrix4x4)
	// Initialize fourth column of viewing Matrix
	cameraToWorld[0][3] = pos.X
	cameraToWorld[1][3] = pos.Y
	cameraToWorld[2][3] = pos.Z
	cameraToWorld[3][3] = 1

	// Initialize first three columns of viewing Matrix
	dir := look.Sub(pos).Normalized()
	if up.Normalized().Cross(dir).Length() == 0 {
		return nil, ErrSameDirectionVectors
	}
	right := up.Normalized().Cross(dir).Normalized()
	newUp := dir.Cross(right)
	cameraToWorld[0][0] = right.X
	cameraToWorld[1][0] = right.Y
	cameraToWorld[2][0] = right.Z
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
		return nil, errors.Wrap(err, "inverting camera")
	}
	return &Transform{Matrix: cameraToWorld, MatrixInverse: inv}, nil
}

func Orthographic(zNear, zFar float64) *Transform {
	return Scale(1, 1, 1/(zFar-zNear)).Mul(Translate(&Vector3f{0, 0, -zNear}))
}

func Perspective(fov, n, f float64) *Transform {
	persp := &Matrix4x4{
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, f / (f - n), -f * n / (f - n)},
		{0, 0, 1, 0},
	}

	invTanAng := 1.0 / math.Tan(math.Radians(fov)/2)
	return Scale(invTanAng, invTanAng, 1).Mul(NewTransform(persp))
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

func NewAnimatedTransform(start, end *Transform, startTime, endTime float64) *AnimatedTransform {
	actuallyAnimated := start != end
	at := &AnimatedTransform{
		startTransform:   start,
		endTransform:     end,
		startTime:        startTime,
		endTime:          endTime,
		actuallyAnimated: actuallyAnimated,
	}

	if !actuallyAnimated {
		return at
	}

	// TODO
	//at.startT, at.startR, at.startS = Decompose(start.Matrix)
	//at.endT, at.endR, at.endS = Decompose(end.Matrix)

	// Flip r if needed to select shortest path
	if at.startR.Dot(at.endR) < 0 {
		at.endR = at.endR.MulScalar(-1)
	}

	at.hasRotation = at.startR.Dot(at.endR) < 0.9995

	// compute terms of motion derivative function
	if at.hasRotation {
		// TODO
	}

	return at
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
			scale[i][j] = math.Lerp(dt, at.startS[i][j], at.endS[i][j])
		}
	}

	// compute interpolated Matrix as product of interpolated components
	return Translate(trans).Mul(rotate.ToTransform()).Mul(NewTransform(scale))
}

func (at *AnimatedTransform) MotionBounds(b *Bounds3) *Bounds3 {
	if !at.actuallyAnimated {
		return at.startTransform.TransformBounds(b)
	}

	// TODO:
	return nil
}

func (at *AnimatedTransform) TransformRay(r *Ray) *Ray {
	if !at.actuallyAnimated || r.Time <= at.startTime {
		ray, _, _ := at.startTransform.TransformRay(r)
		return ray
	} else if r.Time >= at.endTime {
		ray, _, _ := at.endTransform.TransformRay(r)
		return ray
	}

	ray, _, _ := at.Interpolate(r.Time).TransformRay(r)
	return ray
}

func (at *AnimatedTransform) TransformPointAtTime(p *Point3f, time float64) *Point3f {
	pError := new(Vector3f)
	if !at.actuallyAnimated || time <= at.startTime {
		pt, _ := at.startTransform.TransformPoint(p, pError)
		return pt
	} else if time >= at.endTime {
		pt, _ := at.endTransform.TransformPoint(p, pError)
		return pt
	}

	pt, _ := at.Interpolate(time).TransformPoint(p, pError)
	return pt
}

func (at *AnimatedTransform) TransformVectorAtTime(v *Vector3f, time float64) *Vector3f {
	pError := new(Vector3f)
	if !at.actuallyAnimated || time <= at.startTime {
		pt, _ := at.startTransform.TransformPoint(v, pError)
		return pt
	} else if time >= at.endTime {
		pt, _ := at.endTransform.TransformPoint(v, pError)
		return pt
	}

	pt, _ := at.Interpolate(time).TransformPoint(v, pError)
	return pt
}
