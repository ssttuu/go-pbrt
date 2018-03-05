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
			r[i][j] = m[i][0] * other[0][j] + m[i][1] * other[1][j] + m[i][2] * other[2][j] + m[i][3] * m[3][j]
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
		matrix: m,
		matrixInverse: mInv,
	}
}

func (t *Transform) Mul(other *Transform) *Transform {
	return &Transform{
		matrix: t.matrix.Mul(other.matrix),
		matrixInverse: t.matrixInverse.Mul(other.matrixInverse),
	}
}

func Translate(delta *Vector3) *Transform {
	return &Transform{
		matrix:&Matrix4x4{
			{1, 0, 0, delta.x},
			{0, 1, 0, delta.y},
			{0, 0, 1, delta.z},
			{0, 0, 0, 1},
		},
		matrixInverse:&Matrix4x4{
			{1, 0, 0, -delta.x},
			{0, 1, 0, -delta.y},
			{0, 0, 1, -delta.z},
			{0, 0, 0, 1},
		},
	}
}

func Scale(x, y, z float64) *Transform {
	return &Transform{
		matrix:&Matrix4x4{
			{x, 0, 0, 0},
			{0, y, 0, 0},
			{0, 0, z, 0},
			{0, 0, 0, 1},
		},
		matrixInverse:&Matrix4x4{
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
		matrix:m,
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
		matrix:m,
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
		matrix:m,
		matrixInverse: m.Transpose(),
	}
}

func Rotate(theta float64, axis Vector3) *Transform {
	axis.Normalize()
	sinTheta := math.Sin(Radians(theta))
	cosTheta := math.Cos(Radians(theta))
	var m *Matrix4x4

	// Compute rotation of first basis vector
	m[0][0] = axis.x * axis.x + (1 - axis.x * axis.x) * cosTheta
	m[0][1] = axis.x * axis.y * (1 - cosTheta) - axis.z * sinTheta
	m[0][2] = axis.x * axis.z * (1 - cosTheta) + axis.y * sinTheta
	m[0][3] = 0

	// Compute rotations of second and third basis vectors
	m[1][0] = axis.x * axis.y * (1 - cosTheta) + axis.z * sinTheta
	m[1][1] = axis.y * axis.y + (1 - axis.y * axis.y) * cosTheta
	m[1][2] = axis.y * axis.z * (1 - cosTheta) - axis.x * sinTheta
	m[1][3] = 0

	m[2][0] = axis.x * axis.z * (1 - cosTheta) - axis.y * sinTheta
	m[2][1] = axis.y * axis.z * (1 - cosTheta) + axis.x * sinTheta
	m[2][2] = axis.z * axis.z + (1 - axis.z * axis.z) * cosTheta
	m[2][3] = 0

	return &Transform{matrix:m, matrixInverse: m.Transpose()}
}

func LookAt(pos *Point3, look *Point3, up *Vector3) (*Transform, error) {
	var cameraToWorld *Matrix4x4
	// Initialize fourth column of viewing matrix
	cameraToWorld[0][3] = pos.x
	cameraToWorld[1][3] = pos.y
	cameraToWorld[2][3] = pos.z
	cameraToWorld[3][3] = 1

	// Initialize first three columns of viewing matrix
	dir := look.Sub(pos).Normalized()
	if up.Normalized().Cross(dir).Length() == 0 {
		return nil, ErrSameDirectionVectors
	}
	left := up.Normalized().Cross(dir).Normalized()
	newUp := dir.Cross(left)
	cameraToWorld[0][0] = left.x
	cameraToWorld[1][0] = left.y
	cameraToWorld[2][0] = left.z
	cameraToWorld[3][0] = 0.
	cameraToWorld[0][1] = newUp.x
	cameraToWorld[1][1] = newUp.y
	cameraToWorld[2][1] = newUp.z
	cameraToWorld[3][1] = 0.
	cameraToWorld[0][2] = dir.x
	cameraToWorld[1][2] = dir.y
	cameraToWorld[2][2] = dir.z
	cameraToWorld[3][2] = 0.

	inv, err := cameraToWorld.Inverse()
	if err != nil {
		return nil, err
	}
	return &Transform{matrix:inv, matrixInverse:cameraToWorld}, nil
}

func Orthographic(zNear, zFar float64) *Transform {
	return Scale(1, 1, 1 / (zFar - zNear)).Mul(Translate(&Vector3{0, 0, -zNear}))
}

func Perspective(fov, n, f float64) *Transform {
	persp := Matrix4x4{
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, f / (f - n), -f * n / (f - n)},
		{0, 0, 1, 0},
	}

	invTanAng := 1.0 / math.Tan(Radians(fov) / 2)
	return Scale(invTanAng, invTanAng, 1).Mul(NewTransform(&persp))
}

