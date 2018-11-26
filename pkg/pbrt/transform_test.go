package pbrt

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestMatrix4x4_Equals(t *testing.T) {
	m := NewMatrix4x4()
	assert.True(t, m.Equals(NewMatrix4x4()))

	m[3][3] = 10.0
	assert.False(t, m.Equals(NewMatrix4x4()))
}

func TestMatrix4x4_Inverse(t *testing.T) {
	m := &Matrix4x4{
		{1, 0, 0, 0},
		{0, 1, 1, 0},
		{0, 0, 1, 0},
		{0, 0, 2, 1},
	}

	expected := &Matrix4x4{
		{1, 0, 0, 0},
		{0, 1, -1, 0},
		{0, 0, 1, 0},
		{0, 0, -2, 1},
	}

	inverse, err := m.Inverse()
	assert.Nil(t, err)
	assert.Equal(t, expected, inverse)

}

func TestMatrix4x4_Mul(t *testing.T) {

}

func TestMatrix4x4_NotEquals(t *testing.T) {

}

func TestMatrix4x4_Transpose(t *testing.T) {

}

func TestTransform_Inverse(t *testing.T) {

}

func TestTransform_IsIdentity(t *testing.T) {

}

func TestTransform_Mul(t *testing.T) {

}

func TestTransform_TransformBounds(t *testing.T) {

}

func TestTransform_TransformPoint(t *testing.T) {
	xform := NewTransform(NewMatrix4x4())
	p, pError := xform.TransformPoint(&Point3f{}, new(Vector3f))
	assert.Equal(t, &Point3f{}, p)
	assert.Equal(t, &Point3f{}, pError)

	xform = Translate(&Vector3f{5, 4, 3})
	p, pError = xform.TransformPoint(&Point3f{}, new(Vector3f))
	assert.Equal(t, &Point3f{5, 4, 3}, p)
}

func TestTransform_TransformRay(t *testing.T) {
	xform := RotateY(90).Mul(Translate(&Vector3f{5, 4, 3}))
	ray, _, _ := xform.TransformRay(NewRay(&Point3f{}, &Vector3f{1, 0, 0}, 0))
	assert.Equal(t, NewRay(&Point3f{3.0000000000000004, 4, -5}, &Vector3f{6.123233995736757e-17, 0, -1}, 0), ray)
}
