package accelerator

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

var primitives = []pbrt.Primitive{
	pbrt.NewGeometricPrimitive(
		pbrt.NewSphereShape(
			"prim1",
			pbrt.Translate(&pbrt.Vector3f{0, 0, 5}),
			false,
			1.0,
		),
		nil,
	),
	pbrt.NewGeometricPrimitive(
		pbrt.NewSphereShape(
			"prim2",
			pbrt.Translate(&pbrt.Vector3f{0, 0, 10}),
			false,
			1.0,
		),
		nil,
	),
	pbrt.NewGeometricPrimitive(
		pbrt.NewSphereShape(
			"prim3",
			pbrt.Translate(&pbrt.Vector3f{10, 10, 10}),
			true,
			1.0,
		),
		nil,
	),
}

func TestSimple_Intersect(t *testing.T) {
	direction := &pbrt.Vector3f{-1, -1, -1}
	direction.Normalize()
	ray := pbrt.NewRay(&pbrt.Point3f{15, 15, 15}, direction, 0)
	si := pbrt.NewSurfaceInteraction()

	s := Simple{primitives: primitives}

	radius := &pbrt.Point3f{1, 1, 1}
	radius.Normalize()
	expected := &pbrt.Point3f{10, 10, 10}
	expected.AddAssign(radius)

	intersects := s.Intersect(ray, si)
	assert.True(t, intersects)
	assert.Equal(t, expected, si.Point)

}

func TestSimple_Intersect_ReturnsFalseWhenNoIntersections(t *testing.T) {
	ray := pbrt.NewRay(&pbrt.Point3f{0, 0, 0}, &pbrt.Vector3f{0, 0, -1.0}, 0)
	si := pbrt.NewSurfaceInteraction()

	s := Simple{primitives: primitives}

	intersects := s.Intersect(ray, si)
	assert.False(t, intersects)
}

func TestSimple_IntersectP(t *testing.T) {
	tests := []struct {
		ray      *pbrt.Ray
		expected bool
	}{
		{
			ray: pbrt.NewRay(
				&pbrt.Point3f{0, 0, 0},
				&pbrt.Vector3f{0, 0, 1.0},
				0),
			expected: true,
		},
		{
			ray: pbrt.NewRay(
				&pbrt.Point3f{0, 0, 0},
				&pbrt.Vector3f{0, 0, -1.0},
				0),
			expected: false,
		},
		{
			ray: pbrt.NewRay(
				&pbrt.Point3f{0, 0, 500},
				&pbrt.Vector3f{0, 0, -1.0},
				0),
			expected: true,
		},
		{
			ray: pbrt.NewRay(
				&pbrt.Point3f{10, 10, 500},
				&pbrt.Vector3f{0, 0, -1.0},
				0),
			expected: true,
		},
	}

	s := Simple{primitives: primitives}
	for _, test := range tests {
		assert.Equal(t, test.expected, s.IntersectP(test.ray))
	}
}
