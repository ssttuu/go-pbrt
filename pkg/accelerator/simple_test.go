package accelerator

import (
	"testing"

	"github.com/golang/mock/gomock"
	"github.com/stretchr/testify/assert"
	"github.com/stupschwartz/go-pbrt/pkg/pbrt"
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
	s := Simple{primitives: primitives}
	radius := &pbrt.Point3f{1,1,1}

	//
	//ray := pbrt.NewRay(&pbrt.Point3f{0, 0, 0}, &pbrt.Vector3f{0, 0, 1.0}, 0)
	//si := pbrt.NewSurfaceInteraction()
	//intersects := s.Intersect(ray, si)
	//assert.True(t, intersects)
	//assert.Equal(t, &pbrt.Point3f{0.0, 0.0, 4.0}, si.Point)
	//
	////
	//ray = pbrt.NewRay(&pbrt.Point3f{0, 0, 0}, &pbrt.Vector3f{1, 1, 1}, 0)
	//si = pbrt.NewSurfaceInteraction()
	//intersects = s.Intersect(ray, si)
	//assert.True(t, intersects)
	//
	//expected := &pbrt.Point3f{10, 10, 10}
	//assert.Equal(t, expected.Sub(radius.Normalized()), si.Point)

	//
	dir := &pbrt.Vector3f{-1, -1, -1}
	ray := pbrt.NewRay(&pbrt.Point3f{15, 15, 15}, dir.Normalized(), 0)
	si := pbrt.NewSurfaceInteraction()
	intersects := s.Intersect(ray, si)
	assert.True(t, intersects)

	expected := &pbrt.Point3f{10, 10, 10}
	assert.Equal(t, expected.Add(radius.Normalized()), si.Point)

}

func TestSimple_Intersect_ReturnsFalseWhenNoIntersections(t *testing.T) {
	ctrl := gomock.NewController(t)
	defer ctrl.Finish()

	prim := pbrt.NewMockPrimitive(ctrl)

	ray := pbrt.NewRay(&pbrt.Point3f{0, 0, 0}, &pbrt.Vector3f{0, 0, 1.0}, 0)
	si := pbrt.NewSurfaceInteraction()

	prim.
		EXPECT().
		Intersect(ray, gomock.Any()).
		Times(1).
		Return(false)

	s := Simple{primitives: []pbrt.Primitive{prim}}

	intersects := s.Intersect(ray, si)
	assert.False(t, intersects)
}

func TestSimple_IntersectP(t *testing.T) {
	ctrl := gomock.NewController(t)
	defer ctrl.Finish()

	prim := pbrt.NewMockPrimitive(ctrl)
	ray := pbrt.NewRay(&pbrt.Point3f{0, 0, 0}, &pbrt.Vector3f{0, 0, 1.0}, 0)

	prim.
		EXPECT().
		IntersectP(ray).
		Times(1).
		Return(false)

	s := Simple{primitives: []pbrt.Primitive{prim}}

	assert.False(t, s.IntersectP(ray))

}
