package accelerator

import (
	"testing"
	"github.com/golang/mock/gomock"
	"github.com/stupschwartz/go-pbrt/pkg/pbrt"
	"github.com/stretchr/testify/assert"
)

func TestSimple_Intersect(t *testing.T) {
	ctrl := gomock.NewController(t)
	defer ctrl.Finish()

	prim1 := pbrt.NewMockPrimitive(ctrl)
	prim2 := pbrt.NewMockPrimitive(ctrl)

	ray := pbrt.NewRay(&pbrt.Point3f{0, 0, 0}, &pbrt.Vector3f{0, 0, 1.0}, 0)
	si := pbrt.NewSurfaceInteraction()

	prim1.
		EXPECT().
		Intersect(ray, gomock.Any()).
		Times(1).
		DoAndReturn(func(r *pbrt.Ray, si *pbrt.SurfaceInteraction) bool {
			si.Point.Set(0.0, 0.0, 5.0)
			return true
		})

	prim2.
		EXPECT().
		Intersect(ray, gomock.Any()).
		Times(1).
		DoAndReturn(func(r *pbrt.Ray, si *pbrt.SurfaceInteraction) bool {
		si.Point.Set(0.0, 0.0, 4.0)
		return true
	})

	s := Simple{primitives:[]pbrt.Primitive{prim1, prim2}}

	intersects := s.Intersect(ray, si)

	assert.True(t, intersects)
	assert.Equal(t, &pbrt.Point3f{0.0, 0.0, 4.0}, si.Point)
}

func TestSimple_Intersect_ReturnsFalseWhenNoIntersections(t *testing.T)  {
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

	s := Simple{primitives:[]pbrt.Primitive{prim}}

	intersects := s.Intersect(ray, si)
	assert.False(t, intersects)
}
