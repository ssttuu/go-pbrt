package pbrt_test

import (
	"testing"

	"github.com/golang/mock/gomock"
	"github.com/ssttuu/go-pbrt/pkg/pbrt"
	"github.com/stretchr/testify/assert"
)

func TestScene_Intersect(t *testing.T) {
	ctrl := gomock.NewController(t)
	defer ctrl.Finish()

	aggregate := pbrt.NewMockAggregate(ctrl)
	aggregate.
		EXPECT().
		WorldBound().
		Times(1).
		Return(&pbrt.Bounds3{})

	scene := pbrt.NewScene(aggregate, []pbrt.Light{}, []pbrt.Light{})

	ray := pbrt.NewRay(&pbrt.Point3f{}, &pbrt.Vector3f{}, 0)
	si := pbrt.NewSurfaceInteraction()

	aggregate.
		EXPECT().
		Intersect(ray, si).
		Times(1).
		Return(true)

	intersects := scene.Intersect(ray, si)
	assert.True(t, intersects)
}

func TestScene_IntersectP(t *testing.T) {
	ctrl := gomock.NewController(t)
	defer ctrl.Finish()

	aggregate := pbrt.NewMockAggregate(ctrl)

	aggregate.
		EXPECT().
		WorldBound().
		Times(1).
		Return(&pbrt.Bounds3{})

	scene := pbrt.NewScene(aggregate, []pbrt.Light{}, []pbrt.Light{})

	ray := pbrt.NewRay(&pbrt.Point3f{}, &pbrt.Vector3f{}, 0)

	aggregate.
		EXPECT().
		IntersectP(ray).
		Times(1).
		Return(true)

	intersects := scene.IntersectP(ray)
	assert.True(t, intersects)
}

func TestScene_IntersectTr(t *testing.T) {
	ctrl := gomock.NewController(t)
	defer ctrl.Finish()

	aggregate := pbrt.NewMockAggregate(ctrl)
	primitive := pbrt.NewMockPrimitive(ctrl)
	material := pbrt.NewMockMaterial(ctrl)

	aggregate.
		EXPECT().
		WorldBound().
		Times(1).
		Return(&pbrt.Bounds3{})

	scene := pbrt.NewScene(aggregate, []pbrt.Light{}, []pbrt.Light{})

	ray := pbrt.NewRay(&pbrt.Point3f{}, &pbrt.Vector3f{}, 0)
	si := pbrt.NewSurfaceInteraction()
	si.Primitive = primitive
	sampler := pbrt.NewMockSampler(ctrl)

	aggregate.
		EXPECT().
		Intersect(ray, si).
		Times(1).
		Return(true)

	primitive.
		EXPECT().
		GetMaterial().
		Times(1).
		Return(material)

	transmittance := pbrt.NewSpectrum(0.1)
	intersects := scene.IntersectTr(ray, si, sampler, transmittance)
	assert.True(t, intersects)
	assert.Equal(t, pbrt.NewSpectrum(1.0), transmittance)
}
