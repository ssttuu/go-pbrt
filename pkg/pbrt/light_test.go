package pbrt

import (
	"testing"

	"github.com/golang/mock/gomock"
	"github.com/stretchr/testify/assert"
)

func TestVisibilityTester_Unoccluded(t *testing.T) {
	ctrl := gomock.NewController(t)
	defer ctrl.Finish()

	mediumAccessor := &MediumAccessor{Inside: nil, Outside: nil}
	p0 := NewInteraction(&Point3f{}, &Point3f{}, &Normal3f{}, nil, 0, mediumAccessor)
	p1 := NewInteraction(&Point3f{10, 0, 0}, &Point3f{}, &Normal3f{}, nil, 0, mediumAccessor)
	vt := NewVisibilityTester(p0, p1)

	scene := NewMockScene(ctrl)

	ray := &Ray{
		Origin:    &Point3f{0, 0, 0},
		Direction: &Vector3f{10, 0, 0},
		TMax:      0.9999,
		Time:      0.0,
		Medium:    nil,
	}

	scene.
		EXPECT().
		IntersectP(ray).
		Times(1).
		Return(true)

	assert.False(t, vt.Unoccluded(scene))

	scene.
		EXPECT().
		IntersectP(ray).
		Times(1).
		Return(false)

	assert.True(t, vt.Unoccluded(scene))
}
