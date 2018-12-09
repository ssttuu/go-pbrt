package accelerator

import (
	"fmt"
	"testing"

	"github.com/stupschwartz/go-pbrt/pkg/pbrt"

	"github.com/stretchr/testify/assert"
)

func Swap(v []*MortonPrimitive) {
	t := make([]*MortonPrimitive, len(v))
	for i := 0; i < len(v); i++ {
		t[i] = &MortonPrimitive{5, 20}
	}
	v = t
}

func TestRadixSortInPlace(t *testing.T) {
	var mps []*MortonPrimitive
	for i := 0; i < 10; i++ {
		mps = append(mps, &MortonPrimitive{i, uint64(10 - (i + 1))})
	}

	RadixSortInPlace(&mps)

	expected := []*MortonPrimitive{
		{primitiveIndex: 9, mortonCode: 0},
		{primitiveIndex: 8, mortonCode: 1},
		{primitiveIndex: 7, mortonCode: 2},
		{primitiveIndex: 6, mortonCode: 3},
		{primitiveIndex: 5, mortonCode: 4},
		{primitiveIndex: 4, mortonCode: 5},
		{primitiveIndex: 3, mortonCode: 6},
		{primitiveIndex: 2, mortonCode: 7},
		{primitiveIndex: 1, mortonCode: 8},
		{primitiveIndex: 0, mortonCode: 9},
	}
	assert.Equal(t, expected, mps)
}

func TestBVH_Intersect(t *testing.T) {
	tests := []struct {
		ray                *pbrt.Ray
		intersects         bool
		surfaceInteraction *pbrt.SurfaceInteraction
	}{
		{
			ray: pbrt.NewRay(
				&pbrt.Point3f{0, 0, 0},
				&pbrt.Vector3f{0, 0, 1.0},
				0),
			intersects: true,
			surfaceInteraction: &pbrt.SurfaceInteraction{
				Primitive: primitives[0],
			},
		},
		{
			ray: pbrt.NewRay(
				&pbrt.Point3f{0, 0, 0},
				&pbrt.Vector3f{0, 0, -1.0},
				0),
			intersects:         false,
			surfaceInteraction: pbrt.NewSurfaceInteraction(),
		},
		{
			ray: pbrt.NewRay(
				&pbrt.Point3f{0, 0, 500},
				&pbrt.Vector3f{0, 0, -1.0},
				0),
			intersects: true,
			surfaceInteraction: &pbrt.SurfaceInteraction{
				Primitive: primitives[1],
			},
		},
		{
			ray: pbrt.NewRay(
				&pbrt.Point3f{10, 10, 500},
				&pbrt.Vector3f{0, 0, -1.0},
				0),
			intersects: true,
			surfaceInteraction: &pbrt.SurfaceInteraction{
				Primitive: primitives[2],
			},
		},
	}

	s := NewBVH(primitives, 255, SplitSAH)
	for _, test := range tests {
		t.Run(fmt.Sprintf("%v->%v: %v", test.ray.Origin, test.ray.Direction, test.intersects), func(t *testing.T) {
			si := pbrt.NewSurfaceInteraction()
			intersects := s.Intersect(test.ray, si)
			assert.Equal(t, test.intersects, intersects)
			assert.Equal(t, test.surfaceInteraction.Primitive, si.Primitive)
		})
	}
}

func TestBVH_IntersectP(t *testing.T) {
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

	s := NewBVH(primitives, 2, SplitSAH)
	for _, test := range tests {
		t.Run(fmt.Sprintf("%v->%v: %v", test.ray.Origin, test.ray.Direction, test.expected), func(t *testing.T) {
			assert.Equal(t, test.expected, s.IntersectP(test.ray))
		})
	}
}

func TestPartitionPrimitiveInfoAt(t *testing.T) {
	tests := []struct {
		in                []*BVHPrimitiveInfo
		start, end, pivot int64
		f                 func(a, b *BVHPrimitiveInfo) bool
		expected          []*BVHPrimitiveInfo
	}{
		{
			in: []*BVHPrimitiveInfo{
				{
					primitiveNumber: 5,
					centroid:        &pbrt.Point3f{5, 0, 0},
				},
				{
					primitiveNumber: 4,
					centroid:        &pbrt.Point3f{4, 0, 0},
				},
				{
					primitiveNumber: 3,
					centroid:        &pbrt.Point3f{3, 0, 0},
				},
				{
					primitiveNumber: 2,
					centroid:        &pbrt.Point3f{2, 0, 0},
				},
				{
					primitiveNumber: 1,
					centroid:        &pbrt.Point3f{1, 0, 0},
				},
			},
			start: 0,
			end:   4,
			pivot: 2,
			f: func(a, b *BVHPrimitiveInfo) bool {
				return a.centroid.Index(0) < b.centroid.Index(0)
			},
			expected: []*BVHPrimitiveInfo{
				{
					primitiveNumber: 1,
					centroid:        &pbrt.Point3f{1, 0, 0},
				},
				{
					primitiveNumber: 2,
					centroid:        &pbrt.Point3f{2, 0, 0},
				},
				{
					primitiveNumber: 3,
					centroid:        &pbrt.Point3f{3, 0, 0},
				},
				{
					primitiveNumber: 4,
					centroid:        &pbrt.Point3f{4, 0, 0},
				},
				{
					primitiveNumber: 5,
					centroid:        &pbrt.Point3f{5, 0, 0},
				},
			},
		},
		{
			in: []*BVHPrimitiveInfo{
				{
					primitiveNumber: 5,
					centroid:        &pbrt.Point3f{5, 0, 0},
				},
				{
					primitiveNumber: 1,
					centroid:        &pbrt.Point3f{1, 0, 0},
				},
				{
					primitiveNumber: 2,
					centroid:        &pbrt.Point3f{2, 0, 0},
				},
				{
					primitiveNumber: 4,
					centroid:        &pbrt.Point3f{4, 0, 0},
				},
				{
					primitiveNumber: 3,
					centroid:        &pbrt.Point3f{3, 0, 0},
				},
			},
			start: 0,
			end:   4,
			pivot: 1,
			f: func(a, b *BVHPrimitiveInfo) bool {
				return a.centroid.Index(0) < b.centroid.Index(0)
			},
			expected: []*BVHPrimitiveInfo{
				{
					primitiveNumber: 1,
					centroid:        &pbrt.Point3f{1, 0, 0},
				},
				{
					primitiveNumber: 3,
					centroid:        &pbrt.Point3f{3, 0, 0},
				},
				{
					primitiveNumber: 2,
					centroid:        &pbrt.Point3f{2, 0, 0},
				},
				{
					primitiveNumber: 4,
					centroid:        &pbrt.Point3f{4, 0, 0},
				},
				{
					primitiveNumber: 5,
					centroid:        &pbrt.Point3f{5, 0, 0},
				},
			},
		},
	}

	for i, test := range tests {
		t.Run(fmt.Sprintf("Test%d", i), func(t *testing.T) {
			PartitionPrimitiveInfoAt(test.in, test.start, test.end, test.pivot, test.f)

			fmt.Printf("%+v\n", test.in)
			assert.Equal(t, test.expected, test.in)
		})
	}
}
