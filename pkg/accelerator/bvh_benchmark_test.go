package accelerator

import (
	"testing"

	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

func benchmarkBVH_Intersect(b *testing.B, n int) {
	s := NewBVH(spheres(n), 255, SplitSAH)
	dir := &pbrt.Vector3f{-1, -1, -1}
	ray := pbrt.NewRay(&pbrt.Point3f{15, 15, 15}, dir.Normalized(), 0)
	si := pbrt.NewSurfaceInteraction()

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		result = s.Intersect(ray, si)
	}
}

func BenchmarkBVH_Intersect1(b *testing.B)    { benchmarkBVH_Intersect(b, 1) }
func BenchmarkBVH_Intersect10(b *testing.B)   { benchmarkBVH_Intersect(b, 10) }
func BenchmarkBVH_Intersect100(b *testing.B)  { benchmarkBVH_Intersect(b, 100) }
func BenchmarkBVH_Intersect1000(b *testing.B) { benchmarkBVH_Intersect(b, 1000) }
