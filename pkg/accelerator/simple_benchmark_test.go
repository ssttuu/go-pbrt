package accelerator

import (
	"testing"
)
import "github.com/ssttuu/go-pbrt/pkg/pbrt"

var result bool

func spheres(n int) []pbrt.Primitive {
	var prims []pbrt.Primitive
	for i := 0; i < n; i++ {
		prims = append(prims, pbrt.NewGeometricPrimitive(
			pbrt.NewSphereShape(
				"",
				pbrt.Translate(&pbrt.Vector3f{0, 0, float64(i * 2)}),
				false,
				0.5,
			),
			nil,
		))
	}
	return prims
}

func benchmarkSimple_Intersect(b *testing.B, n int) {
	s := Simple{primitives: spheres(n)}
	dir := &pbrt.Vector3f{-1, -1, -1}
	ray := pbrt.NewRay(&pbrt.Point3f{15, 15, 15}, dir.Normalized(), 0)
	si := pbrt.NewSurfaceInteraction()

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		result = s.Intersect(ray, si)
	}
}

func BenchmarkSimple_Intersect1(b *testing.B)    { benchmarkSimple_Intersect(b, 1) }
func BenchmarkSimple_Intersect10(b *testing.B)   { benchmarkSimple_Intersect(b, 10) }
func BenchmarkSimple_Intersect100(b *testing.B)  { benchmarkSimple_Intersect(b, 100) }
func BenchmarkSimple_Intersect1000(b *testing.B) { benchmarkSimple_Intersect(b, 1000) }
