//go:generate mockgen -source=material.go -destination=material.mock.go -package=pbrt

package pbrt

import "github.com/ssttuu/go-pbrt/pkg/math"

type TransportMode int

const (
	Radiance TransportMode = iota + 1
	Importance
)

type Material interface {
	ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool)
}

func Bump(m Material, d FloatTexture, si *SurfaceInteraction) {
	// compute offset positions and evaluate displacement texture
	siEval := *si

	du := 0.5 * (math.Abs(si.dudx) + math.Abs(si.dudy))

	// The most common reason for du to be zero is for ray that start from
	// light sources, where no differentials are available. In this case,
	// we try to choose a small enough du so that we still get a decently
	// accurate bump value.

	if du == 0 {
		du = 0.0005
	}

	siEval.Point = si.Point.Add(si.Shading.dpdu.MulScalar(du))
}
