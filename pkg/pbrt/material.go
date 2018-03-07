package pbrt

import "math"

type TransportMode int

const (
	Radiance TransportMode = iota + 1
	Importance
)

type Material struct {

}

func (m *Material) Bump(d Texture, si *SurfaceInteraction) {
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

	siEval.point = si.point.Add(si.shading.dpdu.MulScalar(du))

}
