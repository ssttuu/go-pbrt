package pbrt

import (
	"math"
	"fmt"
)

type TransportMode int

const (
	Radiance   TransportMode = iota + 1
	Importance
)

type Materialer interface {
	ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool)
}

type Material struct {
}

func (m *Material) Bump(d FloatTexture, si *SurfaceInteraction) {
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

type MatteMaterial struct {
	*Material

	Kd             SpectrumTexture
	sigma, bumpMap FloatTexture
}

func NewMatteMaterial(Kd SpectrumTexture, sigma, bumpMap FloatTexture) *MatteMaterial {
	return &MatteMaterial{
		Material: &Material{},
		Kd: Kd,
		sigma: sigma,
		bumpMap: bumpMap,
	}
}

func (m *MatteMaterial) ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool) {
	if m.bumpMap != nil {
		m.Bump(m.bumpMap, si)
	}

	// evaluate textures for MatteMaterial and allocate BRDF
	si.bsdf = NewBSDF(si, 1.0)
	fmt.Println(si.bsdf.ss)
	r := m.Kd.Evaluate(si).Clamp(0, math.Inf(1))
	sig := Clamp(m.sigma.Evaluate(si), 0, 90)
	if !r.IsBlack() {
		if sig == 0 {
			si.bsdf.Add(NewLambertianReflection(r))
		} else {
			si.bsdf.Add(NewOrenNayar(r, sig))
		}
	}

}