package materials

import (
	"github.com/ssttuu/go-pbrt/pkg/math"
	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

type MatteMaterial struct {
	Kd             pbrt.SpectrumTexture
	sigma, bumpMap pbrt.FloatTexture
}

func NewMatteMaterial(Kd pbrt.SpectrumTexture, sigma, bumpMap pbrt.FloatTexture) *MatteMaterial {
	return &MatteMaterial{
		Kd:      Kd,
		sigma:   sigma,
		bumpMap: bumpMap,
	}
}

func (m *MatteMaterial) ComputeScatteringFunctions(si *pbrt.SurfaceInteraction, mode pbrt.TransportMode, allowMultipleLobes bool) {
	if m.bumpMap != nil {
		pbrt.Bump(m, m.bumpMap, si)
	}

	// evaluate textures for MatteMaterial and allocate BRDF
	si.BSDF = pbrt.NewBSDF(si, 1.0)
	r := m.Kd.Evaluate(si).Clamp(0, math.Inf(1))
	sig := math.Clamp(m.sigma.Evaluate(si), 0, 90)
	if !r.IsBlack() {
		if sig == 0 {
			si.BSDF.Add(pbrt.NewLambertianReflection(r))
		} else {
			si.BSDF.Add(pbrt.NewOrenNayar(r, sig))
		}
	}
}
