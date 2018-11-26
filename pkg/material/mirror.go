package material

import (
	"math"

	"github.com/stupschwartz/go-pbrt/pkg/pbrt"
)

func NewMirror() *Mirror {
	return &Mirror{
		Kr:      pbrt.NewConstantSpectrumTexture(pbrt.NewSpectrum(0.9)),
		bumpMap: nil,
	}
}

type Mirror struct {
	Kr      pbrt.SpectrumTexture
	bumpMap pbrt.FloatTexture
}

func (m *Mirror) ComputeScatteringFunctions(si *pbrt.SurfaceInteraction, mode pbrt.TransportMode, allowMultipleLobes bool) {
	if m.bumpMap != nil {
		pbrt.Bump(m, m.bumpMap, si)
	}

	si.BSDF = pbrt.NewBSDF(si, 1.0)
	r := m.Kr.Evaluate(si).Clamp(0.0, math.Inf(1))

	if !r.IsBlack() {
		si.BSDF.Add(pbrt.NewSpecularReflection(r, &pbrt.FresnelNoOp{}))
	}
}
