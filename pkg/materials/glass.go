package materials

import (
	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

type Glass struct {
	Kr, Kt                *pbrt.ConstantSpectrumTexture
	uRoughness, vRoughness pbrt.FloatTexture
	index                  pbrt.FloatTexture
	bumpMap                pbrt.FloatTexture
	remapRoughness         bool
}

func NewGlass(Kr, Kt *pbrt.ConstantSpectrumTexture, uRoughness, vRoughness, index, bumpMap pbrt.FloatTexture) pbrt.Material {
	return &Glass{
		Kr:             Kr,
		Kt:             Kt,
		uRoughness:     uRoughness,
		vRoughness:     vRoughness,
		index:          index,
		bumpMap:        bumpMap,
		remapRoughness: false,
	}
}

func (g *Glass) ComputeScatteringFunctions(si *pbrt.SurfaceInteraction, mode pbrt.TransportMode, allowMultipleLobes bool) {
	if g.bumpMap != nil {
		pbrt.Bump(g, g.bumpMap, si)
	}

	eta := g.index.Evaluate(si)
	uRough := g.uRoughness.Evaluate(si)
	vRough := g.vRoughness.Evaluate(si)
	R := g.Kr.Evaluate(si).Clamp(0, 1)
	T := g.Kt.Evaluate(si).Clamp(0, 1)
	// initialize bsdf for smooth or rough dielectric
	si.BSDF = pbrt.NewBSDF(si, eta)

	if R.IsBlack() && T.IsBlack() {
		return
	}

	isSpecular := uRough == 0 && vRough == 0
	if isSpecular && allowMultipleLobes {
		si.BSDF.Add(pbrt.NewFresnelSpecular(R, T, 1.0, eta, mode))
	} else {
		if g.remapRoughness {
			uRough = pbrt.RoughnessToAlpha(uRough)
			vRough = pbrt.RoughnessToAlpha(vRough)
		}

		var distrib pbrt.Distribution
		if !isSpecular {
			distrib = pbrt.NewTrowbridgeReitz(uRough, vRough)
		}

		if !R.IsBlack() {
			fresnel := pbrt.NewFresnelDielectric(1.0, eta)
			if isSpecular {
				si.BSDF.Add(pbrt.NewSpecularReflection(R, fresnel))
			} else {
				si.BSDF.Add(pbrt.NewMicrofacetReflection(R, distrib, fresnel))
			}
		}
		if !T.IsBlack() {
			if isSpecular {
				si.BSDF.Add(pbrt.NewSpecularTransmission(T, 1.0, eta, mode))
			} else {
				si.BSDF.Add(pbrt.NewMicrofacetTransmission(T, distrib, 1.0, eta, mode))
			}
		}

	}
}
