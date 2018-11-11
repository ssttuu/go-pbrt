package pbrt

type SpectrumTexture interface {
	Evaluate(si *SurfaceInteraction) Spectrum
}

type FloatTexture interface {
	Evaluate(si *SurfaceInteraction) float64
}

type ConstantSpectrumTexture struct {
	value Spectrum
}

func NewConstantSpectrumTexture(s Spectrum) SpectrumTexture {
	return &ConstantSpectrumTexture{
		value: s,
	}
}

func (t *ConstantSpectrumTexture) Evaluate(si *SurfaceInteraction) Spectrum {
	return t.value
}

type ConstantFloatTexture struct {
	value float64
}

func NewConstantFloatTexture(v float64) FloatTexture {
	return &ConstantFloatTexture{
		value: v,
	}
}

func (t *ConstantFloatTexture) Evaluate(si *SurfaceInteraction) float64 {
	return t.value
}
