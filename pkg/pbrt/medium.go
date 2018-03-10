package pbrt

type PhaseFunction interface {
	P(wo, wi *Vector3f) float64
	SampleP(wo, wi *Vector3f, u *Point2f) float64
}

type Mediumer interface {
	Tr(r *Ray, s Sampler) Spectrum
	Sample(r *Ray, s Sampler, mi *MediumInteraction) Spectrum
}

// Named MediumInterface in PBRT book.  Avoiding confusing names
type MediumAccessor struct {
	Inside, Outside Mediumer
}

func (ma *MediumAccessor) IsMediumTransition() bool {
	return ma.Inside == ma.Outside // pointer comparison
}
