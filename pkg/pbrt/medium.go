//go:generate mockgen -source=medium.go -destination=medium.mock.go -package=pbrt

package pbrt

type PhaseFunction interface {
	P(wo, wi *Vector3f) float64
	SampleP(wo, wi *Vector3f, u *Point2f) float64
}

type Medium interface {
	Tr(r *Ray, s Sampler) Spectrum
	Sample(r *Ray, s Sampler, mi *MediumInteraction) Spectrum
}

// Named MediumInterface in PBRT book.  Avoiding confusing names
type MediumAccessor struct {
	Inside, Outside Medium
}

func (ma *MediumAccessor) IsMediumTransition() bool {
	if ma.Inside == nil {
		return false
	}
	return ma.Inside == ma.Outside // pointer comparison
}
