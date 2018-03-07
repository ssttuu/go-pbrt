package pbrt


type Mediumer interface {
	Tr(r *Ray, s *Sampler) Spectrum
	Sample(r *Ray, s *Sampler, arena *MemoryArena, mi *MediumInteraction) Spectrum
}

// Named MediumInterface in PBRT book.  Avoiding confusing names
type MediumAccessor struct {
	Inside, Outside Mediumer
}

func (ma *MediumAccessor) IsMediumTransition() bool {
	return ma.Inside == ma.Outside // pointer comparison
}
