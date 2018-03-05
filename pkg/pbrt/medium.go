package pbrt


type Medium interface {
	Tr() Spectrum
	Sample() Spectrum
}

// Named MediumInterface in PBRT book.  Avoiding confusing names
type MediumAccessor struct {
	Inside, Outside *Medium
}

func (ma *MediumAccessor) IsMediumTransition() bool {
	return ma.Inside == ma.Outside // pointer comparison
}
