package geometry


type Medium interface {
	Tr() Spectrum
	Sample() Spectrum
}
