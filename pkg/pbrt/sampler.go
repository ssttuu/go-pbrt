//go:generate mockgen -source=sampler.go -destination=sampler.mock.go -package=pbrt

package pbrt

type Sampler interface {
	GetSamplesPerPixel() int32
	Clone(seed uint64) Sampler
	Get1D() float64
	Get2D() *Point2f
	Get1DArray(n int32) []float64
	Get2DArray(n int32) []Point2f
	GetCameraSample(pRaster *Point2i) *CameraSample
	StartPixel(p *Point2i)
	StartNextSample() bool
	RoundCount(n int32) int32
	Request2DArray(n int32)
}

func GetCameraSample(s Sampler, pRaster *Point2i) *CameraSample {
	return &CameraSample{
		pFilm: NewPoint2fFromPoint2i(pRaster).Add(s.Get2D()),
		time:  s.Get1D(),
		pLens: s.Get2D(),
	}
}

