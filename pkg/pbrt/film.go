package pbrt

import (
	"sync/atomic"
	"sync"
	"math"
)

type FilmTilePixel struct {
	contribSum      Spectrum
	filterWeightSum float64
}

type Pixel struct {
	value           [3]float64
	filterWeightSum float64
	splat           [3]*atomic.Value
	pad             float64
}

type Film struct {
	FullResolution     *Point2i
	Diagonal           float64
	Filter             Filterer
	CroppedPixelBounds *Bounds2i

	pixels             []*Pixel
	filterTableWidth   int
	filterTable        []float64
	mutex              sync.Mutex
	scale              float64
	maxSampleLuminance float64
}

func NewFilm(resolution *Point2i, cropWindow Bounds2f, filter Filterer, diagonal, scale, maxSampleLuminance float64) *Film {
	croppedPixelBounds := &Bounds2i{
		min: &Point2i{int64(math.Ceil(float64(resolution.X) * cropWindow.min.X)), int64(math.Ceil(float64(resolution.Y) * cropWindow.min.Y))},
		max: &Point2i{int64(math.Ceil(float64(resolution.X) * cropWindow.max.X)), int64(math.Ceil(float64(resolution.Y) * cropWindow.max.Y))},
	}
	filterTableWidth := 16.0
	f := &Film{
		FullResolution:     resolution,
		CroppedPixelBounds: croppedPixelBounds,
		Diagonal:           diagonal * 0.001,
		Filter:             filter,
		scale:              scale,
		maxSampleLuminance: maxSampleLuminance,
		pixels:             make([]*Pixel, croppedPixelBounds.Area()),
		filterTableWidth:   int(filterTableWidth),
		filterTable:        make([]float64, int(filterTableWidth*filterTableWidth)),
	}

	offset := 0
	filterRadius := filter.GetRadius()
	for y := 0; y < f.filterTableWidth; y++ {
		for x := 0; x < f.filterTableWidth; x++ {
			p := &Point2f{
				X: (float64(x) + 0.5) * filterRadius.X / filterTableWidth,
				Y: (float64(y) + 0.5) * filterRadius.Y / filterTableWidth,
			}
			f.filterTable[offset] = filter.Evaluate(p)

			offset++
		}
	}

	return f
}

func (f *Film) getPixel(p *Point2i) *Pixel {
	width := f.CroppedPixelBounds.max.X - f.CroppedPixelBounds.min.X
	offset := (p.X - f.CroppedPixelBounds.min.X) + (p.Y-f.CroppedPixelBounds.min.Y)*width
	return f.pixels[offset]
}

func (f *Film) GetSampleBounds() *Bounds2i {
	min := NewPoint2fFromPoint2i(f.CroppedPixelBounds.min)
	max := NewPoint2fFromPoint2i(f.CroppedPixelBounds.max)
	return &Bounds2i{
		min: NewPoint2iFromPoint2f(min.Add(&Vector2f{0.5, 0.5}).Sub(f.Filter.GetRadius()).Floor()),
		max: NewPoint2iFromPoint2f(max.Sub(&Vector2f{0.5, 0.5}).Add(f.Filter.GetRadius()).Ceil()),
	}
}

func (f *Film) GetPhysicalExtent() *Bounds2f {
	aspect := float64(f.FullResolution.Y) / float64(f.FullResolution.X)
	x := math.Sqrt(f.Diagonal * f.Diagonal / (1.0 + aspect*aspect))
	y := aspect * x
	return &Bounds2f{
		min: &Point2f{-x / 2, -y / 2},
		max: &Point2f{x / 2, y / 2},
	}
}

func (f *Film) GetFilmTile(sampleBounds *Bounds2i) *FilmTile {
	halfPixel := &Vector2f{0.5, 0.5}
	floatBounds := &Bounds2f{NewPoint2fFromPoint2i(sampleBounds.min), NewPoint2fFromPoint2i(sampleBounds.max)}
	p0 := NewPoint2iFromPoint2f(floatBounds.min.Sub(halfPixel).Sub(f.Filter.GetRadius()).Ceil())
	p1 := NewPoint2iFromPoint2f(floatBounds.max.Sub(halfPixel).Add(f.Filter.GetRadius()).Floor()).Add(&Point2i{1, 1})
	tilePixelBounds := f.CroppedPixelBounds.Intersect(&Bounds2i{p0, p1})
	return NewFilmTile(tilePixelBounds, f.Filter.GetRadius(), f.filterTable, f.filterTableWidth, f.maxSampleLuminance)
}

func (f *Film) MergeFilmTile(tile *FilmTile) {
	f.mutex.Lock()
	// TODO

	f.mutex.Unlock()
}

func (f *Film) SetImage(img Spectrum) {
	// TODO
}

func (f *Film) AddSplat(p *Point2f, v Spectrum) {
	// TODO
}

func (f *Film) WriteImage(splatScale float64) {
	// TODO
}

func (f *Film) Clear() {
	// TODO
}

type FilmTile struct {
	pixelBounds                   *Bounds2i
	filterRadius, invFilterRadius *Vector2f
	filterTable                   []float64
	filterTableSize               int
	pixels                        []*FilmTilePixel
	maxSampleLuminance            float64
}

func NewFilmTile(pixelBounds *Bounds2i, filterRadius *Vector2f, filterTable []float64, filterTableSize int, maxSampleLuminance float64) *FilmTile {
	return &FilmTile{
		pixelBounds:        pixelBounds,
		filterRadius:       filterRadius,
		filterTable:        filterTable,
		filterTableSize:    filterTableSize,
		maxSampleLuminance: maxSampleLuminance,
		pixels:             make([]*FilmTilePixel, int64(math.Max(0, float64(pixelBounds.Area())))),
	}
}

func (f *FilmTile) AddSample(pFilm *Point2f, L Spectrum, sampleWeight float64) {

}

func (f *FilmTile) GetPixel(p *Point2i) *FilmTilePixel {
	width := f.pixelBounds.max.X - f.pixelBounds.min.X
	offset := (p.X - f.pixelBounds.min.X) + (p.Y-f.pixelBounds.min.Y)*width
	return f.pixels[offset]
}

func (f *FilmTile) GetPixelBounds() *Bounds2i {
	return f.pixelBounds
}
