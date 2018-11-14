package pbrt

import (
	"sync/atomic"
	"sync"
	"math"
	"os"
	"image/png"
	"image"
	"image/color"
	"log"
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
	Filter             Filter
	CroppedPixelBounds *Bounds2i
	Filename           string

	pixels             []Pixel
	filterTableWidth   int
	filterTable        []float64
	mutex              sync.Mutex
	scale              float64
	maxSampleLuminance float64
}

func NewFilm(filename string, resolution *Point2i, cropWindow *Bounds2f, filter Filter, diagonal, scale, maxSampleLuminance float64) *Film {
	croppedPixelBounds := &Bounds2i{
		Min: &Point2i{int64(math.Ceil(float64(resolution.X) * cropWindow.Min.X)), int64(math.Ceil(float64(resolution.Y) * cropWindow.Min.Y))},
		Max: &Point2i{int64(math.Ceil(float64(resolution.X) * cropWindow.Max.X)), int64(math.Ceil(float64(resolution.Y) * cropWindow.Max.Y))},
	}
	filterTableWidth := 16.0
	f := &Film{
		Filename:           filename,
		FullResolution:     resolution,
		CroppedPixelBounds: croppedPixelBounds,
		Diagonal:           diagonal * 0.001,
		Filter:             filter,
		scale:              scale,
		maxSampleLuminance: maxSampleLuminance,
		pixels:             make([]Pixel, croppedPixelBounds.Area()),
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
	width := f.CroppedPixelBounds.Max.X - f.CroppedPixelBounds.Min.X
	offset := (p.X - f.CroppedPixelBounds.Min.X) + (p.Y-f.CroppedPixelBounds.Min.Y)*width
	return &f.pixels[offset]
}

func (f *Film) GetSampleBounds() *Bounds2i {
	//min := NewPoint2fFromPoint2i(f.CroppedPixelBounds.Min)
	//max := NewPoint2fFromPoint2i(f.CroppedPixelBounds.Max)
	// TODO
	//return &Bounds2i{
	//	Min: NewPoint2iFromPoint2f(min.Add(&Vector2f{0.5, 0.5}).Sub(f.filter.GetRadius()).Floor()),
	//	Max: NewPoint2iFromPoint2f(max.Sub(&Vector2f{0.5, 0.5}).Add(f.filter.GetRadius()).Ceil()),
	//}

	return f.CroppedPixelBounds
}

func (f *Film) GetPhysicalExtent() *Bounds2f {
	aspect := float64(f.FullResolution.Y) / float64(f.FullResolution.X)
	x := math.Sqrt(f.Diagonal * f.Diagonal / (1.0 + aspect*aspect))
	y := aspect * x
	return &Bounds2f{
		Min: &Point2f{-x / 2, -y / 2},
		Max: &Point2f{x / 2, y / 2},
	}
}

func (f *Film) GetFilmTile(sampleBounds *Bounds2i) *FilmTile {
	halfPixel := &Vector2f{0.5, 0.5}
	floatBounds := &Bounds2f{NewPoint2fFromPoint2i(sampleBounds.Min), NewPoint2fFromPoint2i(sampleBounds.Max)}
	p0 := NewPoint2iFromPoint2f(floatBounds.Min.Sub(halfPixel).Sub(f.Filter.GetRadius()).Ceil())
	p1 := NewPoint2iFromPoint2f(floatBounds.Max.Sub(halfPixel).Add(f.Filter.GetRadius()).Floor()).Add(&Point2i{1, 1})
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

	imgPng := image.NewNRGBA(image.Rect(
		int(f.CroppedPixelBounds.Min.X),
		int(f.CroppedPixelBounds.Min.Y),
		int(f.CroppedPixelBounds.Max.X),
		int(f.CroppedPixelBounds.Max.Y),
	))

	for y := f.CroppedPixelBounds.Min.Y; y < f.CroppedPixelBounds.Max.Y; y++ {
		for x := f.CroppedPixelBounds.Min.X; x < f.CroppedPixelBounds.Max.X; x++ {

			pixel := f.getPixel(&Point2i{x, y})

			imgPng.Set(int(x), int(y), color.NRGBA{
				R: uint8(pixel.value[0] * 255),
				G: uint8(pixel.value[1] * 255),
				B: uint8(pixel.value[2] * 255),
				A: 255,
			})
		}
	}

	fname, err := os.Create(f.Filename)
	if err != nil {
		log.Fatal(err)
	}

	if err := png.Encode(fname, imgPng); err != nil {
		fname.Close()
		log.Fatal(err)
	}

	if err := fname.Close(); err != nil {
		log.Fatal(err)
	}

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
	//if L.Y() > f.maxSampleLuminance {
	//	L = L.MulScalar(f.maxSampleLuminance / L.Y())
	//}
	//
	//// compute sample's raster bounds
	//pFilmDiscrete := pFilm.Sub(&Vector2f{0.5, 0.5})
	//p0f := pFilmDiscrete.Sub(f.filterRadius).Ceil()
	//p1f := pFilmDiscrete.Add(f.filterRadius).Floor().Add(&Point2f{1, 1})
	//p0 := &Point2i{int64(math.Max(p0f.X, float64(f.pixelBounds.Min.X))), int64(math.Max(p0f.Y, float64(f.pixelBounds.Min.Y)))}
	//p1 := &Point2i{int64(math.Min(p1f.X, float64(f.pixelBounds.Max.X))), int64(math.Min(p1f.Y, float64(f.pixelBounds.Max.Y)))}

	// loop over filter support and add sample to pixel arrays

	// precompute x and y filter table offsets

}

func (f *FilmTile) GetPixel(p *Point2i) *FilmTilePixel {
	width := f.pixelBounds.Max.X - f.pixelBounds.Min.X
	offset := (p.X - f.pixelBounds.Min.X) + (p.Y-f.pixelBounds.Min.Y)*width
	return f.pixels[offset]
}

func (f *FilmTile) GetPixelBounds() *Bounds2i {
	return f.pixelBounds
}
