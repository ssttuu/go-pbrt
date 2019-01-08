package pbrt

type BSSRDF interface {
	S(pi SurfaceInteraction, wi Vector3f) Spectrum
	SampleS(scene Scene, u1 float64, u2 *Point2f, si *SurfaceInteraction) (s Spectrum, pdf float64)
}

type bssrdf struct {
	po  *SurfaceInteraction
	eta float64
}



