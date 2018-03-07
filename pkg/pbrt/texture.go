package pbrt


type Texture interface {
	Evaluate(si *SurfaceInteraction) float64
}