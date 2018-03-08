package pbrt

type Filterer interface {
	GetRadius() *Point2f
	Evaluate(p *Point2f) float64
}