package pbrt

type Filterer interface {
	GetRadius() *Point2f
	Evaluate(p *Point2f) float64
}

type Filter struct {
	radius *Point2f
}

func (f *Filter) GetRadius() *Point2f {
	return f.radius
}

// TODO: move from pbrt core

type BoxFilter struct {
	*Filter
}

func NewBoxFilter(radius *Point2f) *BoxFilter {
	return &BoxFilter{
		Filter: &Filter{radius: radius},
	}
}

func (f *BoxFilter) Evaluate(p *Point2f) float64 {
	return 1.0
}