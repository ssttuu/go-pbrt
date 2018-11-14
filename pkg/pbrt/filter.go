//go:generate mockgen -source=filter.go -destination=filter.mock.go -package=pbrt

package pbrt

type Filter interface {
	GetRadius() *Point2f
	Evaluate(p *Point2f) float64
}

type filter struct {
	radius *Point2f
}

func (f *filter) GetRadius() *Point2f {
	return f.radius
}

// TODO: move from pbrt core

type BoxFilter struct {
	*filter
}

func NewBoxFilter(radius *Point2f) *BoxFilter {
	return &BoxFilter{
		filter: &filter{radius: radius},
	}
}

func (f *BoxFilter) Evaluate(p *Point2f) float64 {
	return 1.0
}