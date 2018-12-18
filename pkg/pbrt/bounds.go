package pbrt

import "github.com/ssttuu/go-pbrt/pkg/math"

//////////////
// Bounds2f //
//////////////

type Bounds2f struct {
	Min *Point2f
	Max *Point2f
}

func (b *Bounds2f) Diagonal() *Vector2f {
	return b.Max.Sub(b.Min)
}

func (b *Bounds2f) Area() float64 {
	bbox := b.Diagonal()
	return bbox.X * bbox.Y
}

func (b *Bounds2f) MaximumExtent() uint8 {
	d := b.Diagonal()
	if d.X > d.Y {
		return 0
	}
	return 1
}

func (b *Bounds2f) Lerp(other *Point2f) *Point2f {
	return &Point2f{
		math.Lerp(other.X, b.Min.X, b.Max.X),
		math.Lerp(other.Y, b.Min.Y, b.Max.Y),
	}
}

func (b *Bounds2f) Offset(other *Point2f) *Vector2f {
	o := other.Sub(b.Min)
	if b.Max.X > b.Min.X {
		o.X /= b.Max.X - b.Min.X
	}
	if b.Max.Y > b.Min.Y {
		o.Y /= b.Max.Y - b.Min.Y
	}
	return o
}

//////////////
// Bounds2i //
//////////////

type Bounds2i struct {
	Min *Point2i
	Max *Point2i
}

func (b *Bounds2i) Diagonal() *Vector2i {
	return b.Max.Sub(b.Min)
}

func (b *Bounds2i) Area() int64 {
	bbox := b.Diagonal()
	return bbox.X * bbox.Y
}

func (b *Bounds2i) MaximumExtent() uint8 {
	d := b.Diagonal()
	if d.X > d.Y {
		return 0
	}
	return 1
}

func (b *Bounds2i) Lerp(other *Point2i) *Point2i {
	return &Point2i{
		int64(math.Lerp(float64(other.X), float64(b.Min.X), float64(b.Max.X))),
		int64(math.Lerp(float64(other.Y), float64(b.Min.Y), float64(b.Max.Y))),
	}
}

func (b *Bounds2i) Offset(other *Point2i) *Vector2i {
	o := other.Sub(b.Min)
	if b.Max.X > b.Min.X {
		o.X /= b.Max.X - b.Min.X
	}
	if b.Max.Y > b.Min.Y {
		o.Y /= b.Max.Y - b.Min.Y
	}
	return o
}

func (b *Bounds2i) Intersect(other *Bounds2i) *Bounds2i {
	return &Bounds2i{
		Min: &Point2i{int64(math.Max(float64(b.Min.X), float64(other.Min.X))), int64(math.Max(float64(b.Min.Y), float64(other.Min.Y)))},
		Max: &Point2i{int64(math.Min(float64(b.Max.X), float64(other.Max.X))), int64(math.Min(float64(b.Max.Y), float64(other.Max.Y)))},
	}
}

type Bounds3 struct {
	Min *Point3f
	Max *Point3f
}

func (b *Bounds3) BoundingSphere() (center *Point3f, radius float64) {
	center = b.Min.Add(b.Max).DivScalar(2.0)
	radius = 0
	if Inside(center, b) {
		radius = center.Distance(b.Max)
	}
	return center, radius
}

func (b *Bounds3) Corner(corner int) *Point3f {
	return &Point3f{
		X: b.Index(corner & 1).X,
		Y: b.Index((corner & 2) / 2).Y,
		Z: b.Index((corner & 4) / 4).Z,
	}
}

func (b *Bounds3) Index(i int) *Point3f {
	if i == 0 {
		return b.Min
	}
	return b.Max
}

func (b *Bounds3) Diagonal() *Vector3f {
	return b.Max.Sub(b.Min)
}

func (b *Bounds3) SurfaceArea() float64 {
	d := b.Diagonal()
	return 2 * (d.X*d.Y + d.X*d.Z + d.Y*d.Z)
}

func (b *Bounds3) MaximumExtent() int {
	d := b.Diagonal()
	if d.X > d.Y && d.X > d.Z {
		return 0
	} else if d.Y > d.Z {
		return 1
	} else {
		return 2
	}
}

func (b *Bounds3) IntersectP(r *Ray, invDir *Vector3f, directionIsNegative [3]int) bool {
	// check for ray intersection against x and y slabs
	tMin := (b.Index(directionIsNegative[0]).X - r.Origin.X) * invDir.X
	tMax := (b.Index(1-directionIsNegative[0]).X - r.Origin.X) * invDir.X
	tyMin := (b.Index(directionIsNegative[1]).Y - r.Origin.Y) * invDir.Y
	tyMax := (b.Index(1-directionIsNegative[1]).Y - r.Origin.Y) * invDir.Y

	// update tMax and tyMax to ensure robust bounds intersection
	tMax *= 1 + 2*math.Gamma(3)
	tyMax *= 1 + 2*math.Gamma(3)
	if tMin > tyMax || tyMin > tMax {
		return false
	}
	if tyMin > tMin {
		tMin = tyMin
	}
	if tyMax < tMax {
		tMax = tyMax
	}

	// check for ray intersection against z slab
	tzMin := (b.Index(directionIsNegative[2]).Z - r.Origin.Z) * invDir.Z
	tzMax := (b.Index(1-directionIsNegative[2]).Z - r.Origin.Z) * invDir.Z

	// update tzMax to ensure robust bounds intersection
	tzMax *= 1 + 2*math.Gamma(3)
	if tMin > tzMax || tzMin > tMax {
		return false
	}
	if tzMin > tMin {
		tMin = tzMin
	}
	if tzMax < tMax {
		tMax = tzMax
	}
	return tMin < r.TMax && tMax > 0
}

func (b *Bounds3) Lerp(other *Point3f) *Point3f {
	return &Point3f{
		math.Lerp(other.X, b.Min.X, b.Max.X),
		math.Lerp(other.Y, b.Min.Y, b.Max.Y),
		math.Lerp(other.Z, b.Min.Z, b.Max.Z),
	}
}

func (b *Bounds3) Offset(other *Point3f) *Vector3f {
	o := other.Sub(b.Min)
	if b.Max.X > b.Min.X {
		o.X /= b.Max.X - b.Min.X
	}
	if b.Max.Y > b.Min.Y {
		o.Y /= b.Max.Y - b.Min.Y
	}
	if b.Max.Z > b.Min.Z {
		o.Z /= b.Max.Z - b.Min.Z
	}
	return o
}

func (b *Bounds3) UnionPoint(p *Point3f) *Bounds3 {
	if b.Min == nil {
		b.Min = p
	}
	if b.Max == nil {
		b.Max = p
	}
	b.Min = MinPoint(b.Min, p)
	b.Max = MaxPoint(b.Max, p)
	return b
}

func (b *Bounds3) Union(b2 *Bounds3) *Bounds3 {
	if b.Min == nil {
		b.Min = b2.Min
	}

	if b.Max == nil {
		b.Max = b2.Max
	}

	if b2.Min == nil || b2.Max == nil {
		return b
	}

	b.Min = MinPoint(b.Min, b2.Min)
	b.Max = MaxPoint(b.Max, b2.Max)

	return b
}
