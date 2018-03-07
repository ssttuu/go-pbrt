package pbrt

type Bounds2 struct {
	min *Point2f
	max *Point2f
}


func (b *Bounds2) Diagonal() *Vector2f {
	return b.max.Sub(b.min)
}

func (b *Bounds2) Area() float64 {
	bbox := b.Diagonal()
	return bbox.X * bbox.Y
}

func (b *Bounds2) MaximumExtent() uint8 {
	d := b.Diagonal()
	if d.X > d.Y {
		return 0
	}
	return 1
}

func (b *Bounds2) Lerp(other *Point2f) *Point2f {
	return &Point2f{
		Lerp(other.X, b.min.X, b.max.X),
		Lerp(other.Y, b.min.Y, b.max.Y),
	}
}

func (b *Bounds2) Offset(other *Point2f) *Vector2f {
	o := other.Sub(b.min)
	if b.max.X > b.min.X {
		o.X /= b.max.X - b.min.X
	}
	if b.max.Y > b.min.Y {
		o.Y /= b.max.Y - b.min.Y
	}
	return o
}




type Bounds3 struct {
	min *Point3f
	max *Point3f
}


func (b *Bounds3) Diagonal() *Vector3f {
	return b.max.Sub(b.min)
}

func (b *Bounds3) Area() float64 {
	bbox := b.Diagonal()
	return bbox.X * bbox.Y
}

func (b *Bounds3) MaximumExtent() uint8 {
	d := b.Diagonal()
	if d.X > d.Y {
		return 0
	}
	return 1
}

func (b *Bounds3) Lerp(other *Point3f) *Point3f {
	return &Point3f{
		Lerp(other.X, b.min.X, b.max.X),
		Lerp(other.Y, b.min.Y, b.max.Y),
		Lerp(other.Z, b.min.Z, b.max.Z),
	}
}

func (b *Bounds3) Offset(other *Point3f) *Vector3f {
	o := other.Sub(b.min)
	if b.max.X > b.min.X {
		o.X /= b.max.X - b.min.X
	}
	if b.max.Y > b.min.Y {
		o.Y /= b.max.Y - b.min.Y
	}
	return o
}

func (b *Bounds3) Union(p *Point3f) *Bounds3 {
	b.min = MinPoint(b.min, p)
	b.max = MaxPoint(b.min, p)
	return b
}
