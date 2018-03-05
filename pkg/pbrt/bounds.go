package pbrt

type Bounds2 struct {
	min *Point2
	max *Point2
}


func (b *Bounds2) Diagonal() *Vector2 {
	return b.max.Sub(b.min)
}

func (b *Bounds2) Area() float64 {
	bbox := b.Diagonal()
	return bbox.x * bbox.y
}

func (b *Bounds2) MaximumExtent() uint8 {
	d := b.Diagonal()
	if d.x > d.y {
		return 0
	}
	return 1
}

func (b *Bounds2) Lerp(other *Point2) *Point2 {
	return &Point2{
		Lerp(other.x, b.min.x, b.max.x),
		Lerp(other.y, b.min.y, b.max.y),
	}
}

func (b *Bounds2) Offset(other *Point2) *Vector2 {
	o := other.Sub(b.min)
	if b.max.x > b.min.x {
		o.x /= b.max.x - b.min.x
	}
	if b.max.y > b.min.y {
		o.y /= b.max.y - b.min.y
	}
	return o
}




type Bounds3 struct {
	min *Point3
	max *Point3
}


func (b *Bounds3) Diagonal() *Vector3 {
	return b.max.Sub(b.min)
}

func (b *Bounds3) Area() float64 {
	bbox := b.Diagonal()
	return bbox.x * bbox.y
}

func (b *Bounds3) MaximumExtent() uint8 {
	d := b.Diagonal()
	if d.x > d.y {
		return 0
	}
	return 1
}

func (b *Bounds3) Lerp(other *Point2) *Point2 {
	return &Point2{
		Lerp(other.x, b.min.x, b.max.x),
		Lerp(other.y, b.min.y, b.max.y),
	}
}

func (b *Bounds3) Offset(other *Point2) *Vector2 {
	o := other.Sub(b.min)
	if b.max.x > b.min.x {
		o.x /= b.max.x - b.min.x
	}
	if b.max.y > b.min.y {
		o.y /= b.max.y - b.min.y
	}
	return o
}
