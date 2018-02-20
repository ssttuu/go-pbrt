package geometry

type Bounds2 struct {
	min VectorXY
	max VectorXY
}


func (b *Bounds2) Diagonal() VectorXY {
	return b.max.Subtract(&b.min)
}

func (b *Bounds2) Area() float64 {
	bbox := b.max.Subtract(&b.min)
	return bbox.GetX() * bbox.GetY()
}


