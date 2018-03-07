package pbrt

type Ray struct {
	origin    *Point3f
	direction *Vector3f
	tMax      float64
	time      float64
	medium    Mediumer
}

func NewRay(origin *Point3f, direction *Vector3f, time float64) *Ray {
	return &Ray{
		origin: origin,
		direction: direction,
		tMax: Infinity,
		time: time,
		medium: nil,
	}
}

func NewRayWithMedium(origin *Point3f, direction *Vector3f, time float64, medium Mediumer) *Ray {
	return &Ray{
		origin: origin,
		direction: direction,
		tMax: Infinity,
		time: time,
		medium: medium,
	}
}

func (r *Ray) PointAt(t float64) *Point3f {
	return r.origin.Mul(r.direction).MulScalar(t)
}


func OffsetRayOrigin(p *Point3f, pError *Vector3f, n *Normal3f, w *Vector3f) *Point3f {
	d := n.Abs().Dot(pError)
	offset := n.MulScalar(d)
	if w.Dot(n) < 0 {
		offset = offset.MulScalar(-1)
	}
	po := p.Add(offset)

	for i := 0; i < 3; i++ {
		if offset.Index(i) > 0 {
			po.SetIndex(i, NextFloatUp(po.Index(i)))
		} else if offset.Index(i) < 0 {
			po.SetIndex(i, NextFloatDown(po.Index(i)))
		}
	}

	return po

}

type RayDifferential struct {
	*Ray

	hasDifferentials         bool
	rxOrigin, ryOrigin       *Point3f
	rxDirection, ryDirection *Vector3f
}
