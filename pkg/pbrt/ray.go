package pbrt

type Ray struct {
	origin    *Point3
	direction *Vector3
	tMax      float64
	time      float64
	medium    *Medium
}

func OffsetRayOrigin(p *Point3, pError *Vector3, n *Normal3, w *Vector3) *Point3 {
	d := n.Abs().Dot(pError)
	offset := n.MulScalar(d)
	if w.Dot(n.ToVector3()) < 0 {
		offset = offset.MulScalar(-1)
	}
	po := p.Add(offset)

	for i := 0; i < 3; i++ {
		if (offset.Index(i) > 0) {
			po.SetIndex(i, NextFloatUp(po.Index(i)))
		} else if offset.Index(i) < 0 {
			po.SetIndex(i, NextFloatDown(po.Index(i)))
		}
	}

	return po

}

type RayDifferential struct {
	origin                   *Vector3
	direction                *Vector3
	tMax                     float64
	time                     float64
	medium                   Medium

	hasDifferentials         bool
	rxOrigin, ryOrigin       *Point3
	rxDirection, ryDirection *Vector3
}
