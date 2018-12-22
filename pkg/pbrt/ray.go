package pbrt

import "github.com/ssttuu/go-pbrt/pkg/math"

type Ray struct {
	Origin    *Point3f
	Direction *Vector3f
	TMax      float64
	Time      float64
	Medium    Medium

	HasDifferentials         bool
	rxOrigin, ryOrigin       *Point3f
	rxDirection, ryDirection *Vector3f
}

func NewRay(origin *Point3f, direction *Vector3f, time float64) *Ray {
	return &Ray{
		Origin:    origin,
		Direction: direction,
		TMax:      math.Infinity,
		Time:      time,
		Medium:    nil,
	}
}

func NewRayWithMedium(origin *Point3f, direction *Vector3f, time float64, medium Medium) *Ray {
	return &Ray{
		Origin:    origin,
		Direction: direction,
		TMax:      math.Infinity,
		Time:      time,
		Medium:    medium,
	}
}

func NewRayDifferentialFromRay(r *Ray) *Ray {
	return &Ray{
		Origin:    r.Origin,
		Direction: r.Direction,
		TMax:      r.TMax,
		Time:      r.Time,
		Medium:    r.Medium,

		HasDifferentials: false,
		rxOrigin:         new(Point3f),
		ryOrigin:         new(Point3f),
		rxDirection:      new(Vector3f),
		ryDirection:      new(Vector3f),
	}
}

func (r *Ray) PointAt(t float64) *Point3f {
	return r.Origin.Add(r.Direction.MulScalar(t))
}

func OffsetRayOrigin(p *Point3f, pError *Vector3f, n *Normal3f, w *Vector3f) *Point3f {
	d := n.Abs().Dot(pError) * 1024.0 // TODO:
	offset := n.MulScalar(d)
	if w.Dot(n) < 0 {
		offset = offset.MulScalar(-1)
	}
	po := p.Add(offset)

	for i := 0; i < 3; i++ {
		if offset.Index(i) > 0 {
			po.SetIndex(i, math.NextFloatUp(po.Index(i)))
		} else if offset.Index(i) < 0 {
			po.SetIndex(i, math.NextFloatDown(po.Index(i)))
		}
	}

	return po
}

func (r *Ray) ScaleDifferentials(s float64) {
	subbin := r.rxOrigin.Sub(r.Origin)
	mulscalarin := subbin.MulScalar(s)
	r.rxOrigin = r.Origin.Add(mulscalarin)
	r.ryOrigin = r.Origin.Add(r.ryOrigin.Sub(r.Origin).MulScalar(s))
	r.rxDirection = r.Direction.Add(r.rxDirection.Sub(r.Direction).MulScalar(s))
	r.ryDirection = r.Direction.Add(r.ryDirection.Sub(r.Direction).MulScalar(s))
}
