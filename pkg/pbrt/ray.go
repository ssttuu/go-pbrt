package pbrt

import "github.com/stupschwartz/go-pbrt/pkg/math"

type Ray struct {
	Origin    *Point3f
	Direction *Vector3f
	TMax      float64
	Time      float64
	Medium    Medium
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

type RayDifferential struct {
	*Ray

	hasDifferentials         bool
	rxOrigin, ryOrigin       *Point3f
	rxDirection, ryDirection *Vector3f
}

func NewRayDifferential() *RayDifferential {
	return &RayDifferential{}
}

func NewRayDifferentialFromRay(r *Ray) *RayDifferential {
	return &RayDifferential{
		Ray:              r,
		hasDifferentials: false,
		rxOrigin:         new(Point3f),
		ryOrigin:         new(Point3f),
		rxDirection:      new(Vector3f),
		ryDirection:      new(Vector3f),
	}
}

func (rd *RayDifferential) ScaleDifferentials(s float64) {
	subbin := rd.rxOrigin.Sub(rd.Origin)
	mulscalarin := subbin.MulScalar(s)
	rd.rxOrigin = rd.Origin.Add(mulscalarin)
	rd.ryOrigin = rd.Origin.Add(rd.ryOrigin.Sub(rd.Origin).MulScalar(s))
	rd.rxDirection = rd.Direction.Add(rd.rxDirection.Sub(rd.Direction).MulScalar(s))
	rd.ryDirection = rd.Direction.Add(rd.ryDirection.Sub(rd.Direction).MulScalar(s))
}
