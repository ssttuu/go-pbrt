package shapes

import (
	"github.com/stupschwartz/go-pbrt/pkg/math"
	"github.com/stupschwartz/go-pbrt/pkg/pbrt"
)

type shape struct {
	objectToWorld, worldToObject *pbrt.Transform
	reverseOrientation           bool
	transformSwapsHandedness     bool
}

type Disk struct {
	shape       *shape
	height      float64
	radius      float64
	innerRadius float64
	phiMax      float64
}

func NewDisk(objectToWorld *pbrt.Transform, height, radius, innerRadius, phiMax float64) pbrt.Shape {
	return &Disk{
		shape: &shape{
			objectToWorld:            objectToWorld,
			worldToObject:            objectToWorld.Inverse(),
			reverseOrientation:       false,
			transformSwapsHandedness: false,
		},
		height:      height,
		radius:      radius,
		innerRadius: innerRadius,
		phiMax:      math.Radians(math.Clamp(phiMax, 0, 360)) ,
	}
}

func (d *Disk) GetName() string {
	return "Disk"
}

func (d *Disk) ObjectBound() *pbrt.Bounds3 {
	return &pbrt.Bounds3{
		Min: &pbrt.Point3f{
			X: -d.radius,
			Y: -d.radius,
			Z: d.height,
		},
		Max: &pbrt.Point3f{
			X: d.radius,
			Y: d.radius,
			Z: d.height,
		},
	}
}
func (d *Disk) WorldBound() *pbrt.Bounds3 {
	return d.shape.objectToWorld.TransformBounds(d.ObjectBound())
}
func (d *Disk) ReverseOrientation() bool {
	return d.shape.reverseOrientation
}
func (d *Disk) TransformSwapsHandedness() bool {
	return d.shape.transformSwapsHandedness
}
func (d *Disk) Intersect(r *pbrt.Ray, si *pbrt.SurfaceInteraction, testAlphaTexture bool) (intersects bool, tHit float64) {
	// transform ray to object space
	ray, _, _ := d.shape.worldToObject.TransformRay(r)

	// compute plane intersection for disk

	// reject disk intersections for rays parallel to the disk's plane
	if ray.Direction.Z == 0 {
		return false, 0
	}
	tShapeHit := (d.height - ray.Origin.Z) / ray.Direction.Z
	if tShapeHit <= 0 || tShapeHit >= ray.TMax {
		return false, 0
	}

	// see if hit point is inside disk radii and phimax
	pHit := ray.PointAt(tShapeHit)
	dist2 := pHit.X*pHit.X + pHit.Y*pHit.Y
	if dist2 > d.radius*d.radius || dist2 < d.innerRadius*d.innerRadius {
		return false, 0
	}

	// test disk phi value against phiMax
	phi := math.Atan2(pHit.Y, pHit.X)
	if phi < 0 {
		phi += 2 * math.Pi
	}
	if phi > d.phiMax {
		return false, 0
	}

	// find parametric representation of disk hit
	u := phi / d.phiMax
	rHit := math.Sqrt(dist2)
	oneMinusV := (rHit - d.innerRadius) / (d.radius - d.innerRadius)
	v := 1 - oneMinusV
	dpdu := &pbrt.Vector3f{-d.phiMax * pHit.Y, d.phiMax * pHit.X, 0}
	dpdv := &pbrt.Vector3f{pHit.X, pHit.Y, 0}
	dpdv = dpdv.MulScalar((d.radius - d.innerRadius) / rHit)
	dndu := new(pbrt.Vector3f)
	dndv := new(pbrt.Vector3f)

	// refine disk intersection point
	pHit.Z = d.height

	// compute error bounds for disk intersection
	*si = *d.shape.objectToWorld.TransformSurfaceInteraction(pbrt.NewSurfaceInteractionWith(
		pHit,
		new(pbrt.Vector3f),
		&pbrt.Point2f{u, v},
		ray.Direction.MulScalar(-1),
		dpdu,
		dpdv,
		dndu,
		dndv,
		ray.Time,
		d,
		0,
	))

	return true, tShapeHit

}
func (d *Disk) IntersectP(r *pbrt.Ray, testAlphaTexture bool) (intersects bool) {
	// transform ray to object space
	ray, _, _ := d.shape.worldToObject.TransformRay(r)

	// compute plane intersection for disk

	// reject disk intersections for rays parallel to the disk's plane
	if ray.Direction.Z == 0 {
		return false
	}
	tShapeHit := (d.height - ray.Origin.Z) / ray.Direction.Z
	if tShapeHit <= 0 || tShapeHit >= ray.TMax {
		return false
	}

	// see if hit point is inside disk radii and phimax
	pHit := ray.PointAt(tShapeHit)
	dist2 := pHit.X*pHit.X + pHit.Y*pHit.Y
	if dist2 > d.radius*d.radius || dist2 < d.innerRadius*d.innerRadius {
		return false
	}

	// test disk phi value against phiMax
	phi := math.Atan2(pHit.Y, pHit.X)
	if phi < 0 {
		phi += 2 * math.Pi
	}
	if phi > d.phiMax {
		return false
	}

	return true
}
func (d *Disk) Sample(u *pbrt.Point2f) (i pbrt.Interaction, pdf float64) {
	pd := pbrt.ConcentricSampleDisk(u)
	pObj := &pbrt.Point3f{pd.X * d.radius, pd.Y * d.radius, d.height}
	it := pbrt.NewSurfaceInteraction()
	it.Normal = d.shape.objectToWorld.TransformNormal(&pbrt.Normal3f{0, 0, 1})
	if d.shape.reverseOrientation {
		it.Normal = it.Normal.MulScalar(-1)
	}
	it.Point, it.PointError = d.shape.objectToWorld.TransformPoint(pObj, new(pbrt.Vector3f))
	return it, 1 / d.Area()
}
func (d *Disk) SampleAtInteraction(ref pbrt.Interaction, u *pbrt.Point2f) (i pbrt.Interaction, pdf float64) {
	return pbrt.SampleAtInteraction(d, ref, u)
}
func (d *Disk) Pdf(ref pbrt.Interaction) float64 {
	return pbrt.Pdf(d, ref)
}
func (d *Disk) PdfWi(ref pbrt.Interaction, wi *pbrt.Vector3f) float64 {
	return pbrt.PdfWi(d, ref, wi)
}
func (d *Disk) SolidAngle(p *pbrt.Point3f, nSamples int) float64 {
	return pbrt.SolidAngle(d, p, nSamples)
}
func (d *Disk) Area() float64 {
	return d.phiMax * 0.5 * (d.radius*d.radius - d.innerRadius*d.innerRadius)
}
