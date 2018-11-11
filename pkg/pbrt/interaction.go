package pbrt

import (
	"math"
)

type Interaction interface {
	//IsSurfaceInteraction() bool
	//IsMediumInteraction() bool
	//
	SpawnRay(direction *Vector3f) *Ray
	SpawnRayToPoint(p *Point3f) *Ray
	SpawnRayToInteraction(to Interaction) *Ray

	GetPoint() *Point3f
	//GetPointError() *Vector3f
	GetTime() float64
	GetNormal() *Normal3f
	//GetMedium(*Vector3f) *Mediumer
	SetMediumAccessor(accessor *MediumAccessor)
}

type interaction struct {
	point          *Point3f
	time           float64
	wo             *Vector3f
	normal         *Normal3f
	mediumAccessor *MediumAccessor
}

func NewInteraction(p *Point3f, n *Normal3f, wo *Vector3f, time float64, mediumAccessor *MediumAccessor) *interaction {
	return &interaction{
		point:          p,
		normal:         n,
		wo:             wo,
		time:           time,
		mediumAccessor: mediumAccessor,
	}
}

func (i *interaction) GetPoint() *Point3f {
	return i.point
}

func (i *interaction) GetNormal() *Normal3f {
	return i.normal
}

func (i *interaction) GetTime() float64 {
	return i.time
}

func (i *interaction) SetMediumAccessor(accessor *MediumAccessor) {
	i.mediumAccessor = accessor
}

func (i *interaction) SpawnRay(direction *Vector3f) *Ray {
	return &Ray{i.point, direction, Infinity, i.time, i.GetMedium(direction)}
}

func (i *interaction) SpawnRayToPoint(p *Point3f) *Ray {
	direction := p.Sub(i.point)
	return &Ray{i.point, direction, 1 - ShadowEpsilon, i.time, i.GetMedium(direction)}
}

func (i *interaction) SpawnRayToInteraction(to Interaction) *Ray {
	direction := to.GetPoint().Sub(i.point)
	return &Ray{
		origin:    i.point,
		direction: direction,
		tMax:      1 - ShadowEpsilon,
		time:      i.time,
		medium:    i.GetMedium(direction),
	}
}

func (i *interaction) IsSurfaceInteraction() bool {
	return i.normal != nil
}

func (i *interaction) IsMediumInteraction() bool {
	return !i.IsSurfaceInteraction()
}

func (i *interaction) GetMedium(w *Vector3f) Mediumer {
	if w.Dot(i.normal) > 0 {
		return i.mediumAccessor.Outside
	}
	return i.mediumAccessor.Inside
}

func PhaseHG(cosTheta float64, g float64) float64 {
	denom := 1.0 + g*g + 2*g*cosTheta
	return Inv4Pi * (1.0 - g*g) / (denom * math.Sqrt(denom))
}

type Shading struct {
	normal     *Normal3f
	dpdu, dpdv *Vector3f
	dndu, dndv *Normal3f
}

type SurfaceInteraction struct {
	*interaction

	uv                     *Point2f
	dpdu, dpdv             *Vector3f
	dndu, dndv             *Normal3f
	shape                  Shaper
	shading                *Shading
	primitive              Primitiver
	bsdf                   *BSDF
	bssrdf                 *BSSRDF
	dpdx, dpdy             *Vector3f
	dudx, dvdx, dudy, dvdy float64

	// Shapes can optionally provide a face index with an
	// intersection point for use in Ptex texture lookups.
	// If Ptex isn'Type being used, then this value is ignored.
	faceIndex int
}

func NewSurfaceInteraction(p *Point3f, pError *Vector3f, uv *Point2f, wo *Vector3f, dpdu, dpdv *Vector3f, dndu, dndv *Normal3f, time float64, shape Shaper, faceIndex int) *SurfaceInteraction {
	normal := dpdu.Cross(dpdv).Normalized()

	return &SurfaceInteraction{
		interaction: &interaction{
			point:  p,
			time:   time,
			wo:     wo,
			normal: normal,
		},
		//bssrdf: new(BSSRDF),
		uv:     uv,
		dpdu:   dpdu,
		dpdv:   dpdv,
		dndu:   dndu,
		dndv:   dndv,
		shape:  shape,
		shading: &Shading{
			normal: normal,
			dpdu:   dpdu,
			dpdv:   dpdv,
			dndu:   dndu,
			dndv:   dndv,
		},
		faceIndex: faceIndex,
	}
}

func (si *SurfaceInteraction) Le(w *Vector3f) Spectrum {
	area := si.primitive.GetAreaLight()
	if area != nil {
		return area.L(si, w)
	}
	return NewSpectrum(0.0)
}

func (si *SurfaceInteraction) ComputeScatteringFunctions(ray *RayDifferential, allowMultipleLobes bool, mode TransportMode) {
	si.ComputeDifferentials(ray)
	si.primitive.ComputeScatteringFunctions(si, mode, allowMultipleLobes)
}

func (si *SurfaceInteraction) ComputeDifferentials(ray *RayDifferential) {
	if ray.hasDifferentials {
		// estimate screen space change in pt and (u,v)

		// compute auxiliary intersection points with plane
		d := si.normal.Dot(si.point)
		tx := -(si.normal.Dot(ray.rxOrigin) - d) / si.normal.Dot(ray.rxDirection)
		if math.IsNaN(tx) {
			//goto Failed
			return
		}
		px := ray.rxOrigin.Add(ray.rxDirection.MulScalar(tx))

		ty := -(si.normal.Dot(ray.ryOrigin) - d) / si.normal.Dot(ray.ryDirection)
		if math.IsNaN(tx) {
			//goto Failed
			return
		}
		py := ray.ryOrigin.Add(ray.ryDirection.MulScalar(ty))

		si.dpdx = px.Sub(si.point)
		si.dpdy = py.Sub(si.point)

		// compute (u,v) offsets at auxiliary points

		// choose two dimesnions to use for ray offset computation
		var dim [2]int
		if math.Abs(si.normal.X) > math.Abs(si.normal.Y) && math.Abs(si.normal.X) > math.Abs(si.normal.Z) {
			dim[0] = 1
			dim[1] = 2
		} else if math.Abs(si.normal.Y) > math.Abs(si.normal.Z) {
			dim[0] = 0
			dim[1] = 2
		} else {
			dim[0] = 0
			dim[1] = 1
		}

		// initialize A, Bx, and By matrices for offset computation
		A := [2][2]float64{
			{si.dpdu.Index(dim[0]), si.dpdv.Index(dim[0])},
			{si.dpdu.Index(dim[1]), si.dpdv.Index(dim[1])},
		}
		Bx := [2]float64{px.Index(dim[0]) - si.point.Index(dim[0]), px.Index(dim[1]) - si.point.Index(dim[1])}
		By := [2]float64{py.Index(dim[0]) - si.point.Index(dim[0]), py.Index(dim[1]) - si.point.Index(dim[1])}

		solved, dudx, dvdx := SolveLinearSystem2x2(A, Bx)
		if !solved {
			si.dudx = 0
			si.dvdx = 0
		}
		si.dudx = dudx
		si.dvdx = dvdx

		solved, dudy, dvdy := SolveLinearSystem2x2(A, By)
		if !solved {
			si.dudx = 0
			si.dvdx = 0
		}
		si.dudy = dudy
		si.dvdy = dvdy
	} else {
		si.dudx = 0
		si.dvdx = 0
		si.dudy = 0
		si.dvdy = 0
		si.dpdx = &Vector3f{}
		si.dpdy = &Vector3f{}
	}
}

type MediumInteraction struct {
	*interaction

	phase PhaseFunction
}

func (mi *MediumInteraction) IsValid() bool {
	return mi.phase != nil
}

type HenyeyGreenstein struct {
	g float64
}

func (hg *HenyeyGreenstein) P(wo, wi *Vector3f) float64 {
	return PhaseHG(wo.Dot(wi), hg.g)
}

func (hg *HenyeyGreenstein) SampleP(wo, wi *Vector3f, u *Point2f) float64 {
	var cosTheta float64
	if math.Abs(hg.g) < 1e-3 {
		cosTheta = 1.0 - 2.0*u.X
	} else {
		sqrTerm := (1.0 - hg.g*hg.g) / (1.0 - hg.g + 2*hg.g*u.X)
		cosTheta = (1.0 + hg.g*hg.g - sqrTerm*sqrTerm) / (2.0 * hg.g)
	}

	sinTheta := math.Sqrt(math.Max(0.0, 1.0-cosTheta*cosTheta))
	phi := 2.0 * math.Pi * u.Y
	v1, v2 := CoordinateSystem(wo)
	*wi = *SphericalDirectionXYZ(sinTheta, cosTheta, phi, v1, v2, wo.MulScalar(-1.0))
	return PhaseHG(-cosTheta, hg.g)
}
