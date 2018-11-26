//go:generate mockgen -source=interaction.go -destination=interaction.mock.go -package=pbrt

package pbrt

import "github.com/stupschwartz/go-pbrt/pkg/math"

type Interaction interface {
	//IsSurfaceInteraction() bool
	//IsMediumInteraction() bool
	//
	SpawnRay(direction *Vector3f) *Ray
	SpawnRayToPoint(p *Point3f) *Ray
	SpawnRayToInteraction(to Interaction) *Ray

	GetPoint() *Point3f
	GetPointError() *Vector3f
	GetTime() float64
	GetNormal() *Normal3f
	//GetMedium(*Vector3f) *Medium
	SetMediumAccessor(accessor *MediumAccessor)
}

type interaction struct {
	Point          *Point3f
	PointError     *Vector3f
	time           float64
	wo             *Vector3f
	Normal         *Normal3f
	mediumAccessor *MediumAccessor
}

func NewInteraction(p *Point3f, pError *Vector3f, n *Normal3f, wo *Vector3f, time float64, mediumAccessor *MediumAccessor) *interaction {
	return &interaction{
		Point:          p,
		PointError:     pError,
		Normal:         n,
		wo:             wo,
		time:           time,
		mediumAccessor: mediumAccessor,
	}
}

func (i *interaction) GetPoint() *Point3f {
	return i.Point
}

func (i *interaction) GetPointError() *Vector3f {
	return i.PointError
}

func (i *interaction) GetNormal() *Normal3f {
	return i.Normal
}

func (i *interaction) GetTime() float64 {
	return i.time
}

func (i *interaction) SetMediumAccessor(accessor *MediumAccessor) {
	i.mediumAccessor = accessor
}

func (i *interaction) SpawnRay(direction *Vector3f) *Ray {
	origin := OffsetRayOrigin(i.Point, i.PointError, i.Normal, direction)
	return &Ray{origin, direction, math.Infinity, i.time, i.GetMedium(direction)}
}

func (i *interaction) SpawnRayToPoint(p *Point3f) *Ray {
	direction := p.Sub(i.Point)
	origin := OffsetRayOrigin(i.Point, i.PointError, i.Normal, direction)
	return &Ray{origin, direction, 1 - math.ShadowEpsilon, i.time, i.GetMedium(direction)}
}

func (i *interaction) SpawnRayToInteraction(to Interaction) *Ray {
	origin := OffsetRayOrigin(i.Point, i.PointError, i.Normal, to.GetPoint().Sub(i.Point))
	target := OffsetRayOrigin(to.GetPoint(), to.GetPointError(), to.GetNormal(), origin.Sub(to.GetPoint()))
	direction := target.Sub(origin)
	return &Ray{
		Origin:    i.Point,
		Direction: direction,
		TMax:      1 - math.ShadowEpsilon,
		Time:      i.time,
		Medium:    i.GetMedium(direction),
	}
}

func (i *interaction) IsSurfaceInteraction() bool {
	return i.Normal != nil
}

func (i *interaction) IsMediumInteraction() bool {
	return !i.IsSurfaceInteraction()
}

func (i *interaction) GetMedium(w *Vector3f) Medium {
	if w.Dot(i.Normal) > 0 {
		return i.mediumAccessor.Outside
	}
	return i.mediumAccessor.Inside
}

func PhaseHG(cosTheta float64, g float64) float64 {
	denom := 1.0 + g*g + 2*g*cosTheta
	return math.Inv4Pi * (1.0 - g*g) / (denom * math.Sqrt(denom))
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
	shape                  Shape
	shading                *Shading
	Primitive              Primitive
	BSDF                   *BSDF
	bssrdf                 *BSSRDF
	dpdx, dpdy             *Vector3f
	dudx, dvdx, dudy, dvdy float64

	// Shapes can optionally provide a face index with an
	// intersection Point for use in Ptex texture lookups.
	// If Ptex isn'Type being used, then this value is ignored.
	faceIndex int
}

func NewSurfaceInteraction() *SurfaceInteraction {
	return &SurfaceInteraction{
		interaction: &interaction{
			Point:          new(Point3f),
			wo:             new(Vector3f),
			Normal:         new(Normal3f),
			mediumAccessor: nil,
		},
		uv:    new(Point2f),
		dpdu:  new(Vector3f),
		dpdv:  new(Vector3f),
		dndu:  new(Normal3f),
		dndv:  new(Normal3f),
		shape: nil,
		shading: &Shading{
			normal: new(Normal3f),
			dpdu:   new(Vector3f),
			dpdv:   new(Vector3f),
			dndu:   new(Normal3f),
			dndv:   new(Normal3f),
		},
		Primitive: nil,
		faceIndex: 0,
	}
}

func NewSurfaceInteractionWith(p *Point3f, pError *Vector3f, uv *Point2f, wo *Vector3f, dpdu, dpdv *Vector3f, dndu, dndv *Normal3f, time float64, shape Shape, faceIndex int) *SurfaceInteraction {
	normal := dpdu.Cross(dpdv).Normalized()

	return &SurfaceInteraction{
		interaction: &interaction{
			Point:      p,
			PointError: pError,
			time:       time,
			wo:         wo,
			Normal:     normal,
		},
		//bssrdf: new(BSSRDF),
		uv:    uv,
		dpdu:  dpdu,
		dpdv:  dpdv,
		dndu:  dndu,
		dndv:  dndv,
		shape: shape,
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
	area := si.Primitive.GetAreaLight()
	if area != nil {
		return area.L(si, w)
	}
	return NewSpectrum(0.0)
}

func (si *SurfaceInteraction) ComputeScatteringFunctions(ray *RayDifferential, allowMultipleLobes bool, mode TransportMode) {
	si.ComputeDifferentials(ray)
	if si.Primitive == nil {
		return
	}
	si.Primitive.ComputeScatteringFunctions(si, mode, allowMultipleLobes)
}

func (si *SurfaceInteraction) ComputeDifferentials(ray *RayDifferential) {
	if ray.hasDifferentials {
		// estimate screen space change in pt and (u,v)

		// compute auxiliary intersection points with plane
		if si.interaction == nil {
			return
		}
		d := si.Normal.Dot(si.Point)
		tx := -(si.Normal.Dot(ray.rxOrigin) - d) / si.Normal.Dot(ray.rxDirection)
		if math.IsNaN(tx) {
			//goto Failed
			return
		}
		px := ray.rxOrigin.Add(ray.rxDirection.MulScalar(tx))

		ty := -(si.Normal.Dot(ray.ryOrigin) - d) / si.Normal.Dot(ray.ryDirection)
		if math.IsNaN(tx) {
			//goto Failed
			return
		}
		py := ray.ryOrigin.Add(ray.ryDirection.MulScalar(ty))

		si.dpdx = px.Sub(si.Point)
		si.dpdy = py.Sub(si.Point)

		// compute (u,v) offsets at auxiliary points

		// choose two dimesnions to use for ray offset computation
		var dim [2]int
		if math.Abs(si.Normal.X) > math.Abs(si.Normal.Y) && math.Abs(si.Normal.X) > math.Abs(si.Normal.Z) {
			dim[0] = 1
			dim[1] = 2
		} else if math.Abs(si.Normal.Y) > math.Abs(si.Normal.Z) {
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
		Bx := [2]float64{px.Index(dim[0]) - si.Point.Index(dim[0]), px.Index(dim[1]) - si.Point.Index(dim[1])}
		By := [2]float64{py.Index(dim[0]) - si.Point.Index(dim[0]), py.Index(dim[1]) - si.Point.Index(dim[1])}

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
