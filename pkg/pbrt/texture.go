//go:generate mockgen -source=texture.go -destination=texture.mock.go -package=pbrt

package pbrt

type TextureMapping2D interface {
	Map(si *SurfaceInteraction) (p *Point2f, dstdx, dstdy *Vector2f)
}

func NewUvMapping2D() *UVMapping2D {
	return &UVMapping2D{
		su: 1,
		sv: 1,
		du: 0,
		dv: 0,
	}
}

type UVMapping2D struct {
	su, sv, du, dv float64
}

func (t *UVMapping2D) Map(si *SurfaceInteraction) (p *Point2f, dstdx, dstdy *Vector2f) {
	return &Point2f{t.su*si.uv.X + t.du, t.sv*si.uv.Y + t.dv},
		&Vector2f{t.su * si.dudx, t.sv * si.dvdx},
		&Vector2f{t.su * si.dudy, t.sv * si.dvdy}
}

func NewPlanarMapping2D(vs, vt *Vector3f, ds, dt float64) *PlanarMapping2D {
	return &PlanarMapping2D{
		vs: vs,
		vt: vt,
		ds: ds,
		dt: dt,
	}
}

type PlanarMapping2D struct {
	vs, vt *Vector3f
	ds, dt float64
}

func (m *PlanarMapping2D) Map(si *SurfaceInteraction) (p *Point2f, dstdx, dstdy *Vector2f) {
	return &Point2f{m.ds + si.Point.Dot(m.vs), m.dt + si.Point.Dot(m.vt)},
		&Vector2f{si.dpdx.Dot(m.vs), si.dpdx.Dot(m.vt)},
		&Vector2f{si.dpdy.Dot(m.vs), si.dpdy.Dot(m.vt)}
}

type SpectrumTexture interface {
	Evaluate(si *SurfaceInteraction) Spectrum
}

type FloatTexture interface {
	Evaluate(si *SurfaceInteraction) float64
}

type ConstantSpectrumTexture struct {
	value Spectrum
}

func NewConstantSpectrumTexture(s Spectrum) *ConstantSpectrumTexture {
	return &ConstantSpectrumTexture{
		value: s,
	}
}

func (t *ConstantSpectrumTexture) Evaluate(si *SurfaceInteraction) Spectrum {
	return t.value
}

type ConstantFloatTexture struct {
	value float64
}

func NewConstantFloatTexture(v float64) FloatTexture {
	return &ConstantFloatTexture{
		value: v,
	}
}

func (t *ConstantFloatTexture) Evaluate(si *SurfaceInteraction) float64 {
	return t.value
}
