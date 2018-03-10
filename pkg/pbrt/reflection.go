package pbrt

import "math"

const MaxBxDFs = 8

type BxDFType int

const (
	BSDFReflection   BxDFType = 1 << iota
	BSDFTransmission
	BSDFDiffuse
	BSDFGlossy
	BSDFSpecular
	BSDFAll          = BSDFDiffuse | BSDFGlossy | BSDFSpecular | BSDFReflection | BSDFTransmission
)

func SameHemisphere(w, wp *Vector3f) bool {
	return w.Z*wp.Z > 0
}

func AbsCosTheta(w *Vector3f) float64 {
	return math.Abs(w.Z)
}

type BSDF struct {
	Eta    float64
	ns, ng *Normal3f
	ss, ts *Vector3f
	nBxDFs int
	bxdfs  [MaxBxDFs]BxDFer
}

func (b *BSDF) Add(b2 BxDFer) {
	b.bxdfs[b.nBxDFs] = b2
	b.nBxDFs++
}

func (b *BSDF) WorldToLocal(v *Vector3f) *Vector3f {
	return &Vector3f{v.Dot(b.ss), v.Dot(b.ts), v.Dot(b.ns)}
}

func (b *BSDF) LocalToWorld(v *Vector3f) *Vector3f {
	return &Vector3f{
		X: b.ss.X*v.X + b.ts.X*v.Y + b.ns.X*v.Z,
		Y: b.ss.Y*v.X + b.ts.Y*v.Y + b.ns.Y*v.Z,
		Z: b.ss.Z*v.X + b.ts.Z*v.Y + b.ns.Z*v.Z,
	}
}

func (b *BSDF) NumComponents(flags BxDFType) int {
	num := 0
	for i := 0; i < b.nBxDFs; i++ {
		if b.bxdfs[i].MatchesFlags(flags) {
			num++
		}
	}
	return num
}

func (b *BSDF) F(woW, wiW *Vector3f, flags BxDFType) Spectrum {
	wi := b.WorldToLocal(wiW)
	wo := b.WorldToLocal(wiW)
	if wo.Z == 0.0 {
		return NewSpectrum(0)
	}
	reflect := wiW.Dot(b.ng)*woW.Dot(b.ng) > 0

	f := NewSpectrum(0)
	for i := 0; i < b.nBxDFs; i++ {
		if b.bxdfs[i].MatchesFlags(flags) &&
			((reflect && (b.bxdfs[i].GetType()&BSDFReflection == 1)) ||
				(!reflect && (b.bxdfs[i].GetType()&BSDFTransmission == 1))) {
			f.AddAssign(b.bxdfs[i].F(wo, wi))
		}
	}
	return f
}

func (b *BSDF) SampleF(woWorld *Vector3f, u *Point2f, t BxDFType) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType) {
	matchingComps := b.NumComponents(t)
	if matchingComps == 0 {
		return NewSpectrum(0), new(Vector3f), 0, BxDFType(0)
	}
	comp := math.Min(math.Floor(u.X*float64(matchingComps)), float64(matchingComps)-1)

	// get bxdf pointer for chosen component
	count := comp
	var bxdf BxDFer
	for i := 0; i < b.nBxDFs; i++ {
		if b.bxdfs[i].MatchesFlags(t) {
			count--
			if count == 0 {
				bxdf = b.bxdfs[i]
				break
			}
		}
	}

	// remap sample u to [0,1)^2
	uRemapped := &Point2f{u.X*float64(matchingComps) - comp, OneMinusEpsilon}

	// sample chosen BxDF
	wo := b.WorldToLocal(woWorld)
	if wo.Z == 0.0 {
		return NewSpectrum(0), new(Vector3f), 0, BxDFType(0)
	}

	f, wi, pdf, sampledType := bxdf.SampleF(wo, uRemapped)
	if pdf == 0.0 {
		return NewSpectrum(0), new(Vector3f), 0, BxDFType(0)
	}

	wiWorld := b.LocalToWorld(wi)

	// compute overall PDF with all matching BxDFs
	if (bxdf.GetType()&BSDFSpecular == 0) && matchingComps > 1 {
		for i := 0; i < b.nBxDFs; i++ {
			if b.bxdfs[i] != bxdf && b.bxdfs[i].MatchesFlags(t) {
				pdf += b.bxdfs[i].Pdf(wo, wi)
			}
		}
	}
	if matchingComps > 1 {
		pdf /= float64(matchingComps)
	}

	// compute value of BSDF for sampled direction
	if (bxdf.GetType()&BSDFSpecular == 0) && matchingComps > 1 {
		reflect := wiWorld.Dot(b.ng)*woWorld.Dot(b.ng) > 0
		f.SetAll(0.0)
		for i := 0; i < b.nBxDFs; i++ {
			if b.bxdfs[i].MatchesFlags(t) &&
				((reflect && (b.bxdfs[i].GetType()&BSDFReflection == 1)) ||
					(!reflect && (b.bxdfs[i].GetType()&BSDFTransmission == 1))) {
				f.AddAssign(b.bxdfs[i].F(wo, wi))
			}
		}
	}

	return f, wi, pdf, sampledType
}

func (b *BSDF) Pdf(woWorld, wiWorld *Vector3f, flags BxDFType) float64 {
	if b.nBxDFs == 0.0 {
		return 0
	}

	wo := b.WorldToLocal(woWorld)
	wi := b.WorldToLocal(wiWorld)
	if wo.Z == 0 {
		return 0
	}

	var pdf float64
	var matchingComps int
	for i := 0; i < b.nBxDFs; i++ {
		if b.bxdfs[i].MatchesFlags(flags) {
			matchingComps++
			pdf += b.bxdfs[i].Pdf(wo, wi)
		}
	}
	if matchingComps <= 0 {
		return 0
	}
	return pdf / float64(matchingComps)
}

type BxDFer interface {
	GetType() BxDFType
	MatchesFlags(t BxDFType) bool
	F(woW, wiW *Vector3f) Spectrum
	SampleF(wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType)
	Rho(wo *Vector3f, samples []*Point2f) Spectrum
	RhoSamples(samples1, samples2 []*Point2f) Spectrum
	Pdf(wo, wi *Vector3f) float64
}

type BxDF struct {
	Type BxDFType
}

func NewBxDF(t BxDFType) *BxDF {
	return &BxDF{Type: t}
}

func (b *BxDF) GetType() BxDFType {
	return b.Type
}

func (b *BxDF) MatchesFlags(t BxDFType) bool {
	return (b.Type & t) == b.Type
}

func (b *BxDF) F(wo, wi *Vector3f) Spectrum {
	return nil
}

func (b *BxDF) SampleF(wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType) {
	// cosine sample the hemisphere, flipping the direction if necessary
	wi = CosineSampleHemisphere(sample)
	if wo.Z < 0 {
		wi.Z *= -1
	}

	pdf = b.Pdf(wo, wi)
	return b.F(wo, wi), wi, pdf, 0
}

func (b *BxDF) Rho(w *Vector3f, samples []*Point2f) Spectrum {
	r := NewSpectrum(0)
	nSamples := len(samples)
	for i := 0; i < nSamples; i++ {
		f, wi, pdf, _ := b.SampleF(w, samples[i])
		if pdf > 0 {
			r.AddAssign(f.MulScalar(AbsCosTheta(wi) / pdf))
		}
	}
	return r.DivScalar(float64(nSamples))
}
func (b *BxDF) RhoSamples(samples1, samples2 []*Point2f) Spectrum {
	r := NewSpectrum(0)
	nSamples := len(samples1)
	for i := 0; i < nSamples; i++ {
		// estimate one term of rho
		wo := UniformSampleHemisphere(samples2[i])
		pdfo := UniformHemispherePdf()
		f, wi, pdfi, _ := b.SampleF(wo, samples2[i])
		if pdfi > 0 {
			r.AddAssign(f.MulScalar(AbsCosTheta(wi) * AbsCosTheta(wo) / (pdfo * pdfi)))
		}
	}

	return r.DivScalar(math.Pi * float64(nSamples))
}

func (b *BxDF) Pdf(wo, wi *Vector3f) float64 {
	if SameHemisphere(wo, wi) {
		return AbsCosTheta(wi) * InvPi
	}
	return 0
}

type ScaledBxDF struct {
	*BxDF

	bxdf  BxDFer
	scale Spectrum
}

func NewScaledBxDF(bxdf BxDFer, scale Spectrum) *ScaledBxDF {
	return &ScaledBxDF{
		BxDF:  NewBxDF(bxdf.GetType()),
		bxdf:  bxdf,
		scale: scale,
	}
}

func (b *ScaledBxDF) F(wo, wi *Vector3f) Spectrum {
	return b.scale.Mul(b.bxdf.F(wo, wi))
}

func (b *ScaledBxDF) SampleF(wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType) {
	f, wi, pdf, sampledType := b.bxdf.SampleF(wo, sample)
	return b.scale.Mul(f), wi, pdf, sampledType
}

func (b *ScaledBxDF) Pdf(wo, wi *Vector3f) float64 {
	return b.bxdf.Pdf(wo, wi)
}
