//go:generate mockgen -source=reflection.go -destination=reflection.mock.go -package=pbrt

package pbrt

import "github.com/stupschwartz/go-pbrt/pkg/math"

const MaxBxDFs = 8

type BxDFType int

const (
	BSDFReflection BxDFType = 1 << iota
	BSDFTransmission
	BSDFDiffuse
	BSDFGlossy
	BSDFSpecular
	BSDFAll = BSDFDiffuse | BSDFGlossy | BSDFSpecular | BSDFReflection | BSDFTransmission
)

func SameHemisphere(w, wp *Vector3f) bool {
	return w.Z*wp.Z > 0
}

func CosTheta(w *Vector3f) float64 {
	return w.Z
}

func Cos2Theta(w *Vector3f) float64 {
	return w.Z * w.Z
}

func AbsCosTheta(w *Vector3f) float64 {
	return math.Abs(w.Z)
}

func Sin2Theta(w *Vector3f) float64 {
	return math.Max(0, 1-Cos2Theta(w))
}

func SinTheta(w *Vector3f) float64 {
	return math.Sqrt(Sin2Theta(w))
}

func TanTheta(w *Vector3f) float64 {
	return SinTheta(w) / CosTheta(w)
}

func Tan2Theta(w *Vector3f) float64 {
	return Sin2Theta(w) / Cos2Theta(w)
}

func CosPhi(w *Vector3f) float64 {
	sinTheta := SinTheta(w)
	if sinTheta == 0 {
		return 1
	}

	return math.Clamp(w.X/sinTheta, -1, 1)
}

func SinPhi(w *Vector3f) float64 {
	sinTheta := SinTheta(w)
	if sinTheta == 0 {
		return 0
	}

	return math.Clamp(w.Y/sinTheta, -1, 1)
}

func Cos2Phi(w *Vector3f) float64 {
	return CosPhi(w) * CosPhi(w)
}

func Sin2Phi(w *Vector3f) float64 {
	return SinPhi(w) * SinPhi(w)
}

type BSDF struct {
	Eta    float64
	ns, ng *Normal3f
	ss, ts *Vector3f
	nBxDFs int
	bxdfs  []BxDF
}

func NewBSDF(si *SurfaceInteraction, eta float64) *BSDF {
	ns := si.shading.normal
	ss := si.shading.dpdu.Normalized()
	return &BSDF{
		Eta:    eta,
		ns:     ns,
		ng:     si.Normal,
		ss:     ss,
		ts:     ns.Cross(ss),
		nBxDFs: 0,
		bxdfs:  []BxDF{},
	}
}

func (b *BSDF) Add(b2 BxDF) {
	b.bxdfs = append(b.bxdfs, b2)
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
		if MatchesFlags(b.bxdfs[i].GetType(), flags) {
			num++
		}
	}
	return num
}

func (b *BSDF) F(woW, wiW *Vector3f, flags BxDFType) Spectrum {
	wi := b.WorldToLocal(wiW)
	wo := b.WorldToLocal(woW)
	if wo.Z == 0.0 {
		return NewSpectrum(0)
	}
	reflect := wiW.Dot(b.ng)*woW.Dot(b.ng) > 0

	f := NewSpectrum(0)
	for i := 0; i < b.nBxDFs; i++ {
		if MatchesFlags(b.bxdfs[i].GetType(), flags) &&
			((reflect && (b.bxdfs[i].GetType()&BSDFReflection > 0)) ||
				(!reflect && (b.bxdfs[i].GetType()&BSDFTransmission > 0))) {
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
	var bxdf BxDF
	for i := 0; i < b.nBxDFs; i++ {
		if MatchesFlags(b.bxdfs[i].GetType(), t) {
			count--
			if count == 0 {
				bxdf = b.bxdfs[i]
				break
			}
		}
	}

	// remap sample u to [0,1)^2
	uRemapped := &Point2f{X: math.Min(u.X*float64(matchingComps)-comp, math.OneMinusEpsilon), Y: u.Y}

	// sample chosen bxDF
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
			if b.bxdfs[i] != bxdf && MatchesFlags(b.bxdfs[i].GetType(), t) {
				pdf += b.bxdfs[i].Pdf(wo, wi)
			}
		}
	}
	if matchingComps > 1 {
		pdf /= float64(matchingComps)
	}

	// compute value of BSDF for sampled Direction
	if (bxdf.GetType()&BSDFSpecular == 0) && matchingComps > 1 {
		reflect := wiWorld.Dot(b.ng)*woWorld.Dot(b.ng) > 0
		f.SetAll(0.0)
		for i := 0; i < b.nBxDFs; i++ {
			if MatchesFlags(b.bxdfs[i].GetType(), t) &&
				((reflect && (b.bxdfs[i].GetType()&BSDFReflection > 0)) ||
					(!reflect && (b.bxdfs[i].GetType()&BSDFTransmission > 0))) {
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
		if MatchesFlags(b.bxdfs[i].GetType(), flags) {
			matchingComps++
			pdf += b.bxdfs[i].Pdf(wo, wi)
		}
	}
	if matchingComps <= 0 {
		return 0
	}
	return pdf / float64(matchingComps)
}

type BxDF interface {
	GetType() BxDFType
	F(woW, wiW *Vector3f) Spectrum
	SampleF(wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType)
	Rho(wo *Vector3f, samples []*Point2f) Spectrum
	RhoSamples(samples1, samples2 []*Point2f) Spectrum
	Pdf(wo, wi *Vector3f) float64
}

type bxDF struct {
	Type BxDFType
}

func NewBxDF(t BxDFType) *bxDF {
	return &bxDF{Type: t}
}

func (b *bxDF) GetType() BxDFType {
	return b.Type
}

func MatchesFlags(t, flags BxDFType) bool {
	return (t & flags) == t
}

//func (b *bxDF) F(wo, wi *Vector3f) Spectrum {
//	return nil
//}

func sampleF(b BxDF, wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType) {
	// cosine sample the hemisphere, flipping the Direction if necessary
	wi = CosineSampleHemisphere(sample)
	if wo.Z < 0 {
		wi.Z *= -1
	}

	pdf = b.Pdf(wo, wi)
	return b.F(wo, wi), wi, pdf, 0
}

func rho(b BxDF, w *Vector3f, samples []*Point2f) Spectrum {
	r := NewSpectrum(0)
	nSamples := len(samples)
	for i := 0; i < nSamples; i++ {
		f, wi, pdf, _ := sampleF(b, w, samples[i])
		if pdf > 0 {
			r.AddAssign(f.MulScalar(AbsCosTheta(wi) / pdf))
		}
	}
	return r.DivScalar(float64(nSamples))
}
func rhoSamples(b BxDF, samples1, samples2 []*Point2f) Spectrum {
	r := NewSpectrum(0)
	nSamples := len(samples1)
	for i := 0; i < nSamples; i++ {
		// estimate one term of rho
		wo := UniformSampleHemisphere(samples2[i])
		pdfo := UniformHemispherePdf()
		f, wi, pdfi, _ := sampleF(b, wo, samples2[i])
		if pdfi > 0 {
			r.AddAssign(f.MulScalar(AbsCosTheta(wi) * AbsCosTheta(wo) / (pdfo * pdfi)))
		}
	}

	return r.DivScalar(math.Pi * float64(nSamples))
}

func pdf(wo, wi *Vector3f) float64 {
	if SameHemisphere(wo, wi) {
		return AbsCosTheta(wi) * math.InvPi
	}
	return 0
}

type ScaledBxDF struct {
	*bxDF

	bxdf  BxDF
	scale Spectrum
}

func NewScaledBxDF(bxdf BxDF, scale Spectrum) *ScaledBxDF {
	return &ScaledBxDF{
		bxDF:  NewBxDF(bxdf.GetType()),
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

type Fresnel interface {
	Evaluate(cosI float64) Spectrum
}

type FresnelNoOp struct {
}

func (f *FresnelNoOp) Evaluate(cosI float64) Spectrum {
	return NewSpectrum(1.0)
}

func NewSpecularReflection(r Spectrum, fresnel Fresnel) *SpecularReflection {
	return &SpecularReflection{
		bxDF:    NewBxDF(BSDFReflection | BSDFDiffuse),
		r:       r,
		fresnel: fresnel,
	}
}

type SpecularReflection struct {
	*bxDF

	r       Spectrum
	fresnel Fresnel
}

func (l *SpecularReflection) F(woW, wiW *Vector3f) Spectrum {
	return NewSpectrum(0.0)
}

func (l *SpecularReflection) SampleF(wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType) {
	// Compute perfect specular reflection Direction
	wi = &Vector3f{X: -wo.X, Y: -wo.Y, Z: wo.Z}
	pdf = 1.0
	return l.fresnel.Evaluate(CosTheta(wi)).Mul(l.r).DivScalar(AbsCosTheta(wi)), wi, pdf, 0
}

func (l *SpecularReflection) Rho(wo *Vector3f, samples []*Point2f) Spectrum {
	return l.r
}

func (l *SpecularReflection) RhoSamples(samples1, samples2 []*Point2f) Spectrum {
	return l.r
}

func (l *SpecularReflection) Pdf(wo, wi *Vector3f) float64 {
	return 0.0
}

func NewLambertianReflection(r Spectrum) *LambertianReflection {
	return &LambertianReflection{
		bxDF: NewBxDF(BSDFReflection | BSDFDiffuse),
		r:    r,
	}
}

type LambertianReflection struct {
	*bxDF

	r Spectrum
}

func (l *LambertianReflection) F(woW, wiW *Vector3f) Spectrum {
	return l.r.MulScalar(math.InvPi)
}

func (l *LambertianReflection) SampleF(wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType) {
	return sampleF(l, wo, sample)
}

func (l *LambertianReflection) Rho(wo *Vector3f, samples []*Point2f) Spectrum {
	return l.r
}

func (l *LambertianReflection) RhoSamples(samples1, samples2 []*Point2f) Spectrum {
	return l.r
}

func (l *LambertianReflection) Pdf(wo, wi *Vector3f) float64 {
	return pdf(wo, wi)
}

type OrenNayar struct {
	*bxDF

	r    Spectrum
	a, b float64
}

func NewOrenNayar(r Spectrum, sigma float64) *OrenNayar {
	sigma = math.Radians(sigma)
	sigma2 := sigma * sigma
	return &OrenNayar{
		bxDF: NewBxDF(BSDFReflection | BSDFDiffuse),
		r:    r,
		a:    1.0 - (sigma2 / (2.0 * (sigma2 + 0.33))),
		b:    0.45 * sigma2 / (sigma2 * 0.09),
	}
}

func (o *OrenNayar) F(wo, wi *Vector3f) Spectrum {
	sinThetaI := SinTheta(wi)
	sinThetaO := SinTheta(wo)

	// compute cosine term of oren-nayar model
	maxCos := 0.0
	if sinThetaI > 1e-4 && sinThetaO > 1e-4 {
		sinPhiI := SinPhi(wi)
		cosPhiI := CosPhi(wi)
		sinPhiO := SinPhi(wo)
		cosPhiO := CosPhi(wo)
		dCos := cosPhiI*cosPhiO + sinPhiI*sinPhiO
		maxCos = math.Max(0.0, dCos)
	}

	// compute sine and tangent terms of Oren-Nayar model
	var sinAlpha, tanBeta float64
	if AbsCosTheta(wi) > AbsCosTheta(wo) {
		sinAlpha = sinThetaO
		tanBeta = sinThetaO / AbsCosTheta(wo)
	} else {
		sinAlpha = sinThetaI
		tanBeta = sinThetaO / AbsCosTheta(wo)
	}
	return o.r.MulScalar(math.InvPi * (o.a + o.b*maxCos*sinAlpha*tanBeta))
}

func (o *OrenNayar) SampleF(wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType) {
	return sampleF(o, wo, sample)
}

func (o *OrenNayar) Rho(wo *Vector3f, samples []*Point2f) Spectrum {
	return rho(o, wo, samples)
}

func (o *OrenNayar) RhoSamples(samples1, samples2 []*Point2f) Spectrum {
	return rhoSamples(o, samples1, samples2)
}

func (o *OrenNayar) Pdf(wo, wi *Vector3f) float64 {
	return pdf(wo, wi)
}
