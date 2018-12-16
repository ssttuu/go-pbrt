package pbrt

import (
	"github.com/ssttuu/go-pbrt/pkg/math"
)

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

func FrDielectric(cosThetaI, etaI, etaT float64) float64 {
	cosThetaI = math.Clamp(cosThetaI, -1, 1)
	// potentially swap indices of refraction
	entering := cosThetaI > 0
	if !entering {
		etaI, etaT = etaT, etaI
		cosThetaI = math.Abs(cosThetaI)
	}

	// compute cosThetaT using Snell's law
	sinThetaI := math.Sqrt(math.Max(0, 1-cosThetaI*cosThetaI))
	sinThetaT := etaI / etaT * sinThetaI

	// handle total internal reflection
	if sinThetaT >= 1 {
		return 1
	}
	cosThetaT := math.Sqrt(math.Max(0, 1-sinThetaT*sinThetaT))
	Rparl := ((etaT * cosThetaI) - (etaI * cosThetaT)) / ((etaT * cosThetaI) + (etaI * cosThetaT))
	Rperp := ((etaI * cosThetaI) - (etaT * cosThetaT)) / ((etaI * cosThetaI) + (etaT * cosThetaT))
	return (Rparl*Rparl + Rperp*Rperp) / 2
}

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

func Reflect(wo, n *Vector3f) *Vector3f {
	return wo.MulScalar(-1).Add(n.MulScalar(2 * wo.Dot(n)))
}

func Refract(wi *Vector3f, n *Normal3f, eta float64) (*Vector3f, bool) {
	// compute cos theta_roman using Snell's law
	cosThetaI := n.Dot(wi)
	sin2ThetaI := math.Max(0, 1-cosThetaI*cosThetaI)
	sin2ThetaT := eta * eta * sin2ThetaI

	// handle total internal reflection for transmission
	if sin2ThetaT >= 1 {
		return nil, false
	}
	cosThetaT := math.Sqrt(1 - sin2ThetaT)
	return wi.MulScalar(-eta).Add(n.MulScalar(eta*cosThetaI - cosThetaT)), true
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
			if count == 0 {
				bxdf = b.bxdfs[i]
				break
			}
			count--
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

//func (b *bxDF) F(Wo, wi *Vector3f) Spectrum {
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

func NewFresnelDielectric(etaI, etaT float64) *FresnelDielectric {
	return &FresnelDielectric{
		etaI: etaI,
		etaT: etaT,
	}
}

type FresnelDielectric struct {
	etaI float64
	etaT float64
}

func (f *FresnelDielectric) Evaluate(cosI float64) Spectrum {
	return NewSpectrum(FrDielectric(cosI, f.etaI, f.etaT))
}

func NewSpecularTransmission(T Spectrum, etaA, etaB float64, mode TransportMode) *SpecularTransmission {
	return &SpecularTransmission{
		bxDF:    NewBxDF(BSDFTransmission | BSDFSpecular),
		T:       T,
		etaA:    etaA,
		etaB:    etaB,
		fresnel: NewFresnelDielectric(etaA, etaB),
		mode:    mode,
	}
}

type SpecularTransmission struct {
	*bxDF
	T          Spectrum
	etaA, etaB float64
	fresnel    *FresnelDielectric
	mode       TransportMode
}

func (t *SpecularTransmission) F(woW, wiW *Vector3f) Spectrum {
	return NewSpectrum(0)
}

func (t *SpecularTransmission) SampleF(wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType) {
	entering := CosTheta(wo) > 0
	var etaI, etaT float64
	if entering {
		etaI = t.etaA
		etaT = t.etaB
	} else {
		etaI = t.etaB
		etaT = t.etaA
	}

	// compute ray direction for specular transmission
	wi, refracts := Refract(wo, FaceForward(&Normal3f{0, 0, 1}, wo), etaI/etaT)
	if !refracts {
		return NewSpectrum(0), wi, 0, 0
	}

	ft := t.T.Mul(NewSpectrum(1).Sub(t.fresnel.Evaluate(CosTheta(wi))))
	// account for non-symmetry with transmission to different medium
	if t.mode == Radiance {
		ft = ft.MulScalar((etaI * etaI) / (etaT * etaT))
	}
	return ft.DivScalar(AbsCosTheta(wi)), wi, 1, 0
}

func (t *SpecularTransmission) Rho(wo *Vector3f, samples []*Point2f) Spectrum {
	return rho(t, wo, samples)
}

func (t *SpecularTransmission) RhoSamples(samples1, samples2 []*Point2f) Spectrum {
	return rhoSamples(t, samples1, samples2)
}

func (t *SpecularTransmission) Pdf(wo, wi *Vector3f) float64 {
	return 0
}

func NewFresnelSpecular(R, T Spectrum, etaA, etaB float64, mode TransportMode) *FresnelSpecular {
	return &FresnelSpecular{
		bxDF: NewBxDF(BSDFReflection | BSDFTransmission | BSDFSpecular),
		R:    R,
		T:    T,
		etaA: etaA,
		etaB: etaB,
		mode: mode,
	}
}

type FresnelSpecular struct {
	*bxDF

	R    Spectrum
	T    Spectrum
	etaA float64
	etaB float64
	mode TransportMode
}

func (f *FresnelSpecular) F(woW, wiW *Vector3f) Spectrum {
	return NewSpectrum(0.0)
}
func (f *FresnelSpecular) SampleF(wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType) {
	F := FrDielectric(CosTheta(wo), f.etaA, f.etaB)
	if sample.X < F {
		// compute specular reflection for FresnelSpecular
		// compute perfect specular reflection direction
		wi := &Vector3f{-wo.X, -wo.Y, wo.Z}
		return f.R.MulScalar(F).DivScalar(AbsCosTheta(wi)), wi, F, BSDFSpecular | BSDFReflection
	}

	// compute specular transmission for FresnelSpecular

	// figure out which eta is incident and which is transmitted
	var etaI, etaT float64
	entering := CosTheta(wo) > 0
	if entering {
		etaI = f.etaA
		etaT = f.etaB
	} else {
		etaI = f.etaB
		etaT = f.etaA
	}

	// compute ray direction for specular transmission
	wi, refracts := Refract(wo, FaceForward(&Normal3f{0, 0, 1}, wo), etaI/etaT)
	if !refracts {
		return NewSpectrum(0), wi, 0, 0
	}

	ft := f.T.MulScalar(1 - F)

	// account for non-symmetry with transmission to different medium
	if f.mode == Radiance {
		ft = ft.MulScalar((etaI * etaI) / (etaT / etaT))
	}
	return ft.DivScalar(AbsCosTheta(wi)), wi, 1 - F, BSDFSpecular | BSDFTransmission
}

func (f *FresnelSpecular) Rho(wo *Vector3f, samples []*Point2f) Spectrum {
	return rho(f, wo, samples)
}

func (f *FresnelSpecular) RhoSamples(samples1, samples2 []*Point2f) Spectrum {
	return rhoSamples(f, samples1, samples2)
}

func (f *FresnelSpecular) Pdf(wo, wi *Vector3f) float64 {
	return 0
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

func NewMicrofacetReflection(R Spectrum, d Distribution, f Fresnel) *MicrofacetReflection {
	return &MicrofacetReflection{
		bxDF:         NewBxDF(BSDFReflection | BSDFGlossy),
		R:            R,
		distribution: d,
		fresnel:      f,
	}
}

type MicrofacetReflection struct {
	*bxDF
	R            Spectrum
	distribution Distribution
	fresnel      Fresnel
}

func (m *MicrofacetReflection) GetType() BxDFType {
	return m.bxDF.Type
}

func (m *MicrofacetReflection) F(wo, wi *Vector3f) Spectrum {
	cosTheta0 := AbsCosTheta(wo)
	cosTheta1 := AbsCosTheta(wi)
	wh := wi.Add(wo)
	// handle degenerate cases for microfacet reflection
	if cosTheta1 == 0 || cosTheta0 == 0 {
		return NewSpectrum(0)
	}
	if wh.X == 0 && wh.Y == 0 && wh.Z == 0 {
		return NewSpectrum(0)
	}
	wh.Normalize()
	F := m.fresnel.Evaluate(wi.Dot(wh))
	return m.R.Mul(F).MulScalar(m.distribution.D(wh) * m.distribution.G(wo, wi) / (4 * cosTheta1 * cosTheta0))
}

func (m *MicrofacetReflection) SampleF(wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType) {
	// sample microfacet orientation wh and reflected direction wi
	if wo.Z == 0 {
		return NewSpectrum(0), nil, 0, 0
	}
	wh := m.distribution.SampleWH(wo, sample)
	wi = Reflect(wo, wh)
	if !SameHemisphere(wo, wi) {
		return NewSpectrum(0), nil, 0, 0
	}

	// compute PDF of wi for microfacet reflection
	pdf = m.distribution.Pdf(wo, wh) / (4 * wo.Dot(wh))
	return m.F(wo, wi), wi, pdf, 0
}

func (m *MicrofacetReflection) Rho(wo *Vector3f, samples []*Point2f) Spectrum {
	return rho(m, wo, samples)
}

func (m *MicrofacetReflection) RhoSamples(samples1, samples2 []*Point2f) Spectrum {
	return rhoSamples(m, samples1, samples2)
}

func (m *MicrofacetReflection) Pdf(wo, wi *Vector3f) float64 {
	if !SameHemisphere(wo, wi) {
		return 0
	}
	wh := wo.Add(wi).Normalized()
	return m.distribution.Pdf(wo, wh) / (4 * wo.Dot(wh))
}

func NewMicrofacetTransmission(T Spectrum, d Distribution, etaA, etaB float64, mode TransportMode) *MicrofacetTransmission {
	return &MicrofacetTransmission{
		bxDF:         NewBxDF(BSDFTransmission | BSDFGlossy),
		T:            T,
		distribution: d,
		etaA:         etaA,
		etaB:         etaB,
		fresnel:      NewFresnelDielectric(etaA, etaB),
	}
}

type MicrofacetTransmission struct {
	*bxDF
	T            Spectrum
	distribution Distribution
	etaA, etaB   float64
	fresnel      *FresnelDielectric
	mode         TransportMode
}

func (mt *MicrofacetTransmission) F(wo, wi *Vector3f) Spectrum {
	if !SameHemisphere(wo, wi) {
		return NewSpectrum(0) // transmission only
	}

	cosThetaO := CosTheta(wo)
	cosThetaI := CosTheta(wi)
	if cosThetaI == 0 || cosThetaO == 0 {
		return NewSpectrum(0)
	}

	// compute wh from Wo and wi for microfacet transmission
	var eta float64
	if CosTheta(wo) > 0 {
		eta = mt.etaA / mt.etaB
	} else {
		eta = mt.etaB / mt.etaA
	}
	wh := wo.Add(wi.MulScalar(eta))
	if wh.Z < 0 {
		wh = wh.MulScalar(-1)
	}

	F := mt.fresnel.Evaluate(wo.Dot(wh))

	sqrtDenom := wo.Dot(wh) * eta * wi.Dot(wh)
	factor := 1.0
	if mt.mode == Radiance {
		factor = 1 / eta
	}

	return NewSpectrum(1.0).Sub(F).Mul(mt.T).MulScalar(math.Abs(mt.distribution.D(wh) * mt.distribution.G(wo, wi) * eta * eta * wi.AbsDot(wh) * wo.AbsDot(wh) * factor * factor / (cosThetaI * cosThetaO * sqrtDenom * sqrtDenom)))
}

func (mt *MicrofacetTransmission) SampleF(wo *Vector3f, sample *Point2f) (s Spectrum, wi *Vector3f, pdf float64, sampledType BxDFType) {
	if wo.Z == 0 {
		return NewSpectrum(0), nil, 0, 0
	}
	wh := mt.distribution.SampleWH(wo, sample)
	var eta float64
	if CosTheta(wo) > 0 {
		eta = mt.etaA / mt.etaB
	} else {
		eta = mt.etaB / mt.etaA
	}

	wi, refracts := Refract(wo, wh, eta)
	if !refracts {
		return NewSpectrum(0), nil, 0, 0
	}

	return mt.F(wo, wi), wi, mt.Pdf(wo, wi), 0
}

func (mt *MicrofacetTransmission) Rho(wo *Vector3f, samples []*Point2f) Spectrum {
	return rho(mt, wo, samples)
}

func (mt *MicrofacetTransmission) RhoSamples(samples1, samples2 []*Point2f) Spectrum {
	return rhoSamples(mt, samples1, samples2)
}

func (mt *MicrofacetTransmission) Pdf(wo, wi *Vector3f) float64 {
	if SameHemisphere(wo, wi) {
		return 0
	}
	var eta float64
	if CosTheta(wo) > 0 {
		eta = mt.etaA / mt.etaB
	} else {
		eta = mt.etaB / mt.etaA
	}
	wh := wo.Add(wi.MulScalar(eta))

	sqrtDenom := wo.Dot(wh) + eta*wi.Dot(wh)
	dwhDwi := math.Abs((eta * eta * wi.Dot(wh)) / (sqrtDenom * sqrtDenom))
	return mt.distribution.Pdf(wo, wh) * dwhDwi
}
