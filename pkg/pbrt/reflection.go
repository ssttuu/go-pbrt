package pbrt

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


type BSDF struct {
	eta float64
	ns, ng *Normal3f
	ss, ts *Vector3f
	nBxDFs int
	bxdfs [MaxBxDFs]BxDF
}

type BxDF struct {
	t BxDFType
}