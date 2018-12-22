package pbrt

type LightSampleStrategy int32

const (
	Uniform LightSampleStrategy = 1 << iota
	Power
	Spatial
)

func CreateLightSampleDistribution(strategy LightSampleStrategy, scene Scene) LightDistribution {
	switch strategy {
	case Uniform:
		return NewUniformLightDistribution(scene)
	case Power:
		return NewPowerLightDistribution(scene)
	}
	return nil
}

type LightDistribution interface {
	Lookup(p *Point3f) *Distribution1D
}

func NewUniformLightDistribution(scene Scene) *UniformLightDistribution {
	f := make([]float64, len(scene.Lights()))
	for i := 0; i < len(scene.Lights()); i++ {
		f[i] = 1.0
	}

	return &UniformLightDistribution{
		distribution: NewDistribution1D(f),
	}
}

type UniformLightDistribution struct {
	distribution *Distribution1D
}

func (d *UniformLightDistribution) Lookup(p *Point3f) *Distribution1D {
	return d.distribution
}

func NewPowerLightDistribution(scene Scene) *PowerLightDistribution {
	return &PowerLightDistribution{
		distribution: ComputeLightPowerDistribution(scene),
	}
}

type PowerLightDistribution struct {
	distribution *Distribution1D
}

func (d *PowerLightDistribution) Lookup(p *Point3f) *Distribution1D {
	return d.distribution
}

func ComputeLightPowerDistribution(scene Scene) *Distribution1D {
	if len(scene.Lights()) == 0 {
		return nil
	}

	lightPower := make([]float64, len(scene.Lights()))
	for _, light := range scene.Lights() {
		lightPower = append(lightPower, light.Power().Y())
	}
	return NewDistribution1D(lightPower)
}
