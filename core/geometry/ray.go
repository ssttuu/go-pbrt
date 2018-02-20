package geometry

type Ray struct {
	origin    VectorXYZ
	direction VectorXYZ
	tMax      float64
	time      float64
	medium    Medium
}