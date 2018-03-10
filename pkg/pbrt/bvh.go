package pbrt
//
//import (
//	"math"
//	"log"
//	"sync/atomic"
//)
//
//type SplitMethod int
//
//const (
//	SplitSAH         SplitMethod = iota + 1
//	SplitHLBVH
//	SplitMiddle
//	SplitEqualCounts
//)
//
//type PrimitiveInfo struct {
//	primitiveNumber int
//	bounds          *Bounds3
//	centroid        *Point3f
//}
//
//func NewPrimitiveInfo(primitiveNumber int, bounds *Bounds3) *PrimitiveInfo {
//	return &PrimitiveInfo{
//		primitiveNumber: primitiveNumber,
//		bounds:          bounds,
//		centroid:        bounds.min.MulScalar(0.5).Add(bounds.max.MulScalar(0.5))
//	}
//}
//
//type MortonPrimitive struct {
//	primitiveIndex int
//	mortonCode     uint32
//}
//
//type BVHBuildNode struct {
//	bounds                                  *Bounds3
//	children                                [2]*BVHBuildNode
//	splitAxis, firstPrimOffset, nPrimitives int
//}
//
//func (bn *BVHBuildNode) InitLeaf(first, n int, b *Bounds3) {
//	bn.firstPrimOffset = first
//	bn.nPrimitives = n
//	bn.bounds = b
//	bn.children[0] = nil
//	bn.children[1] = nil
//}
//
//func (bn *BVHBuildNode) InitInterior(axis int, c0, c1 *BVHBuildNode) {
//	bn.children[0] = c0
//	bn.children[1] = c1
//	bn.bounds = c0.bounds.Union(c1.bounds)
//	bn.splitAxis = axis
//	bn.nPrimitives = 0
//}
//
//type LBVHTreelet struct {
//	startIndex, nPrimitives int
//	buildNodes              []*BVHBuildNode
//}
//
//type BVHAccel struct {
//	*Aggregate
//
//	maxPrimsInNode int
//	splitMethod    SplitMethod
//	primitives     []Primitive
//	nodes          []*LinearBVHNode
//}
//
//func NewBVHAccel(primitives []Primitive, maxPrimsInNode int, splitMethod SplitMethod) *BVHAccel {
//	bvh := &BVHAccel{
//		maxPrimsInNode: int(math.Min(255, float64(maxPrimsInNode))),
//		splitMethod:    splitMethod,
//		primitives:     primitives,
//	}
//
//	nPrimitives := len(primitives)
//
//	if nPrimitives == 0 {
//		return bvh
//	}
//
//	// build BVH for primitives
//
//	// initialize primitiveInfo array for primitives
//	primitiveInfo := make([]*PrimitiveInfo, nPrimitives)
//
//	for i := 0; i < nPrimitives; i++ {
//		primitiveInfo[i] = NewPrimitiveInfo(i, primitives[i].WorldBound())
//	}
//
//	// build BVH tree for primitives using primitiveInfo
//	var totalNodes int
//	var orderedPrims []Primitive
//
//	root := &BVHBuildNode{}
//	if splitMethod == SplitHLBVH {
//		root, totalNodes, orderedPrims = bvh.HLBVHBuild(primitiveInfo)
//	} else {
//		root, totalNodes, orderedPrims = bvh.RecursiveBuild(primitiveInfo, 0, nPrimitives)
//	}
//
//	bvh.primitives = orderedPrims
//
//	offset := bvh.flattenBVHTree(root)
//
//	return bvh
//}
//
//func (b *BVHAccel) HLBVHBuild(info []PrimitiveInfo) (*BVHBuildNode, int, []Primitive) {
//	var bounds *Bounds3
//	for _, pi := range info {
//		bounds.UnionPoint(pi.centroid)
//	}
//
//	// compute morton indices of primitives
//	mortonPrims := make([]*MortonPrimitive, len(info))
//	// TODO: parallel
//	for i, pi := range info {
//		var mortonBits uint = 10
//		mortonScale := 1 << mortonBits
//		mortonPrims[i].primitiveIndex = pi.primitiveNumber
//		centroidOffset := bounds.Offset(pi.centroid)
//		mortonPrims[i].mortonCode = EncodeMorton3(centroidOffset.MulScalar(float64(mortonScale)))
//	}
//
//	RadixSort(mortonPrims)
//
//	// create LBVH treelets at bottom of BVH
//
//	// find intervals of primitives for each treelet
//	var treeletsToBuild []*LBVHTreelet
//	start := 0
//	for end := 0; end <= len(mortonPrims); end++ {
//
//		var mask uint32 = 0x3ffc0000
//
//		if end == len(mortonPrims) || ((mortonPrims[start].mortonCode & mask) != (mortonPrims[end].mortonCode & mask)) {
//			nPrimitives := end - start
//			maxBVHNodes := 2 * nPrimitives
//			nodes := make([]*BVHBuildNode, maxBVHNodes)
//			treeletsToBuild = append(treeletsToBuild, &LBVHTreelet{start, nPrimitives, nodes})
//			start = end
//		}
//	}
//
//	// Create LBVHs for treelets in parallel
//	var totalNodes, orderedPrimsOffset int64
//	run := func(i int) {
//		// generate ith LBVH treelet
//		var nodesCreated int64 = 0
//		firstBitIndex := 29 - 12
//		tr := treeletsToBuild[i]
//		tr.buildNodes = b.emitLBVH()
//		atomic.AddInt64(&totalNodes, nodesCreated)
//	}
//
//
//	// create and return SAH BVH from LBVH treelets
//	finishedTreelets := make([]*BVHBuildNode, len(treeletsToBuild))
//	for i, treelet := range treeletsToBuild {
//		finishedTreelets = append(finishedTreelets, treelet.buildNodes...)
//	}
//
//	return bn.BuildUpperSAH(finishedTreelets, 0, len(finishedTreelets), totalNodes), totalNodes, primit
//
//
//}
//
//func (b *BVHBuildNode) emitLBVH(buildNodes )
//
//func (a *BVHAccel) WorldBound() Bounds3 {
//
//}
//func (a *BVHAccel) Intersect(r *Ray) (bool, *SurfaceInteraction) {
//
//}
//func (a *BVHAccel) IntersectP(r *Ray) bool {
//
//}
//
//func (a *BVHAccel) GetAreaLight() *AreaLight {
//	log.Panic("Aggregate.GetAreaLight called; should have gone to Geometric Primitive")
//	return nil
//}
//
//func (a *BVHAccel) GetMaterial() *Material {
//	log.Panic("Aggregate.GetMaterial called; should have gone to GeometricPrimitive")
//	return nil
//}
//
//func (a *BVHAccel) ComputeScatteringFunctions(si *SurfaceInteraction, arena *MemoryArena, mode TransportMode, allowMultipleLobes bool) {
//	log.Panic("Aggregate.ComputeScatteringFunctions called; should have gone to GeometricPrimitive")
//}
