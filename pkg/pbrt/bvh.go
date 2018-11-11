package pbrt

import (
	"math"
	"sync/atomic"
	"fmt"
)

type SplitMethod int

const (
	SplitSAH         SplitMethod = iota + 1
	SplitHLBVH
	SplitMiddle
	SplitEqualCounts
)

type BVHPrimitiveInfo struct {
	primitiveNumber int
	bounds          *Bounds3
	centroid        *Point3f
}

func NewBVHPrimitiveInfo(primitiveNumber int, bounds *Bounds3) *BVHPrimitiveInfo {
	return &BVHPrimitiveInfo{
		primitiveNumber: primitiveNumber,
		bounds:          bounds,
		centroid:        bounds.Min.MulScalar(0.5).Add(bounds.Max.MulScalar(0.5)),
	}
}

type BVHBuildNode struct {
	bounds          *Bounds3
	children        [2][]*BVHBuildNode
	splitAxis       int
	firstPrimOffset int64
	nPrimitives     int64
}

func (bn *BVHBuildNode) InitLeaf(first, n int64, b Bounds3) {
	bn.firstPrimOffset = first
	bn.nPrimitives = n
	bn.bounds = &b
	bn.children[0] = nil
	bn.children[1] = nil
}

func (bn *BVHBuildNode) InitInterior(axis int, c0, c1 []*BVHBuildNode) {
	bn.children[0] = c0
	bn.children[1] = c1
	bn.bounds = c0[0].bounds.Union(c1[0].bounds)
	bn.splitAxis = axis
	bn.nPrimitives = 0
}

type MortonPrimitive struct {
	primitiveIndex int
	mortonCode     uint64
}

type LBVHTreelet struct {
	startIndex  int64
	nPrimitives int64
	buildNodes  []*BVHBuildNode
}

type LinearBVHNode struct {
	bounds *Bounds3

	primitiveOffset, secondChildOffset int64
	nPrimitives                        uint64
	axis                               uint8
	pad                                [1]uint8
}

func LeftShift3(x uint64) uint64 {
	if x == 1<<10 {
		x--
	}

	// x = ---- --98 ---- ---- ---- ---- 7654 3210
	x = (x | (x << 16)) & 0x30000ff
	// x = ---- --98 ---- ---- 7654 ---- ---- 3210
	x = (x | (x << 8)) & 0x300f00f
	// x = ---- --98 ---- 76-- --54 ---- 32-- --10
	x = (x | (x << 4)) & 0x30c30c3
	// x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
	x = (x | (x << 2)) & 0x9249249

	return x
}

func EncodeMorton3(v *Vector3f) uint64 {
	return (LeftShift3(uint64(v.Z)) << 2) | (LeftShift3(uint64(v.Y)) << 1) | LeftShift3(uint64(v.X))
}

func RadixSortInPlace(v *[]*MortonPrimitive) {
	tempVector := make([]*MortonPrimitive, len(*v))

	// bitsPerPass must evenly divide nBits
	const (
		bitsPerPass = 6
		nBits       = 30
		nBuckets    = 1 << bitsPerPass
		bitMask     = (1 << bitsPerPass) - 1
	)

	var nPasses uint64 = nBits / bitsPerPass
	var pass uint64
	for pass = 0; pass < nPasses; pass++ {
		// perform one pass of radix sort, sorting by bitsPerPass bits
		var lowbit uint64 = pass * bitsPerPass

		// set in and out vector pointers for radix sort pass
		var in, out []*MortonPrimitive
		if pass&1 == 1 {
			in = tempVector
			out = *v
		} else {
			in = *v
			out = tempVector
		}

		// count number of zero bits in array for current radix sort bit
		var bucketCount [nBuckets]int
		for _, mp := range in {
			bucket := (mp.mortonCode >> lowbit) & bitMask
			bucketCount[bucket]++
		}

		// compute starting index in output array for each bucket
		var outIndex [nBuckets]int
		for i := 1; i < nBuckets; i++ {
			outIndex[i] = outIndex[i-1] + bucketCount[i-1]
		}

		// store sorted values in output array
		for _, mp := range in {
			bucket := (mp.mortonCode >> lowbit) & bitMask
			out[outIndex[bucket]] = mp
			outIndex[bucket]++
		}
	}

	if nPasses&1 == 1 {
		*v = tempVector
	}
}

func PartitionPrimitiveInfo(in []*BVHPrimitiveInfo, f func(pi *BVHPrimitiveInfo) bool) int64 {
	var first int64
	var last = int64(len(in) - 1)

OuterLoop:
	for first != last {
		for f(in[first]) {
			first++
			if first == last {
				break OuterLoop
			}
		}
		for !f(in[last]) {
			last--
			if first == last {
				break OuterLoop
			}
		}
		in[first], in[last] = in[last], in[first]
		first++
	}

	return first
}

func PartitionPrimitiveInfoAt(in []*BVHPrimitiveInfo, pivot int64, f func(a, b *BVHPrimitiveInfo) bool) int {
	first := 0
	last := len(in) - 1
	pivotValue := in[pivot]
	in[pivot], in[last] = in[last], in[pivot]
	//storeIndex := first

	for i := first; i < last; i++ {
		if f(in[i], pivotValue) {
			in[first], in[i] = in[i], in[first]
			first++
		}
	}
	in[last], in[first] = in[first], in[last]
	return first
}

func PartitionBuildNode(in []*BVHBuildNode, f func(pi *BVHBuildNode) bool) int64 {
	var first int64
	var last = int64(len(in) - 1)

OuterLoop:
	for first != last {
		for f(in[first]) {
			first++
			if first == last {
				break OuterLoop
			}
		}
		for !f(in[last]) {
			last--
			if first == last {
				break OuterLoop
			}
		}
		in[first], in[last] = in[last], in[first]
		first++
	}

	return first
}

func PartitionBuildNodeAt(in []*BVHBuildNode, pivot int64, f func(a, b *BVHBuildNode) bool) int {
	first := 0
	last := len(in) - 1
	pivotValue := in[pivot]
	in[pivot], in[last] = in[last], in[pivot]
	//storeIndex := first

	for i := first; i < last; i++ {
		if f(in[i], pivotValue) {
			in[first], in[i] = in[i], in[first]
			first++
		}
	}
	in[last], in[first] = in[first], in[last]
	return first
}

type BVHAccel struct {
	*Aggregate

	maxPrimsInNode int64
	splitMethod    SplitMethod
	primitives     []Primitiver
	nodes          []*LinearBVHNode
}

func NewBVHAccel(primitives []Primitiver, maxPrimsInNode int, splitMethod SplitMethod) *BVHAccel {
	bvh := &BVHAccel{
		maxPrimsInNode: int64(math.Min(255, float64(maxPrimsInNode))),
		splitMethod:    splitMethod,
		primitives:     primitives,
	}

	nPrimitives := len(primitives)

	if nPrimitives == 0 {
		return bvh
	}

	// build BVH for primitives

	// initialize primitiveInfo array for primitives
	primitiveInfo := make([]*BVHPrimitiveInfo, nPrimitives)

	for i := 0; i < nPrimitives; i++ {
		primitiveInfo[i] = NewBVHPrimitiveInfo(i, primitives[i].WorldBound())
	}

	// build BVH tree for primitives using primitiveInfo
	var totalNodes int64
	var orderedPrims []Primitiver

	var nodes []*BVHBuildNode
	if splitMethod == SplitHLBVH {
		nodes = bvh.HLBVHBuild(primitiveInfo, &totalNodes, &orderedPrims)
	} else {
		nodes = bvh.RecursiveBuild(primitiveInfo, &totalNodes, &orderedPrims)
	}
	fmt.Println(nodes)

	bvh.primitives = orderedPrims

	//offset := bvh.flattenBVHTree(nodes)

	return bvh
}

type BucketInfo struct {
	count  int
	bounds *Bounds3
}

func (b *BVHAccel) RecursiveBuild(primInfo []*BVHPrimitiveInfo, totalNodes *int64, orderedPrims *[]Primitiver) []*BVHBuildNode {
	nodes := make([]*BVHBuildNode, 1)
	*totalNodes++

	// compute bounds of all primitives in BVH node
	var bounds Bounds3
	for _, pi := range primInfo {
		bounds.Union(pi.bounds)
	}

	nPrimitives := int64(len(primInfo))

	if nPrimitives == 1 {
		// create leaf BVHBuildNode
		firstPrimOffset := int64(len(*orderedPrims))
		for _, pi := range primInfo {
			*orderedPrims = append(*orderedPrims, b.primitives[pi.primitiveNumber])
		}
		nodes[0].InitLeaf(firstPrimOffset, nPrimitives, bounds)
		return nodes
	}

	// compute bound of primitive centroids, choose split dimension dim
	var centroidBounds Bounds3
	for _, pi := range primInfo {
		centroidBounds.UnionPoint(pi.centroid)
	}
	dim := centroidBounds.MaximumExtent()

	// partition primitives into two sets and build children
	mid := nPrimitives / 2
	if centroidBounds.Max.Index(dim) == centroidBounds.Min.Index(dim) {
		// create leaf BVHBuildNode
		firstPrimOffset := int64(len(*orderedPrims))
		for _, pi := range primInfo {
			*orderedPrims = append(*orderedPrims, b.primitives[pi.primitiveNumber])
		}
		nodes[0].InitLeaf(firstPrimOffset, nPrimitives, bounds)
		return nodes
	}

	// partition primitives based on splitMethod
	switch b.splitMethod {
	case SplitMiddle:
		// partition primitives through node's midpoint
		pmid := (centroidBounds.Min.Index(dim) + centroidBounds.Min.Index(dim)) / 2

		mid := PartitionPrimitiveInfo(primInfo, func(pi *BVHPrimitiveInfo) bool {
			return pi.centroid.Index(dim) < pmid
		})
		// For lots of prims with large overlapping bounding boxes, this
		// may fail to partition; in that case don't break and fall
		// through to EqualCounts.
		if mid != 0 && mid != nPrimitives {
			break
		}
	case SplitEqualCounts:
		PartitionPrimitiveInfoAt(primInfo, mid, func(a, b *BVHPrimitiveInfo) bool {
			return a.centroid.Index(dim) < b.centroid.Index(dim)
		})
		break
	case SplitSAH:
	default:
		// partition primitives using approximate SAH
		if nPrimitives <= 2 {
			PartitionPrimitiveInfoAt(primInfo, mid, func(a, b *BVHPrimitiveInfo) bool {
				return a.centroid.Index(dim) < b.centroid.Index(dim)
			})
			break
		}

		// allocate bucketinfo for SAH partition buckets
		nBuckets := 12
		buckets := make([]*BucketInfo, nBuckets)

		// initialize bucketInfo for SAH parititon buckets
		for _, pi := range primInfo {
			b := nBuckets * int(centroidBounds.Offset(pi.centroid).Index(dim))
			if b == nBuckets {
				b = nBuckets - 1
			}
			buckets[b].count++
			buckets[b].bounds.Union(pi.bounds)
		}

		// compute costs for splitting after each bucket
		cost := make([]float64, nBuckets-1)
		for i := 0; i < nBuckets-1; i++ {
			var b0, b1 Bounds3
			var count0, count1 int
			for j := 0; j <= i; j++ {
				b0.Union(buckets[j].bounds)
				count0 += buckets[j].count
			}
			for j := i + 1; j < nBuckets; j++ {
				b1.Union(buckets[j].bounds)
				count1 += buckets[j].count
			}
			cost[i] = 1.0 + (float64(count0)*b0.SurfaceArea()+float64(count1)*b1.SurfaceArea())/bounds.SurfaceArea()
		}

		// find bucket to split at that minimizes SAH metric
		minCost := cost[0]
		minCostSplitBucket := 0
		for i := 1; i < nBuckets; i++ {
			if cost[i] < minCost {
				minCost = cost[i]
				minCostSplitBucket = i
			}
		}

		// either create leaf or split primitives at selected SAH bucket
		leafCost := float64(nPrimitives)
		if nPrimitives > b.maxPrimsInNode || minCost < leafCost {
			mid = PartitionPrimitiveInfo(primInfo, func(pi *BVHPrimitiveInfo) bool {
				b := nBuckets * int(centroidBounds.Offset(pi.centroid).Index(dim))
				if b == nBuckets {
					b = nBuckets - 1
				}
				return b <= minCostSplitBucket
			})

		} else {
			// create leaf BVHBuildNode
			firstPrimOffset := int64(len(*orderedPrims))
			for _, pi := range primInfo {
				*orderedPrims = append(*orderedPrims, b.primitives[pi.primitiveNumber])
			}
			nodes[0].InitLeaf(firstPrimOffset, nPrimitives, bounds)
			return nodes
		}
		break
	}

	nodes[0].InitInterior(dim,
		b.RecursiveBuild(primInfo[:mid], totalNodes, orderedPrims),
		b.RecursiveBuild(primInfo[mid:], totalNodes, orderedPrims),
	)

	return nodes

}

func (b *BVHAccel) HLBVHBuild(primInfo []*BVHPrimitiveInfo, totalNodes *int64, orderedPrims *[]Primitiver) []*BVHBuildNode {
	var bounds *Bounds3
	for _, pi := range primInfo {
		bounds.UnionPoint(pi.centroid)
	}

	// compute morton indices of primitives
	mortonPrims := make([]*MortonPrimitive, len(primInfo))
	// TODO: parallel
	for i, pi := range primInfo {
		var mortonBits uint = 10
		mortonScale := 1 << mortonBits
		mortonPrims[i].primitiveIndex = pi.primitiveNumber
		centroidOffset := bounds.Offset(pi.centroid)
		mortonPrims[i].mortonCode = EncodeMorton3(centroidOffset.MulScalar(float64(mortonScale)))
	}

	RadixSortInPlace(&mortonPrims)

	// create LBVH treelets at bottom of BVH

	// find intervals of primitives for each treelet
	var treeletsToBuild []*LBVHTreelet
	start := int64(0)
	for end := int64(0); end <= int64(len(mortonPrims)); end++ {

		var mask uint64 = 0x3ffc0000

		if end == int64(len(mortonPrims)) || ((mortonPrims[start].mortonCode & mask) != (mortonPrims[end].mortonCode & mask)) {
			nPrimitives := end - start
			maxBVHNodes := 2 * nPrimitives
			nodes := make([]*BVHBuildNode, maxBVHNodes)
			treeletsToBuild = append(treeletsToBuild, &LBVHTreelet{start, nPrimitives, nodes})
			start = end
		}
	}

	var orderedPrimsOffset int64

	// Create LBVHs for treelets in parallel
	createLBVHTreelet := func(i int) {
		// generate ith LBVH treelet
		var nodesCreated int64
		var firstBitIndex int64 = 29 - 12
		tr := treeletsToBuild[i]
		tr.buildNodes = b.emitLBVH(
			tr.buildNodes,
			primInfo,
			mortonPrims[tr.startIndex:tr.startIndex+tr.nPrimitives],
			tr.nPrimitives,
			&nodesCreated,
			orderedPrims,
			&orderedPrimsOffset,
			firstBitIndex,
		)
		atomic.AddInt64(totalNodes, nodesCreated)
	}
	// TODO: run in parallel
	for i := 0; i < len(treeletsToBuild); i++ {
		createLBVHTreelet(i)
	}

	// create and return SAH BVH from LBVH treelets
	finishedTreelets := make([]*BVHBuildNode, len(treeletsToBuild))
	for _, treelet := range treeletsToBuild {
		finishedTreelets = append(finishedTreelets, treelet.buildNodes...)
	}

	return b.buildUpperSAH(finishedTreelets, totalNodes)
}

func (b *BVHAccel) emitLBVH(
	buildNodes []*BVHBuildNode, primInfo []*BVHPrimitiveInfo, mortonPrims []*MortonPrimitive,
	nPrimitives int64, totalNodes *int64, orderedPrims *[]Primitiver, orderedPrimsOffset *int64, bitIndex int64) []*BVHBuildNode {

	if bitIndex == -1 || nPrimitives < b.maxPrimsInNode {
		// create and return leaf node of LBVH treelet
		*totalNodes++
		node := buildNodes[1]
		var bounds Bounds3
		firstPrimOffset := atomic.AddInt64(orderedPrimsOffset, nPrimitives)
		for i := int64(0); i < nPrimitives; i++ {
			primitiveIndex := mortonPrims[i].primitiveIndex
			(*orderedPrims)[firstPrimOffset+i] = b.primitives[primitiveIndex]
			bounds.Union(primInfo[primitiveIndex].bounds)
		}
		node.InitLeaf(firstPrimOffset, nPrimitives, bounds)
		return buildNodes[1:] // TODO: double check this is right
	}

	var mask uint64 = 1 << uint64(bitIndex)
	// Advance to next subtree level if there's no LBVH split for this bit
	if (mortonPrims[0].mortonCode & mask) == (mortonPrims[nPrimitives-1].mortonCode & mask) {
		return b.emitLBVH(buildNodes, primInfo, mortonPrims, nPrimitives, totalNodes, orderedPrims, orderedPrimsOffset, bitIndex-1)
	}

	// find LBVH split point for this dimension
	var searchStart, searchEnd int64 = 0, nPrimitives-1
	for searchStart+1 != searchEnd {
		mid := nPrimitives / 2
		if (mortonPrims[searchStart].mortonCode & mask) == (mortonPrims[mid].mortonCode & mask) {
			searchStart = mid
		} else {
			searchEnd = mid
		}
	}
	splitOffset := searchEnd

	// create and return interior LBVH node
	*totalNodes++
	lbvh := [2][]*BVHBuildNode{
		b.emitLBVH(buildNodes[1:], primInfo, mortonPrims[:splitOffset], splitOffset, totalNodes, orderedPrims, orderedPrimsOffset, bitIndex-1),
		b.emitLBVH(buildNodes[1:], primInfo, mortonPrims[splitOffset:], nPrimitives-splitOffset, totalNodes, orderedPrims, orderedPrimsOffset, bitIndex-1),
	}
	axis := bitIndex % 3
	buildNodes[0].InitInterior(int(axis), lbvh[0], lbvh[1])

	return buildNodes
}

func (b *BVHAccel) buildUpperSAH(treeletRoots []*BVHBuildNode, totalNodes *int64) []*BVHBuildNode {
	nNodes := len(treeletRoots)
	if nNodes == 1 {
		return treeletRoots
	}
	*totalNodes++

	// compute bounds of all nodes under this HLBVH node
	var bounds Bounds3
	for _, tr := range treeletRoots {
		bounds.Union(tr.bounds)
	}

	// compute bound of HLBVH node centroids, choose split dimension
	var centroidBounds Bounds3
	for _, tr := range treeletRoots {
		centroidBounds.UnionPoint(tr.bounds.Min.Add(tr.bounds.Max).MulScalar(0.5))
	}
	dim := centroidBounds.MaximumExtent()

	// allocate BucketInfo for SAH partition buckets
	nBuckets := 12
	buckets := make([]*BucketInfo, nBuckets)

	// initialize BucketInfo for HLBVH SAH partition buckets
	for _, tr := range treeletRoots {
		centroid := (tr.bounds.Min.Index(dim) + tr.bounds.Max.Index(dim)) * 0.5
		b := nBuckets * int((centroid-centroidBounds.Min.Index(dim))/centroidBounds.Max.Index(dim)-centroidBounds.Min.Index(dim))
		if b == nBuckets {
			b = nBuckets - 1
		}

		buckets[b].count++
		buckets[b].bounds.Union(tr.bounds)
	}

	// compute costs for splitting after each bucket
	cost := make([]float64, nBuckets-1)
	for i := 0; i < nBuckets; i++ {
		var b0, b1 Bounds3
		var count0, count1 int
		for j := 0; j <= i; j++ {
			b0.Union(buckets[j].bounds)
			count0 += buckets[j].count
		}
		for j := i + 1; j < nBuckets; j++ {
			b1.Union(buckets[j].bounds)
			count1 += buckets[j].count
		}
		cost[i] = 0.125 + (float64(count0)*b0.SurfaceArea()+float64(count1)*b1.SurfaceArea())/bounds.SurfaceArea()
	}

	// find bucket to split at that minimizes SAH metric
	minCost := cost[0]
	minCostSplitBucket := 0
	for i := 1; i < nBuckets-1; i++ {
		if cost[i] < minCost {
			minCost = cost[i]
			minCostSplitBucket = i
		}
	}

	// split nodes and create interior HLBVH SAH node
	mid := PartitionBuildNode(treeletRoots, func(node *BVHBuildNode) bool {
		centroid := (node.bounds.Min.Index(dim) + node.bounds.Max.Index(dim)) * 0.5
		b := nBuckets * int((centroid-centroidBounds.Min.Index(dim))/(centroidBounds.Max.Index(dim)-centroidBounds.Min.Index(dim)))
		if b == nBuckets {
			b = nBuckets - 1
		}

		return b <= minCostSplitBucket
	})

	treeletRoots[0].InitInterior(
		dim,
		b.buildUpperSAH(treeletRoots[:mid], totalNodes),
		b.buildUpperSAH(treeletRoots[mid:], totalNodes),
	)
	return treeletRoots
}

//func (bn *BVHBuildNode) emitLBVH(buildNodes )

//func (b *BVHAccel) WorldBound() Bounds3 {
//
//}
//func (b *BVHAccel) Intersect(r *Ray) (bool, *SurfaceInteraction) {
//
//}
//func (b *BVHAccel) IntersectP(r *Ray) bool {
//
//}
//
//func (b *BVHAccel) GetAreaLight() *AreaLight {
//	log.Panic("Aggregate.GetAreaLight called; should have gone to Geometric Primitiver")
//	return nil
//}
//
//func (b *BVHAccel) GetMaterial() *Material {
//	log.Panic("Aggregate.GetMaterial called; should have gone to GeometricPrimitive")
//	return nil
//}
//
//func (b *BVHAccel) ComputeScatteringFunctions(si *SurfaceInteraction, arena *MemoryArena, mode TransportMode, allowMultipleLobes bool) {
//	log.Panic("Aggregate.ComputeScatteringFunctions called; should have gone to GeometricPrimitive")
//}
