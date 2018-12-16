package accelerator

import (
	"context"
	"log"
	"math"
	"sync/atomic"

	"golang.org/x/sync/errgroup"

	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

type SplitMethod int

const (
	SplitSAH SplitMethod = iota + 1
	SplitHLBVH
	SplitMiddle
	SplitEqualCounts
)

type BVHPrimitiveInfo struct {
	primitiveNumber int
	bounds          *pbrt.Bounds3
	centroid        *pbrt.Point3f
}

func NewBVHPrimitiveInfo(primitiveNumber int, bounds *pbrt.Bounds3) *BVHPrimitiveInfo {
	return &BVHPrimitiveInfo{
		primitiveNumber: primitiveNumber,
		bounds:          bounds,
		centroid:        bounds.Min.MulScalar(0.5).Add(bounds.Max.MulScalar(0.5)),
	}
}

type BVHBuildNode struct {
	bounds          pbrt.Bounds3
	children        [2]*BVHBuildNode
	splitAxis       uint8
	firstPrimOffset int64
	nPrimitives     int64
}

func (bn *BVHBuildNode) InitLeaf(first, n int64, b pbrt.Bounds3) {
	bn.firstPrimOffset = first
	bn.nPrimitives = n
	bn.bounds = b

	bn.children[0] = nil
	bn.children[1] = nil
}

func (bn *BVHBuildNode) InitInterior(axis uint8, c0, c1 *BVHBuildNode) {
	bn.children[0] = c0
	bn.children[1] = c1

	bn.bounds = pbrt.Bounds3{
		Min: c0.bounds.Min,
		Max: c0.bounds.Max,
	}
	bn.bounds.Union(&c1.bounds)

	bn.splitAxis = axis
	bn.nPrimitives = 0
}

type MortonPrimitive struct {
	primitiveIndex int
	mortonCode     uint64
}

type LBVHTreelet struct {
	startIndex   int64
	nPrimitives  int64
	buildNodes   []BVHBuildNode
	finishedNode *BVHBuildNode
}

type LinearBVHNode struct {
	bounds pbrt.Bounds3

	primitiveOffset, secondChildOffset uint64
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

func EncodeMorton3(v *pbrt.Vector3f) uint64 {
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

	if (nPasses & 1) == 1 {
		*v = tempVector
	}
}

func PartitionPrimitiveInfoAt(in []*BVHPrimitiveInfo, start, end, pivot int64, f func(a, b *BVHPrimitiveInfo) bool) int64 {
	pivotValue := in[pivot]
	in[pivot], in[end] = in[end], in[pivot]

	for i := start; i < end; i++ {
		if f(in[i], pivotValue) {
			in[start], in[i] = in[i], in[start]
			start++
		}
	}
	in[end], in[start] = in[start], in[end]
	return start
}

func PartitionBuildNode(in []BVHBuildNode, f func(pi *BVHBuildNode) bool) int64 {
	var first int64
	var last = int64(len(in) - 1)

OuterLoop:
	for first != last {
		for f(&in[first]) {
			first++
			if first == last {
				break OuterLoop
			}
		}
		for !f(&in[last]) {
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

func PartitionBuildNodeAt(in []*BVHBuildNode, start, end, pivot int64, f func(a, b *BVHBuildNode) bool) int64 {
	pivotValue := in[pivot]
	in[pivot], in[end] = in[end], in[pivot]

	for i := start; i < end; i++ {
		if f(in[i], pivotValue) {
			in[start], in[i] = in[i], in[start]
			start++
		}
	}
	in[end], in[start] = in[start], in[end]
	return start
}

type BVH struct {
	splitMethod    SplitMethod
	maxPrimsInNode uint8
	primitives     []pbrt.Primitive
	nodes          []LinearBVHNode
}

func NewBVH(primitives []pbrt.Primitive, maxPrimsInNode int, splitMethod SplitMethod) *BVH {
	bvh := &BVH{
		maxPrimsInNode: uint8(math.Min(255, float64(maxPrimsInNode))),
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
	var orderedPrims []pbrt.Primitive

	var root *BVHBuildNode
	if splitMethod == SplitHLBVH {
		root = bvh.HLBVHBuild(primitiveInfo, &totalNodes, &orderedPrims)
	} else {
		root = bvh.RecursiveBuild(primitiveInfo, 0, int64(len(primitives)), &totalNodes, &orderedPrims)
	}

	bvh.primitives = orderedPrims

	// compute representation of depth-first traversal of BVH tree

	bvh.nodes = make([]LinearBVHNode, totalNodes)
	var offset uint64
	bvh.flattenBVHTree(root, &offset)

	return bvh
}

type BucketInfo struct {
	count  int
	bounds pbrt.Bounds3
}

func (b *BVH) RecursiveBuild(primitiveInfo []*BVHPrimitiveInfo, start, end int64, totalNodes *int64, orderedPrims *[]pbrt.Primitive) *BVHBuildNode {
	node := &BVHBuildNode{}
	*totalNodes++

	// compute bounds of all primitives in BVH node
	var bounds pbrt.Bounds3
	for i := start; i < end; i++ {
		bounds.Union(primitiveInfo[i].bounds)
	}

	nPrimitives := end - start
	if nPrimitives == 1 {
		// create leaf BVHBuildNode
		firstPrimOffset := int64(len(*orderedPrims))
		for i := start; i < end; i++ {
			primitiveNumber := primitiveInfo[i].primitiveNumber
			*orderedPrims = append(*orderedPrims, b.primitives[primitiveNumber])
		}
		node.InitLeaf(firstPrimOffset, nPrimitives, bounds)
		return node
	}

	// compute bound of Primitive centroids, choose split dimension dim
	var centroidBounds pbrt.Bounds3
	for i := start; i < end; i++ {
		centroidBounds.UnionPoint(primitiveInfo[i].centroid)
	}
	dim := centroidBounds.MaximumExtent()

	// partition primitives into two sets and build children
	mid := (start + end) / 2
	if centroidBounds.Max.Index(dim) == centroidBounds.Min.Index(dim) {
		// create leaf BVHBuildNode
		firstPrimOffset := int64(len(*orderedPrims))
		for i := start; i < end; i++ {
			*orderedPrims = append(*orderedPrims, b.primitives[primitiveInfo[i].primitiveNumber])
		}
		node.InitLeaf(firstPrimOffset, nPrimitives, bounds)
		return node
	}

	// partition primitives based on splitMethod
	switch b.splitMethod {
	case SplitMiddle:
		// partition primitives through node's midpoint
		pmid := (centroidBounds.Min.Index(dim) + centroidBounds.Min.Index(dim)) / 2

		mid := PartitionPrimitiveInfoAt(primitiveInfo, start, end-1, end-1, func(a, b *BVHPrimitiveInfo) bool {
			return a.centroid.Index(dim) < pmid
		})
		// For lots of prims with large overlapping bounding boxes, this
		// may fail to partition; in that case don't break and fall
		// through to EqualCounts.
		if mid != 0 && mid != nPrimitives {
			break
		}

		fallthrough
	case SplitEqualCounts:
		PartitionPrimitiveInfoAt(primitiveInfo, start, end-1, mid, func(a, b *BVHPrimitiveInfo) bool {
			return a.centroid.Index(dim) < b.centroid.Index(dim)
		})
		break
	case SplitSAH:
		// partition primitives using approximate SAH
		if nPrimitives <= 2 {
			PartitionPrimitiveInfoAt(primitiveInfo, start, end-1, mid, func(a, b *BVHPrimitiveInfo) bool {
				return a.centroid.Index(dim) < b.centroid.Index(dim)
			})
			break
		}

		// allocate bucketinfo for SAH partition buckets
		nBuckets := 12
		buckets := make([]BucketInfo, nBuckets)
		// initialize bucketInfo for SAH partition buckets
		for i := start; i < end; i++ {
			b := nBuckets * int(centroidBounds.Offset(primitiveInfo[i].centroid).Index(dim))
			if b == nBuckets {
				b = nBuckets - 1
			}
			buckets[b].count++
			buckets[b].bounds.Union(primitiveInfo[i].bounds)
		}

		// compute costs for splitting after each bucket
		cost := make([]float64, nBuckets-1)
		for i := 0; i < nBuckets-1; i++ {
			var b0, b1 pbrt.Bounds3
			var count0, count1 int
			for j := 0; j <= i; j++ {
				b0.Union(&buckets[j].bounds)
				count0 += buckets[j].count
			}
			for j := i + 1; j < nBuckets; j++ {
				b1.Union(&buckets[j].bounds)
				count1 += buckets[j].count
			}
			cost[i] = 1.0 + (float64(count0)*b0.SurfaceArea()+float64(count1)*b1.SurfaceArea())/bounds.SurfaceArea()
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

		// either create leaf or split primitives at selected SAH bucket
		leafCost := float64(nPrimitives)
		if nPrimitives > int64(b.maxPrimsInNode) || minCost < leafCost {
			mid = PartitionPrimitiveInfoAt(primitiveInfo, start, end-1, end-1, func(a, b *BVHPrimitiveInfo) bool {
				c := nBuckets * int(centroidBounds.Offset(a.centroid).Index(dim))
				if c == nBuckets {
					c = nBuckets - 1
				}
				return c <= minCostSplitBucket
			})
		} else {
			// create leaf BVHBuildNode
			firstPrimOffset := int64(len(*orderedPrims))
			for i := start; i < end; i++ {
				*orderedPrims = append(*orderedPrims, b.primitives[primitiveInfo[i].primitiveNumber])
			}
			node.InitLeaf(firstPrimOffset, nPrimitives, bounds)
			return node
		}
	}

	node.InitInterior(uint8(dim),
		b.RecursiveBuild(primitiveInfo, start, mid, totalNodes, orderedPrims),
		b.RecursiveBuild(primitiveInfo, mid, end, totalNodes, orderedPrims),
	)

	return node

}

func (b *BVH) HLBVHBuild(primInfo []*BVHPrimitiveInfo, totalNodes *int64, orderedPrims *[]pbrt.Primitive) *BVHBuildNode {
	var bounds pbrt.Bounds3
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
	for start, end := 0, 0; end <= len(mortonPrims); end++ {
		var mask uint64 = 0x3ffc0000

		if end == len(mortonPrims) || ((mortonPrims[start].mortonCode & mask) != (mortonPrims[end].mortonCode & mask)) {
			nPrimitives := end - start
			maxBVHNodes := 2 * nPrimitives
			treeletsToBuild = append(treeletsToBuild, &LBVHTreelet{
				startIndex:  int64(start),
				nPrimitives: int64(nPrimitives),
				buildNodes:  make([]BVHBuildNode, maxBVHNodes),
			})
			start = end
		}
	}

	var orderedPrimsOffset int64
	*orderedPrims = make([]pbrt.Primitive, len(b.primitives))

	// Create LBVHs for treelets in parallel
	g, _ := errgroup.WithContext(context.Background())
	for i := 0; i < len(treeletsToBuild); i++ {
		g.Go(func(i int) func() error {
			return func() error {
				// generate ith LBVH treelet
				var nodesCreated int64
				var firstBitIndex int64 = 29 - 12
				tr := treeletsToBuild[i]
				tr.finishedNode = b.emitLBVH(
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
				return nil
			}
		}(i))
	}

	err := g.Wait()
	if err != nil {
		// TODO
	}

	// create and return SAH BVH from LBVH treelets
	finishedTreelets := make([]BVHBuildNode, len(treeletsToBuild))
	for _, treelet := range treeletsToBuild {
		finishedTreelets = append(finishedTreelets, treelet.buildNodes...)
	}

	return b.buildUpperSAH(finishedTreelets, 0, int64(len(finishedTreelets)), totalNodes)
}

func (b *BVH) emitLBVH(
	buildNodes []BVHBuildNode,
	primInfo []*BVHPrimitiveInfo,
	mortonPrims []*MortonPrimitive,
	nPrimitives int64,
	totalNodes *int64,
	orderedPrims *[]pbrt.Primitive,
	orderedPrimsOffset *int64,
	bitIndex int64) *BVHBuildNode {

	if bitIndex == -1 || nPrimitives < int64(b.maxPrimsInNode) {
		// create and return leaf node of LBVH treelet
		*totalNodes++
		node := buildNodes[0]
		bounds := *primInfo[0].bounds
		firstPrimOffset := atomic.AddInt64(orderedPrimsOffset, nPrimitives)
		for i := int64(0); i < nPrimitives; i++ {
			primitiveIndex := mortonPrims[i].primitiveIndex
			(*orderedPrims)[firstPrimOffset+i] = b.primitives[primitiveIndex]
			bounds.Union(primInfo[primitiveIndex].bounds)
		}
		node.InitLeaf(firstPrimOffset, nPrimitives, bounds)
		return &node
	}

	var mask uint64 = 1 << uint64(bitIndex)
	// Advance to next subtree level if there's no LBVH split for this bit
	if (mortonPrims[0].mortonCode & mask) == (mortonPrims[nPrimitives-1].mortonCode & mask) {
		return b.emitLBVH(buildNodes, primInfo, mortonPrims, nPrimitives, totalNodes, orderedPrims, orderedPrimsOffset, bitIndex-1)
	}

	// find LBVH split Point for this dimension
	var searchStart, searchEnd int64 = 0, nPrimitives - 1
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
	node := buildNodes[0]
	*totalNodes++
	lbvh := [2]*BVHBuildNode{
		b.emitLBVH(buildNodes[1:], primInfo, mortonPrims[:splitOffset], splitOffset, totalNodes, orderedPrims, orderedPrimsOffset, bitIndex-1),
		b.emitLBVH(buildNodes[1:], primInfo, mortonPrims[splitOffset:], nPrimitives-splitOffset, totalNodes, orderedPrims, orderedPrimsOffset, bitIndex-1),
	}
	axis := uint8(bitIndex % 3)
	node.InitInterior(axis, lbvh[0], lbvh[1])

	return &node
}

func (b *BVH) buildUpperSAH(treeletRoots []BVHBuildNode, start, end int64, totalNodes *int64) *BVHBuildNode {
	nNodes := end - start
	if nNodes == 1 {
		return &treeletRoots[start]
	}
	*totalNodes++
	node := &BVHBuildNode{}

	// compute bounds of all nodes under this HLBVH node
	bounds := treeletRoots[start].bounds
	for i := start; i < end; i++ {
		bounds.Union(&treeletRoots[i].bounds)
	}

	// compute bound of HLBVH node centroids, choose split dimension
	var centroidBounds pbrt.Bounds3
	for i := start; i < end; i++ {
		centroidBounds.UnionPoint(treeletRoots[i].bounds.Min.Add(treeletRoots[i].bounds.Max).MulScalar(0.5))
	}
	dim := centroidBounds.MaximumExtent()

	// allocate BucketInfo for SAH partition buckets
	nBuckets := 12
	buckets := make([]BucketInfo, nBuckets)

	// initialize BucketInfo for HLBVH SAH partition buckets
	for i := start; i < end; i++ {
		centroid := (treeletRoots[i].bounds.Min.Index(dim) + treeletRoots[i].bounds.Max.Index(dim)) * 0.5
		b := nBuckets * int((centroid-centroidBounds.Min.Index(dim))/centroidBounds.Max.Index(dim)-centroidBounds.Min.Index(dim))
		if b == nBuckets {
			b = nBuckets - 1
		}

		buckets[b].count++
		buckets[b].bounds.Union(&treeletRoots[i].bounds)
	}

	// compute costs for splitting after each bucket
	cost := make([]float64, nBuckets-1)
	for i := 0; i < nBuckets; i++ {
		var b0, b1 pbrt.Bounds3
		var count0, count1 int
		for j := 0; j <= i; j++ {
			b0.Union(&buckets[j].bounds)
			count0 += buckets[j].count
		}
		for j := i + 1; j < nBuckets; j++ {
			b1.Union(&buckets[j].bounds)
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
	mid := PartitionBuildNode(treeletRoots[start:end], func(node *BVHBuildNode) bool {
		centroid := (node.bounds.Min.Index(dim) + node.bounds.Max.Index(dim)) * 0.5
		b := nBuckets * int((centroid-centroidBounds.Min.Index(dim))/(centroidBounds.Max.Index(dim)-centroidBounds.Min.Index(dim)))
		if b == nBuckets {
			b = nBuckets - 1
		}

		return b <= minCostSplitBucket
	})

	node.InitInterior(
		uint8(dim),
		b.buildUpperSAH(treeletRoots, start, mid, totalNodes),
		b.buildUpperSAH(treeletRoots, mid, end, totalNodes),
	)
	return node
}

func (b *BVH) flattenBVHTree(node *BVHBuildNode, offset *uint64) uint64 {
	linearNode := &b.nodes[*offset]
	linearNode.bounds = node.bounds

	myOffset := *offset
	*offset++

	if node.nPrimitives > 0 {
		linearNode.primitiveOffset = uint64(node.firstPrimOffset)
		linearNode.nPrimitives = uint64(node.nPrimitives)
	} else {
		// create interior flattened BVH node
		linearNode.axis = node.splitAxis
		linearNode.nPrimitives = 0
		b.flattenBVHTree(node.children[0], offset)
		linearNode.secondChildOffset = b.flattenBVHTree(node.children[1], offset)
	}

	return myOffset
}

func (b *BVH) WorldBound() *pbrt.Bounds3 {
	if len(b.nodes) == 0 {
		return &pbrt.Bounds3{}
	}
	return &b.nodes[0].bounds
}
func (b *BVH) Intersect(ray *pbrt.Ray, si *pbrt.SurfaceInteraction) bool {
	if len(b.nodes) == 0 {
		return false
	}

	hit := false
	invDir := &pbrt.Vector3f{1 / ray.Direction.X, 1 / ray.Direction.Y, 1 / ray.Direction.Z}
	directionIsNegative := invDir.IsNegative()

	// Follow ray through BVH nodes to find primitive intersections
	var toVisitOffset, currentNodeIndex uint64
	nodesToVisit := [64]uint64{}

	for {
		node := b.nodes[currentNodeIndex]
		// check ray against BVH node
		if node.bounds.IntersectP(ray, invDir, directionIsNegative) {
			if node.nPrimitives > 0 {
				// intersect ray with primitives in leaf BVH node
				for i := uint64(0); i < node.nPrimitives; i++ {
					if b.primitives[node.primitiveOffset+i].Intersect(ray, si) {
						hit = true
					}
				}
				if toVisitOffset == 0 {
					break
				}
				toVisitOffset--
				currentNodeIndex = nodesToVisit[toVisitOffset]
			} else {
				// put far BVH node on nodesToVisit stack, advance to near node
				if directionIsNegative[node.axis] == 1 {
					nodesToVisit[toVisitOffset] = currentNodeIndex + 1
					currentNodeIndex = node.secondChildOffset
					toVisitOffset++
				} else {
					nodesToVisit[toVisitOffset] = node.secondChildOffset
					currentNodeIndex = currentNodeIndex + 1
					toVisitOffset++
				}
			}

		} else {
			if toVisitOffset == 0 {
				break
			}

			toVisitOffset--
			currentNodeIndex = nodesToVisit[toVisitOffset]
		}
	}

	return hit
}
func (b *BVH) IntersectP(ray *pbrt.Ray) bool {
	if len(b.nodes) == 0 {
		return false
	}

	invDir := &pbrt.Vector3f{1 / ray.Direction.X, 1 / ray.Direction.Y, 1 / ray.Direction.Z}
	directionIsNegative := invDir.IsNegative()

	// Follow ray through BVH nodes to find primitive intersections
	var toVisitOffset, currentNodeIndex uint64
	nodesToVisit := [64]uint64{}

	for {
		node := b.nodes[currentNodeIndex]
		// check ray against BVH node
		if node.bounds.IntersectP(ray, invDir, directionIsNegative) {
			if node.nPrimitives > 0 {
				// intersect ray with primitives in leaf BVH node
				for i := uint64(0); i < node.nPrimitives; i++ {
					if b.primitives[node.primitiveOffset+i].IntersectP(ray) {
						return true
					}
				}
				if toVisitOffset == 0 {
					return false
				}
				toVisitOffset--
				currentNodeIndex = nodesToVisit[toVisitOffset]
			} else {
				// put far BVH node on nodesToVisit stack, advance to near node
				if directionIsNegative[node.axis] == 1 {
					nodesToVisit[toVisitOffset] = currentNodeIndex + 1
					currentNodeIndex = node.secondChildOffset
					toVisitOffset++
				} else {
					nodesToVisit[toVisitOffset] = node.secondChildOffset
					currentNodeIndex++
					toVisitOffset++
				}
			}

		} else {
			if toVisitOffset == 0 {
				return false
			}

			toVisitOffset--
			currentNodeIndex = nodesToVisit[toVisitOffset]
			continue
		}

	}
}

func (b *BVH) GetAreaLight() pbrt.AreaLighter {
	log.Panic("Aggregate.GetAreaLight called; should have gone to Geometric Primitiver")
	return nil
}

func (b *BVH) GetMaterial() pbrt.Material {
	log.Panic("Aggregate.GetMaterial called; should have gone to GeometricPrimitive")
	return nil
}

func (b *BVH) ComputeScatteringFunctions(si *pbrt.SurfaceInteraction, mode pbrt.TransportMode, allowMultipleLobes bool) {
	log.Panic("aggregate.ComputeScatteringFunctions called; should have gone to GeometricPrimitive")
}
