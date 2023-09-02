package geo

import (
	"encoding/binary"
	"errors"
	"fmt"
	"log"
	"math"
	"sort"
	"strconv"
	"strings"
)

var (
	ErrInvalidCoordinates = errors.New("invalid coordinates")
)

type Pair [2]float64
type Rect [2]Pair

var (
	//alameda = Pair{AlaLat, AlaLon}
	//longitudeKilometerPerLatitude [91]float64  // lookup table of longitude to Km per each degree latitude
	lonKmLookup [901]float64 // lookup table of longitude to Km per each degree latitude
)

// GeoType for coordinates with slightly less accuracy
// a float32 has 7 digits of precision, which is within ~11cm
//
// Virtually all geo data is "close enough" using this,
// and for data that heavily comprises these points,
// one can reduce memory footprint by one half (for points)
type GeoType float64

// GeoPoints provides abstraction for slices of data with coordinates
type GeoPoints interface {
	IndexPoint(int) Point
	Len() int
}

const (
	// DegreeToKilometer is a "constant" for latitude but varies for longitude
	DegreeToKilometer     = 111.111 //111.321
	MilesToKilometer      = 1.609344
	EarthRadiusInKM       = 6371.1 // 6378.1
	SquareKmPerSquareMile = 0.386102
)

func init() {
	for i := 0; i < len(lonKmLookup); i++ {
		lat := float64(i) / 10.0
		lonKmLookup[i] = LonKilos(lat)
	}
}

const Radian = math.Pi / 180.0

func deg2rad(d float64) float64 {
	return d * Radian
}

// LookupLonKmPerLat returns the ratio of kilometers to degrees longitude
// at the given latitude.
//
// Accuracy is with 1% under 80 degrees, which is good enough for most work
func LookupLonKmPerLat(lat float64) float64 {
	idx := int(lat * 10)
	return lonKmLookup[idx]
}

// LonKilos returns the kilometers per degree longitude at the given latitude
func LonKilos(lat float64) float64 {
	return math.Cos(deg2rad(lat)) * DegreeToKilometer
}

// LongitudeKilometerDegrees converts the distance given
// in kilometers at that latitude to degrees
func LongitudeKilometerDegrees(lat, kilometers float64) float64 {
	return kilometers / LonKilos(lat)
}

// LongitudeKilometers returns the distance of the degrees lon, at latitude lat
func LongitudeKilometers(lat, lon float64) float64 {
	return LonKilos(lat) * lon
}

// Expand returns the a box with a radius in kM for
func Expand(lat, lon, radiusKM float64) Rect {
	latx := radiusKM / lon
	lonx := LongitudeKilometerDegrees(lat, radiusKM)
	return Rect{
		Pair{lat - latx, lon - lonx},
		Pair{lat + latx, lon + lonx},
	}
}

// Distance returns the distance in kM between 2 geographic points
// It uses the Haversine formula for spherical calculations
func Distance(lat1, lon1, lat2, lon2 float64) float64 {
	dlat1 := deg2rad(lat1)
	dlon1 := deg2rad(lon1)
	dlat2 := deg2rad(lat2)
	dlon2 := deg2rad(lon2)

	return math.Acos(math.Sin(dlat1)*math.Sin(dlat2)+math.Cos(dlat1)*math.Cos(dlat2)*math.Cos(dlon2-dlon1)) * EarthRadiusInKM
}
func Distance32(lat1, lon1, lat2, lon2 float32) float64 {
	return Distance(float64(lat1), float64(lon1), float64(lat2), float64(lon2))
}

func DistanceGeoType(lat1, lon1, lat2, lon2 GeoType) float64 {
	return Distance(float64(lat1), float64(lon1), float64(lat2), float64(lon2))
}

// ApproximateDistanceGeo returns the approximate distance between 2 points
// It uses the pythagarean distance calc which is meant for 2d operations
// but is "good enough" for shorter distances (which we primarily care about)
// It is about 7 times faster than the "proper way"
func ApproximateDistance(lat1, lon1, lat2, lon2 float64) float64 {
	//	lonDegreeKm := LookupLonKmPerLatInt(round(float64(lat2)))
	//lonDegreeKm := LookupLonKmPerLatInt(int(lat2))
	lonDegreeKm := LookupLonKmPerLat(lat2)
	a := float64(lat2-lat1) * DegreeToKilometer
	b := float64(lon2-lon1) * lonDegreeKm
	return math.Sqrt(math.Pow(a, 2) + math.Pow(b, 2))
}

// ApproximateDistanceGeo returns the approximate distance between 2 points
// It uses the pythagarean distance calc which is meant for 2d operations
// but is "good enough" for shorter distances (which we primarily care about)
// It is about 7 times faster than the "proper way"
func ApproximateDistanceGeo(lat1, lon1, lat2, lon2 GeoType) float64 {
	//	lonDegreeKm := LookupLonKmPerLatInt(round(float64(lat2)))
	lonDegreeKm := LookupLonKmPerLat(float64(lat1)) //LookupLonKmPerLatInt(int(lat2))
	a := float64(lat2-lat1) * DegreeToKilometer
	b := float64(lon2-lon1) * lonDegreeKm
	return math.Sqrt(math.Pow(a, 2) + math.Pow(b, 2))
}

func GeoPoint(lat, lon float64) Point {
	return Point{GeoType(lat), GeoType(lon)}
}

func AreaInKm(lat1, lon1, lat2, lon2 float64) float64 {
	delta := math.Abs(lon2 - lon1)
	d1 := LongitudeKilometers(lat1, delta)
	d2 := LongitudeKilometers(lat2, delta)
	avg := math.Abs(d2+d1) / 2
	h := math.Abs(lat2-lat1) * DegreeToKilometer
	a := h * avg
	return a
}

func SquareKmInMiles(k float64) float64 {
	return k * SquareKmPerSquareMile
}

func AreaInMiles(lat1, lon1, lat2, lon2 float64) float64 {
	return SquareKmInMiles(AreaInKm(lat1, lon1, lat2, lon2))
}

func Coords(query string) ([]GeoType, error) {
	parts := strings.Split(query, ",")
	if len(parts) != 4 {
		parts = strings.Split(query, "/")
		if len(parts) != 4 {
			return nil, ErrInvalidCoordinates
		}
	}
	coords, err := geos(parts...)
	if err != nil {
		return nil, fmt.Errorf("parse failure (%v): %w", err, ErrInvalidCoordinates)
	}
	if err := coordCheck(coords...); err != nil {
		for i := len(coords); i < 4; i++ {
			coords = append(coords, 0)
		}
		return nil, fmt.Errorf(badPrefix+" %w", err)
	}
	return coords, nil
}

type Point struct {
	Lat, Lon GeoType
}

// Less returns true if it is less than the given point
func (p Point) Less(x Point) bool {
	if p.Lat < x.Lat {
		return true
	} else if p.Lat > x.Lat {
		return false
	} else {
		// lon is secondary sort
		return p.Lon < x.Lon
	}
}

func (p Point) Equal(x Point) bool {
	return p.Lat == x.Lat && p.Lon == x.Lon
}

// Label returns a consistent string representation of the coordinates
func (p Point) Label() string {
	return fmt.Sprintf("%010.5f_%010.5f", p.Lat, p.Lon)
}

func (p Point) Distance(x Point) float64 {
	return DistanceGeoType(p.Lat, p.Lon, x.Lat, x.Lon)
}

func (p Point) Approximately(x Point) float64 {
	return ApproximateDistanceGeo(p.Lat, p.Lon, x.Lat, x.Lon)
}

// AreaInRange64 is like AreaInRange but using float64
func AreaInRange64(pt Pair, distance float64) Rect {
	lat := pt[0]
	lon := pt[1]
	deltaLat := (distance / DegreeToKilometer)
	deltaLon := (LongitudeKilometerDegrees(float64(lat), distance))
	min := Pair{lat - deltaLat, lon - deltaLon}
	max := Pair{lat + deltaLat, lon + deltaLon}
	return Rect{min, max}
}

// Rightmost searches a list for the next smallest value in the list
// per https://en.wikipedia.org/wiki/Binary_search_algorithm
func Rightmost(n int, fn func(int) bool) int {
	L := 0
	R := n
	for L < R {
		// m := floor((L + R) / 2)
		m := (L + R) / 2
		fmt.Printf("m:%d L:%d R:%d\n", m, L, R)
		// if A[m] > T {
		if fn(m) {
			R = m
		} else {
			L = m + 1
		}
	}
	return R - 1
}

// Closest searches for a matching point within the distance (in Km)
// of the specified point.
// It returns the index of the closest point and the distance from the target
// If nothing is found, it returns the Len() of the points list and -1 distance
//
// NOTE: this is an adaptation of Bestest, but distances are approximated to
//
//	minimize computational load
//
// TODO: the len return is in line w/ Go sort.Search, but perhaps -1 would be better?
// TODO part too: use distance func to share same routine w/ approx and haversine calcs?
func Closest(g GeoPoints, pt Point, deltaKm float64) (int, float64) {
	// Do a binary search to find the "closest" match

	// The point found is not guaranteed to actually be
	// the shortest distance, as it finds the first point
	// that is equal to or *greater* than what is searched for

	// It's possible for the first hit to be significantly further away,
	// whereas an entry before it, while less than the point, is still
	// closer.

	// The final confounding factor is that the data is order by
	// latitude *then* longitude, so the following point could be
	// very close in latitude, the longitude could be much greater.
	// A subsequent point could be 0.000001 degrees latitude
	// further (0.11 m), but have the longitude diff be much less

	x := sort.Search(g.Len(), func(i int) bool {
		h := g.IndexPoint(i)
		return pt.Less(h)
	})

	// did search fail?
	if x == g.Len() {
		return x, -1 //math.MaxFloat64//closest
	}

	// so we either came in exactly on target (not likely),
	// or somewhat past it.
	// Start by going backwards as we already know we probably overshot
	// if we are greater than this then there's no possibility
	// we will be in range continuing on

	// To minimize work done (calculating distance),
	// calculate the furthest away directly by latidude only,
	// as that is (effectively) invariant
	minLat := pt.Lat - GeoType(deltaKm/DegreeToKilometer)

	best := x //g.Len()

	//closest := deltaKm + 0.0001 // ensure we have something to best
	counter := 0

	// our first hit is guaranteed to be equal to or *greater* than our
	// requested point.
	//
	// we have to check both above and below the point in question to see
	// which has the closed hit
	this := g.IndexPoint(x)
	dist := this.Approximately(pt)
	closest := dist
	debugf("first hit for %v: %v -- %6d/%6d (%f)", pt, this, x, g.Len(), dist)

	lonKmPerDegree := LookupLonKmPerLat(float64(pt.Lat)) //LookupLonKmPerLatInt(int(pt.Lat))
	deltaLon := GeoType(closest / lonKmPerDegree)
	lonOutside := func(lon GeoType) bool {
		maxLon := pt.Lon + deltaLon
		minLon := pt.Lon - deltaLon
		return lon < minLon || lon > maxLon
	}

	// work backwards first, as we likely overshot our target
	for i := x - 1; i > 0; i-- {
		counter++
		this = g.IndexPoint(i)
		if this.Lat < minLat {
			//debugf("%v exceeded minimum possible lat: %v", this, minLat)
			break
		}
		if lonOutside(this.Lon) {
			//debugf("below lon outside: %v", this)
			continue
		}
		if dist := pt.Approximately(this); dist < closest {
			closest = dist
			best = i
			minLat = pt.Lat - GeoType(closest/DegreeToKilometer)
			deltaLon = GeoType(closest / lonKmPerDegree)
			//debugf("(%d) MINLAT: %f", counter, minLat)
		}
	}
	/*
	   so we're within range... now we keep looking for anything closer past
	   the initial hit, which will not last long, as we're moving *away* from our destination,
	   and the only likely improvement is if the longitude was off and we
	   find a hit that's closer

	   The closest *possible* will be the same lon but directly above our point,
	   so we calculate what that lat is for the current minimal distance
	   and once we hit pass that lat we know nothing can be closer and we're
	   done with that sweep.
	*/
	maxLat := pt.Lat + GeoType(closest/DegreeToKilometer)
	for i := x + 1; i < g.Len(); i++ {
		counter++
		this = g.IndexPoint(i)
		if this.Lat > maxLat {
			//debugf("%v exceeds max lat of %v", this, maxLat)
			break
		}
		if lonOutside(this.Lon) {
			//debugf("above lon outside: %v", this)
			continue
		}
		if dist := this.Approximately(pt); dist < closest {
			best = i
			closest = dist
			maxLat = pt.Lat + GeoType(dist/DegreeToKilometer)
		}
	}
	//debugf("Examined %d records", counter)

	return best, closest
}

func nextClosest(gp GeoPoints, idx int, maxLat float64, pt Point) (int, float64) {
	dist := math.MaxFloat64
	same := true
	best := idx
	for idx < gp.Len() {
		this := gp.IndexPoint(idx)
		if this.Lat > pt.Lat && same {
			same = false
			log.Printf("moved past original lat:%f to %f", pt.Lat, this.Lat)
		}
		if this.Equal(pt) {
			log.Printf("hit the jackpot at %d", idx)
			return idx, 0.0
		}
		far := this.Distance(pt)
		if far < dist {
			dist = far
			best = idx
			log.Printf("[%11d:%8.3f:%v] next closer", best, dist, this)
			idx++
			continue
		}
		delta := float64(this.Lat - pt.Lat)
		if delta > maxLat {
			log.Printf("delta %f exceeds max lat diff %f at pt %v", delta, maxLat, this)
			return idx, dist
		}
		if float64(this.Lat) >= maxLat {
			return best, dist
		}
		// now moving away, so our last hit was best
		// return idx - 1, dist
		idx++
	}
	log.Println("shouldn't ever hit here, right")
	return gp.Len(), math.MaxFloat64

}

func priorClosest(g GeoPoints, idx int, firstLat, maxLat float64, pt Point) (int, float64) {
	// first move to beginning of previous lat series (there may be a bunch)
	// idx is 1 less than was found originally
	log.Printf("prior lat to %.12f (for %v)", firstLat, pt)
	if gl := g.Len(); idx >= gl {
		log.Printf("oops idx:%d exceeds len:%d", idx, gl)
		return gl, 999999999999999.9
	} else {
		log.Printf("mkay idx:%d versus len:%d", idx, gl)
	}
	previous := g.IndexPoint(idx)
	size := idx + 1
	x := sort.Search(size, func(i int) bool {
		h := g.IndexPoint(i)
		return pt.Less(h)
	})
	if x == size {
		log.Printf("failed to find beginning of previous")
		return g.Len(), -1
	}
	log.Printf("first rec for prior lat %f -- %d (%v)", firstLat, x, g.IndexPoint(x))

	dist := math.MaxFloat64
	same := true
	for idx >= 0 {
		this := g.IndexPoint(idx)
		log.Printf("scrollback %d -- %v", idx, this)
		if this.Lat == pt.Lat && this.Lon == pt.Lon {
			log.Printf("hot damn a back search hits at %d", idx)
			return idx, 0
		}
		if this.Lat != previous.Lat && same {
			same = false
			log.Printf("at idx %d we slipped from %f to %f", idx, previous.Lat, this.Lat)
		}
		far := this.Distance(pt)
		if far < dist {
			dist = far
			idx--
			continue
		}
		// now moving away, so our last hit was best
		return idx - 1, dist
	}
	log.Println("shouldn't ever hit here, right")
	return g.Len(), math.MaxFloat64

}

// CloseEnough searches for a matching point within the distance (in Km)
// of the specified point.
// It returns the index of the closest point and the distance from the target
// If nothing is found, it returns the Len() of the points list and -1 distance
//
// Possible search scenarios:
/*
	1. Exact match on lat/lon
	2. Exact match on lat, closest lon is in range of that lon
	3. Exact match on lat, but cloest pt is subsequent lat/lon (i.e., only lon in that lat range are too far)
	4. Exact match on lat, but closest point on previous lat/lon
	5. Exact match on lat, but closest point is on *prior to previous* lat
	6. Lat is past point, but that or subsequent records has closest point
	7. Lat is past point but prior lat has closest record
	8. Lat is passt point but record previous to prior lat is closest
*/

// NOTE: this is an adaptation of Bestest, but distances are approximated to
//
//	minimize computational load
//
// TODO: the len return is in line w/ Go sort.Search, but perhaps -1 would be better?
// TODO part too: use distance func to share same routine w/ approx and haversine calcs?
func CloseEnough(g GeoPoints, pt Point, deltaKm float64) (int, float64) {
	// Do a binary search to find the "closest" match

	// The point found is not guaranteed to actually be
	// the shortest distance, as it finds the first point
	// that is equal to or *greater* than what is searched for

	// It's possible for the first hit to be significantly further away,
	// whereas an entry before it, while less than the point, is still
	// closer.

	// The final confounding factor is that the data is order by
	// latitude *then* longitude, so the following point could be
	// very close in latitude, the longitude could be much greater.
	// A subsequent point could be 0.000001 degrees latitude
	// further (0.11 m), but have the longitude diff be much less

	// To minimize work done (calculating distance),
	// calculate the furthest away directly by latidude only,
	// as that is (effectively) invariant
	maxDeltaLat := pt.Lat - GeoType(deltaKm/DegreeToKilometer)

	x := sort.Search(g.Len(), func(i int) bool {
		h := g.IndexPoint(i)
		// return pt.Less(h)
		return pt.Lat <= h.Lat
	})

	// did search fail?
	if x == g.Len() {
		return x, -1 //math.MaxFloat64//closest
	}

	this := g.IndexPoint(x)
	// if we're lucky, maybe we're right on target and our search is complete!
	if this.Lat == pt.Lat && this.Lon == pt.Lon {
		log.Println("that's a bingo:", x)
		return x, 0.0
	}

	latDelta := this.Lat - pt.Lat
	dist := pt.Distance(this)
	log.Printf("[%11d:%8.3f:%v] OUR FIRST HIT: lat delta:%f", x, dist, this, latDelta)
	for i := 0; i < 8; i++ {
		y := x - i
		if y < 0 {
			break
		}
		that := g.IndexPoint(y)
		if d := pt.Distance(that); d < dist {
			x = y
			dist = d
			log.Printf("[%11d:%8.3f:%v] BEFORE FIRST HIT", x, dist, this)
			// log.Printf("BEFORE HIT: %d -- %v (%f)", y, that, pt.Distance(that))
		} else {
			log.Printf("[%11d:%8.3f:%v] NOPES BEFORE", y, d, that)
			// log.Printf("NOPES BEFORE: %d -- %v (%f)", y, that, pt.Distance(that))
			break
		}
	}
	// ok, not a full on match; do we at least have an exact match on latitude?
	// if this.Lat == pt.Lat {
	// if latDelta == 0 {
	latRange := float64(maxDeltaLat)
	lonIdx, nd := nextClosest(g, x, latRange, pt)
	// log.Printf("[%11d:%8.3f:%v] OUR FIRST HIT: lat delta:%f", x, dist, this, latDelta)
	// log.Printf("next closest idx:%d dist:%f pt:%v", lonIdx, dist, g.IndexPoint(lonIdx))
	// if lonIdx != g.Len() && dist == 0.0 {
	// 	return lonIdx, dist
	// }
	if nd == 0 {
		return lonIdx, 0
	}
	if lonIdx != g.Len() && nd < dist {
		log.Printf("[%11d:%8.3f:%v] NEXT CLOSEST YES!", lonIdx, nd, this)
		dist = nd
		x = lonIdx
	} else {
		log.Printf("[%11d:%8.3f:%v] NEXT CLOSEST nah", lonIdx, nd, this)
	}
	if x > 0 {
		pIdx, pDist := priorClosest(g, x-1, float64(this.Lat), latRange, pt)
		if pIdx >= g.Len() {
			log.Printf("WHOA prior faile idx: %d/%d", pIdx, g.Len())
		} else {
			log.Printf("prev closest idx:%d dist:%f pt:%v", pIdx, pDist, g.IndexPoint(pIdx))
		}
		// NOTE: this is a temp hack for initial
		if pIdx != g.Len() && pDist < dist {
			dist = pDist
			x = pIdx
			// return pIdx, dist
		}
	}
	if true {
		return x, dist
	}
	// }
	// log.Println("oh fuckitall we are in the mud")
	// so we either came in exactly on target (not likely),
	// or somewhat past it.
	// Start by going backwards as we already know we probably overshot
	// if we are greater than this then there's no possibility
	// we will be in range continuing on

	best := x //g.Len()

	//closest := deltaKm + 0.0001 // ensure we have something to best
	counter := 0

	// our first hit is guaranteed to be equal to or *greater* than our
	// requested point.
	//
	// we have to check both above and below the point in question to see
	// which has the closed hit
	dist = this.Approximately(pt)
	closest := dist
	debugf("first hit for %v: %v -- %6d/%6d (%f)", pt, this, x, g.Len(), dist)

	lonKmPerDegree := LookupLonKmPerLat(float64(pt.Lat)) //LookupLonKmPerLatInt(int(pt.Lat))
	deltaLon := GeoType(closest / lonKmPerDegree)
	lonOutside := func(lon GeoType) bool {
		maxLon := pt.Lon + deltaLon
		minLon := pt.Lon - deltaLon
		return lon < minLon || lon > maxLon
	}

	// work backwards first, as we likely overshot our target
	for i := x - 1; i > 0; i-- {
		counter++
		this = g.IndexPoint(i)
		if this.Lat < maxDeltaLat {
			//debugf("%v exceeded minimum possible lat: %v", this, minLat)
			break
		}
		if lonOutside(this.Lon) {
			//debugf("below lon outside: %v", this)
			continue
		}
		if dist := pt.Approximately(this); dist < closest {
			closest = dist
			best = i
			maxDeltaLat = pt.Lat - GeoType(closest/DegreeToKilometer)
			deltaLon = GeoType(closest / lonKmPerDegree)
			//debugf("(%d) MINLAT: %f", counter, minLat)
		}
	}
	/*
	   so we're within range... now we keep looking for anything closer past
	   the initial hit, which will not last long, as we're moving *away* from our destination,
	   and the only likely improvement is if the longitude was off and we
	   find a hit that's closer

	   The closest *possible* will be the same lon but directly above our point,
	   so we calculate what that lat is for the current minimal distance
	   and once we hit pass that lat we know nothing can be closer and we're
	   done with that sweep.
	*/
	maxLat := pt.Lat + GeoType(closest/DegreeToKilometer)
	for i := x + 1; i < g.Len(); i++ {
		counter++
		this = g.IndexPoint(i)
		if this.Lat > maxLat {
			//debugf("%v exceeds max lat of %v", this, maxLat)
			break
		}
		if lonOutside(this.Lon) {
			//debugf("above lon outside: %v", this)
			continue
		}
		if dist := this.Approximately(pt); dist < closest {
			best = i
			closest = dist
			maxLat = pt.Lat + GeoType(dist/DegreeToKilometer)
		}
	}
	//debugf("Examined %d records", counter)

	return best, closest
}

func between(check, min, max GeoType) bool {
	return min <= check && check <= max
}

func Within(lat, lon, minLat, minLon, maxLat, maxLon GeoType) bool {
	return between(lat, minLat, maxLat) && between(lon, minLon, maxLon)
}

func geos(ss ...string) ([]GeoType, error) {
	ff := make([]GeoType, 0, len(ss))
	for i, s := range ss {
		f, err := strconv.ParseFloat(s, 64)
		if err != nil {
			return nil, fmt.Errorf("(%d/%d): %q ain't a number: %w", i, len(ss), s, err)
		}
		ff = append(ff, GeoType(f))
	}
	return ff, nil
}

const badPrefix = `bad coordinates --`

func coordCheck(coord ...GeoType) error {
	if len(coord) < 4 {
		return fmt.Errorf("requires 4 coordinates")
	}
	switch {
	// affirm positive lats, negative lons
	case coord[0] == 0 && coord[1] == 0 && coord[2] == 0 && coord[3] == 0:
		return fmt.Errorf(badPrefix + " zero coordinates are not allowed")
	case coord[0] < 0:
		return fmt.Errorf(badPrefix+" latitude %.5f cannot be negative", coord[0])
	case coord[1] > 0:
		return fmt.Errorf(badPrefix+" longitude %.5f must be negative", coord[1])
	case coord[2] < 0:
		return fmt.Errorf(badPrefix+" latitude %.5f cannot be negative", coord[2])
	case coord[3] > 0:
		return fmt.Errorf(badPrefix+" longitude %.5f cannot be negative", coord[3])

	case coord[0] > coord[2]:
		return fmt.Errorf(badPrefix+" latitude %.5f must be less than %.5f", coord[0], coord[2])
		/*
			case coord[0] < AllowedMinLat:
				return fmt.Errorf(badPrefix+" latitude %.5f must be at least %.5f", coord[0], AllowedMinLat)
			case coord[0] > AllowedMaxLat:
				return fmt.Errorf(badPrefix+" latitude %.5f must be less than %.5f", coord[0], AllowedMaxLat)
		*/
	case coord[1] > coord[3]:
		return fmt.Errorf(badPrefix+" longitude %.5f must be less than %.5f", coord[1], coord[3])
		/*
			case coord[1] < AllowedMinLon:
				return fmt.Errorf(badPrefix+" longitude %.5f must be at least %.5f", coord[1], AllowedMinLon)
			case coord[1] > AllowedMaxLon:
				return fmt.Errorf(badPrefix+" longitude %.5f must be less than %.5f", coord[1], AllowedMaxLon)

			case coord[2] > AllowedMaxLat:
				return fmt.Errorf(badPrefix+" latitude %.5f must be less than %.5f", coord[2], AllowedMaxLat)

			case coord[3] > AllowedMaxLon:
				return fmt.Errorf(badPrefix+" longitude %.5f must be less than %.5f", coord[3], AllowedMaxLon)
		*/
	}
	return nil
}

func QueryCoords(s string) (Pair, error) {
	parts := strings.Split(s, ",")
	if len(parts) < 2 {
		parts = strings.Split(s, "/")
		if len(parts) < 2 {
			return Pair{}, ErrInvalidCoordinates
		}
	}
	lat, err := strconv.ParseFloat(parts[0], 32)
	if err != nil {
		return Pair{}, fmt.Errorf("invalid latitude %q -- %w", parts[0], ErrInvalidCoordinates)
	}
	lon, err := strconv.ParseFloat(parts[1], 32)
	if err != nil {
		return Pair{}, fmt.Errorf("invalid longitude %q -- %w", parts[0], ErrInvalidCoordinates)
	}
	return Pair{lat, lon}, nil
}

func QueryPoint(s string) (Point, error) {
	pt, err := QueryCoords(s)
	return Point{GeoType(pt[0]), GeoType(pt[1])}, err
}

// Bestest searches for a matching point within the distance (in Km)
// of the specified point.
// It returns the index of the closest point and the distance from the target
// If nothing is found, it returns the Len() of the points list and -1 distance
//
// NOTE: this is an adaptation of Bestest, but distances are approximated to
//
//	minimize computational load
//
// TODO: the len return is in line w/ Go sort.Search, but perhaps -1 would be better?
// TODO part too: use distance func to share same routine w/ approx and haversine calcs?
func Bestest(g GeoPoints, pt Point, deltaKm float64) (int, float64) {
	// Do a binary search to find the "closest" match

	// The point found is not guaranteed to actually be
	// the shortest distance, as it finds the first point
	// that is equal to or *greater* than what is searched for

	// It's possible for the first hit to be significantly further away,
	// whereas an entry before it, while less than the point, is still
	// closer.

	// The final confounding factor is that the data is order by
	// latitude *then* longitude, so the following point could be
	// very close in latitude, the longitude could be much greater.
	// A subsequent point could be 0.000001 degrees latitude
	// further (0.11 m), but have the longitude diff be much less

	x := sort.Search(g.Len(), func(i int) bool {
		h := g.IndexPoint(i)
		return pt.Less(h)
	})

	// did search fail?
	if x == g.Len() {
		return x, -1 //math.MaxFloat64//closest
	}

	// so we either came in exactly on target (not likely),
	// or somewhat past it.
	// Start by going backwards as we already know we probably overshot
	// if we are greater than this then there's no possibility
	// we will be in range continuing on

	// To minimize work done (calculating distance),
	// calculate the furthest away directly by latidude only,
	// as that is (effectively) invariant
	minLat := pt.Lat - GeoType(deltaKm/DegreeToKilometer)

	best := g.Len()

	closest := deltaKm + 0.0001 // ensure we have something to best
	counter := 0

	// our first hit is guaranteed to be equal to or *greater* than our
	// requested point.
	//
	// we have to check both above and below the point in question to see
	// which has the closed hit
	this := g.IndexPoint(x)
	dist := this.Distance(pt)
	debugf("first hit: %6d/%6d (%f)", x, g.Len(), dist)
	if dist < closest {
		closest = dist
		best = x
	}
	debugf("(%d) PT.LAT:%f MINLAT:%f", counter, this.Lat, minLat)

	// only check if lon is in range as well
	/*
		e.g, 100km per degree lon
		10 km range
		10 km / 100km/degree = 0.1 degree delta lon
	*/
	lonKmPerDegree := LonKilos(float64(pt.Lat)) //LookupLonKmPerLatInt(int(pt.Lat))
	deltaLon := GeoType(closest / lonKmPerDegree)
	lonOutside := func(lon GeoType) bool {
		maxLon := pt.Lon + deltaLon
		minLon := pt.Lon - deltaLon
		return lon < minLon || lon > maxLon
	}

	// work backwards first, as we likely overshot our target
	for i := x - 1; i > 0; i-- {
		counter++
		this = g.IndexPoint(i)
		if this.Lat < minLat {
			debugf("%v exceeded minimum possible lat: %v", this, minLat)
			break
		}
		if lonOutside(this.Lon) {
			continue
		}
		if dist := pt.Distance(this); dist < closest {
			closest = dist
			best = i
			minLat = pt.Lat - GeoType(closest/DegreeToKilometer)
			deltaLon = GeoType(closest / lonKmPerDegree)
			debugf("(%d) MINLAT: %f", counter, minLat)
		}
	}
	/*
	   so we're within range... now we keep looking for anything closer past
	   the initial hit, which will not last long, as we're moving *away* from our destination,
	   and the only likely improvement is if the longitude was off and we
	   find a hit that's closer

	   The closest *possible* will be the same lon but directly above our point,
	   so we calculate what that lat is for the current minimal distance
	   and once we hit pass that lat we know nothing can be closer and we're
	   done with that sweep.
	*/
	maxLat := pt.Lat + GeoType(closest/DegreeToKilometer)
	for i := x + 1; i < g.Len(); i++ {
		counter++
		this = g.IndexPoint(i)
		if this.Lat > maxLat {
			debugf("%v exceeds max lat of %v", this, maxLat)
			break
		}
		if lonOutside(this.Lon) {
			continue
		}
		if dist := this.Distance(pt); dist < closest {
			best = i
			closest = dist
			maxLat = pt.Lat + GeoType(dist/DegreeToKilometer)
		}
	}
	debugf("Examined %d records", counter)

	return best, closest
}

func ToGeoType(value interface{}) (GeoType, error) {
	switch value := value.(type) {
	case float32:
		return GeoType(value), nil
	case float64:
		return GeoType(value), nil
	case int:
		return GeoType(value), nil
	case int32:
		return GeoType(value), nil
	case int64:
		return GeoType(value), nil
	case string:
		f, err := strconv.ParseFloat(value, 32)
		return GeoType(f), err
	}
	return 0, fmt.Errorf("%v is un unsupported type: %T", value, value)
}

func DecodePoint(buf []byte) Point {
	Lat := GeoType(math.Float32frombits(binary.LittleEndian.Uint32(buf)))
	Lon := GeoType(math.Float32frombits(binary.LittleEndian.Uint32(buf[4:])))
	return Point{Lat, Lon}
}

func DecodePair(buf []byte) Pair {
	Lat := math.Float64frombits(binary.LittleEndian.Uint64(buf))
	Lon := math.Float64frombits(binary.LittleEndian.Uint64(buf[8:]))
	return Pair{Lat, Lon}
}

func Exists(it *Iter, pt Point) int {
	for i := 0; i < it.Len(); i++ {
		this := it.IndexPoint(i)
		if pt.Equal(this) {
			return i
		}
		if pt.Less(this) {
			log.Printf("giving up on pt %v at %d (exceeded by %v)", pt, i, this)
			break
		}
	}
	return it.Len()
}

func FuckMe() {
	log.Println("FUCK ME")
}
