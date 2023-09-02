package geo

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
)

const (
	// handy for me
	AlaLat, AlaLon = 37.7703358, -122.2569864
	AlaGeo         = 6001 // GeoID of Alameda

	PortLat, PortLon = 45.557966, -122.867614 // 45.5623746, -122.8675477

	// Houston City Hall
	HouLat, HouLon = 29.6223432, -95.5079226 //29.6223432, -95.5079226

	// SF Transit Center
	SFLat, SFLon = 37.789831923849874, -122.3953253191605

	// Zephyr Cove, Lake Tahoe
	ZepLat, ZepLon = 39.00162305522508, -119.95229429908939

	// https://www.nhc.noaa.gov/gccalc.shtml
	// From SF Transit Center to Zephyr Cove
	SFtoZep  = 252   // km, per the noaa calc service
	SFtoAla  = 12.55 // per Google maps
	SFtoHou  = 2635
	AlaToPDX = 867
)

func TestGeoTypeAccuracy(t *testing.T) {
	roundTrip := float64(GeoType(AlaLat))
	// returns diff of -0.000002
	// which is ~ 0.22m
	delta := AlaLat - roundTrip
	distance := (delta * DegreeToKilometer) * 1000.0
	t.Logf("lost in translation: %f -- %fm", delta, distance)
}

func TestDistance(t *testing.T) {
	const kmExpected = 868.05
	dist := Distance(AlaLat, AlaLon, PortLat, PortLon)
	diff := dist - kmExpected
	t.Logf("delta: %f", diff)
}

func TestExpand(t *testing.T) {
	box := Expand(AlaLat, AlaLon, 1.0)
	area := AreaInKm(box[0][0], box[0][1], box[1][0], box[1][1])
	t.Logf("area: %f", area)
}

/*
func TestAreaInRange(t *testing.T) {
	t.Skip("review this")
	const (
		lat = AlaLat
		lon = AlaLon
	)
	pt := Point{lat, lon}
	const (
		distance    = 5.0
		diffAllowed = 0.0001
	)
	area := AreaInRange64(pt, distance)
	dLat := float64(((area.Max.Lat - area.Min.Lat) * DegreeToKilometer) / 2)
	dLon := float64(LongitudeKilometers(lat, float64(area.Max.Lon-area.Min.Lon)) / 2)

	t.Logf("AREA: %v", area)
	t.Logf("Lat: %f, Lon: %f", dLat, dLon)
	assert.LessOrEqual(t, float64(diffAllowed), float64(distance-dLat))
	assert.LessOrEqual(t, float64(diffAllowed), float64(distance-dLon))
}
*/

func TestAreaInRange64(t *testing.T) {
	t.Skip("wtf?")
	const (
		lat = AlaLat
		lon = AlaLon
	)
	pt := Pair{lat, lon}
	const (
		distance    = 5.0
		diffAllowed = 0.0001
	)
	area := AreaInRange64(pt, distance)
	dLat := ((area[1][0] - area[0][0]) * DegreeToKilometer) / 2
	dLon := LongitudeKilometers(lat, float64(area[1][1]-area[0][1])) / 2

	t.Logf("AREA: %v", area)
	t.Logf("Lat: %f, Lon: %f", dLat, dLon)
	assert.LessOrEqual(t, float64(diffAllowed), float64(distance-dLat))
	assert.LessOrEqual(t, float64(diffAllowed), float64(distance-dLon))
}

type testPoints []Point

func (t testPoints) IndexPoint(i int) Point {
	return t[i]
}

func (t testPoints) Len() int {
	return len(t)
}

func (t testPoints) Less(i, j int) bool {
	return t[i].Less(t[j])
}

func (t testPoints) Swap(i, j int) {
	t[i], t[j] = t[j], t[i]
}

/*
type Helper interface {
	Helper()
}
*/

func searchSample(t Helper, include bool) (Point, GeoPoints) {
	t.Helper()
	xLat := GeoType(AlaLat)
	xLon := GeoType(AlaLon)
	center := Point{xLat, xLon}
	var points testPoints
	for lat := GeoType(0.001); lat < 1.0; lat += 0.005 {
		for lon := GeoType(0.0001); lat < 1.0; lat += 0.0001 {
			points = append(points,
				Point{xLat + lat, xLon + lon},
				Point{xLat - lat, xLon + lon},
				Point{xLat + lat, xLon - lon},
				Point{xLat - lat, xLon - lon},
			)
		}
	}
	if include {
		points = append(points, center)
	}
	sort.Sort(points)
	return center, points
}

func TestClosest(t *testing.T) {
	pt, list := searchSample(t, false)
	const deltaKm = 0.1
	now := time.Now()
	i, dist := Closest(list, pt, deltaKm)
	if i == list.Len() {
		t.Fatalf("nada, baby")
	}
	ptx := list.IndexPoint(i)
	t.Logf("I: %d, DIST:%f (%v) TIME:%s", i, dist, ptx, time.Since(now))
}

func TestBestest(t *testing.T) {
	pt, list := searchSample(t, false)
	const deltaKm = 0.1 //1.0 //0.1
	now := time.Now()
	i, dist := Bestest(list, pt, deltaKm)
	t.Logf("I: %d, DIST:%f TIME:%s", i, dist, time.Since(now))
}

func BenchmarkClosest(b *testing.B) {
	pt, list := searchSample(b, false)
	const deltaKm = 0.1
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Closest(list, pt, deltaKm)
	}
}

func BenchmarkBestest(b *testing.B) {
	pt, list := searchSample(b, false)
	const deltaKm = 0.1
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Bestest(list, pt, deltaKm)
	}
}

func BenchmarkDistance(b *testing.B) {
	const lat1, lon1 = AlaLat, AlaLon
	const lat2, lon2 = PortLat, PortLon
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Distance(lat1, lon1, lat2, lon2)
	}
}

func BenchmarkApproximateDistance(b *testing.B) {
	const lat1, lon1 = AlaLat, AlaLon
	const lat2, lon2 = PortLat, PortLon
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ApproximateDistance(lat1, lon1, lat2, lon2)
	}
}

func TestAccuracy(t *testing.T) {
	const (
		pt1Lat, pt1Lon = 47.7690679, -122.2592744
		pt2Lat, pt2Lon = 47.78978486773513, -122.39617851945805
	)
	real := Distance(pt1Lat, pt1Lon, pt2Lat, pt2Lon)
	fake := ApproximateDistanceGeo(pt1Lat, pt1Lon, pt2Lat, pt2Lon)
	off := (math.Abs(fake-real) / real) * 100.0
	t.Logf("REAL:%f FAKE:%f (%.1f)", real, fake, off)
}

func TestAccuracyPart2(t *testing.T) {
	const (
		pt1Lat, pt1Lon = 47.7690679, -122.2592744
		pt2Lat, pt2Lon = 47.78978486773513, -122.39617851945805
		pt3Lat, pt3Lon = 37.7, -122.2
		pt4Lat, pt4Lon = 37.8, -122.3
	)
	pt1 := GeoPoint(pt1Lat, pt1Lon)
	pt2 := GeoPoint(pt2Lat, pt2Lon)
	pt3 := GeoPoint(pt3Lat, pt3Lon)
	pt4 := GeoPoint(pt4Lat, pt4Lon)
	real := pt1.Distance(pt2)
	fake := pt1.Approximately(pt2)
	off := (math.Abs(fake-real) / real) * 100.0
	t.Logf("REAL:%f FAKE:%f (%.1f)", real, fake, off)

	real = pt3.Distance(pt4)
	fake = pt3.Approximately(pt4)
	off = (math.Abs(fake-real) / real) * 100.0
	t.Logf("REAL:%f FAKE:%f (%.1f)", real, fake, off)

}

func TestAccuracyPart3(t *testing.T) {
	tests := []struct {
		pt1lat, pt1lon, pt2lat, pt2lon, distance float64
	}{
		{SFLat, SFLon, ZepLat, ZepLon, SFtoZep},
		{SFLat, SFLon, AlaLat, AlaLon, SFtoAla},
		{SFLat, SFLon, HouLat, HouLon, SFtoHou},
		{AlaLat, AlaLon, PortLat, PortLon, AlaToPDX},
	}
	for _, tt := range tests {
		real := Distance(tt.pt1lat, tt.pt1lon, tt.pt2lat, tt.pt2lon)
		fake := ApproximateDistance(tt.pt1lat, tt.pt1lon, tt.pt2lat, tt.pt2lon)
		off1 := (math.Abs(tt.distance-real) / tt.distance) * 100.0
		off2 := (math.Abs(tt.distance-fake) / tt.distance) * 100.0
		t.Logf("from (% 9.6f, % 9.6f)->(% 9.6f, % 11.6f) real:% 7.1f calc: %7.1f (%.2f%%) fake: %7.1f (%.2f%%)",
			tt.pt1lat, tt.pt1lon, tt.pt2lat, tt.pt2lon, tt.distance, real, off1, fake, off2,
		)

	}
}

func TestLonLookup(t *testing.T) {
	lat := 37.73
	km1 := LonKilos(lat) //LongitudeKilometers(lat, 1.0)
	km2 := LookupLonKmPerLat(lat)
	off := math.Abs((km1-km2)/km1) * 100.0
	assert.Less(t, off, 1.0)
	//t.Logf("calc: %f vs. look: %f (%.2f)", km1, km2, off)
}

func TestLonAccuracy(t *testing.T) {
	for i := 0.9; i < 90.0; i += 1.0 {
		km1 := LonKilos(i)
		km2 := LookupLonKmPerLat(i)
		//km2 := LookupLonKmPerLatInt(int(lat))
		off := math.Abs(((km1 - km2) / km1)) * 100.0
		assert.Less(t, off, 1.0)
		//assert.Less(t, 1.0, off)
		//fmt.Printf("calc %4.1f: % 12.7f vs. look: % 12.7f (%.2f%%)\n", i, km1, km2, off)
	}
}

func TestLonAccuracyRedux(t *testing.T) {
	for i := 0.05; i < 50.0; i += 0.1 {
		km1 := LonKilos(i)
		km2 := LookupLonKmPerLat(i)
		//km2 := LookupLonKmPerLatInt(int(lat))
		off := math.Abs(((km1 - km2) / km1)) * 100.0
		assert.Less(t, off, 1.0)
		//fmt.Printf("calc %4.1f: % 12.7f vs. look: % 12.7f (%.2f%%)\n", i, km1, km2, off)
	}
}

func BenchmarkLonDistanceCalc(b *testing.B) {
	const lat, lon = 37.73, -122.34
	for i := 0; i < b.N; i++ {
		_ = LongitudeKilometers(lat, lon)
	}
}

/*
func BenchmarkCalculateLongitudeKilometers(b *testing.B) {
	const lat, lon = 37.73, -122.34
	for i := 0; i < b.N; i++ {
		_ = LongitudeKilometers(lat, lon)
	}
}
*/

func BenchmarkLookupLonKmPerLat(b *testing.B) {
	const lat = 37.73
	for i := 0; i < b.N; i++ {
		_ = LookupLonKmPerLat(lat)
	}
}

func BenchmarkDistanceReal(b *testing.B) {
	const (
		pt1Lat, pt1Lon = 37.7690679, -122.2592744
		pt2Lat, pt2Lon = 37.78978486773513, -122.39617851945805
	)
	for i := 0; i < b.N; i++ {
		_ = Distance(pt1Lat, pt1Lon, pt2Lat, pt2Lon)
	}
}

/*
func BenchmarkApproximateDistance(b *testing.B) {
	const (
		pt1Lat, pt1Lon = 37.7690679, -122.2592744
		pt2Lat, pt2Lon = 37.78978486773513, -122.39617851945805
	)
	for i := 0; i < b.N; i++ {
		_ = ApproximateDistanceGeo(pt1Lat, pt1Lon, pt2Lat, pt2Lon)
	}
}
*/

type HeatData struct {
	Lat, Lon float64
	Index    int
}

type Heated []HeatData

func (h Heated) Len() int {
	return len(h)
}

func (h Heated) IndexPoint(i int) Point {
	v := h[i]
	return GeoPoint(v.Lat, v.Lon)
}

func (h *HeatData) Import(ss []string) error {
	lat, err := strconv.ParseFloat(ss[0], 64)
	if err != nil {
		return fmt.Errorf("bad lat: %w", err)
	}
	lon, err := strconv.ParseFloat(ss[1], 64)
	if err != nil {
		return fmt.Errorf("bad lat: %w", err)
	}
	idx, err := strconv.Atoi(ss[2])
	if err != nil {
		return fmt.Errorf("bad idx: %w", err)
	}
	h.Lat = lat
	h.Lon = lon
	h.Index = idx
	return nil
}

const heatFile = "testdata/heat.csv.gz"

func TestHeated(t *testing.T) {
	heated, err := loadHeated(heatFile)
	if err != nil {
		t.Fatal(err)
	}
	for i, h := range heated {
		t.Logf("%v", h)
		if i > 10 {
			break
		}
	}
}

type Helper interface {
	Helper()
	Fatal(...interface{})
}

func testHeat(t Helper) Heated {
	heated, err := loadHeated(heatFile)
	if err != nil {
		t.Fatal(err)
	}
	return heated
}

func TestClosestHeat(t *testing.T) {
	heated := testHeat(t)
	pt := GeoPoint(AlaLat, AlaLon)
	const within = 10.0
	idx, dist := Closest(heated, pt, within)
	h := heated[idx]
	t.Logf("%d/%d:(%f) %v", idx, len(heated), dist, h)
}

func TestBestestHeat(t *testing.T) {
	heated := testHeat(t)
	pt := GeoPoint(AlaLat, AlaLon)
	const within = 10.0
	idx, dist := Bestest(heated, pt, within)
	h := heated[idx]
	t.Logf("%d/%d:(%f) %v", idx, len(heated), dist, h)
}

func TestClosestAllocs(t *testing.T) {
	heated := testHeat(t)
	pt := GeoPoint(AlaLat, AlaLon)
	const within = 10.0
	allocs := testing.AllocsPerRun(100, func() {
		Closest(heated, pt, within)
	})
	t.Logf("%.0f", allocs)
}

func TestNearest(t *testing.T) {
	pt := GeoPoint(AlaLat, AlaLon)
	info, err := Nearest(heatFile, pt, true)
	if err != nil {
		t.Fatal(err)
	}
	t.Log(info)
}
func loadHeated(filename string) (Heated, error) {
	count := 0
	var heated Heated
	fn := func(ss []string) error {
		count++
		if count == 1 {
			return nil
		}

		var h HeatData
		if err := (&h).Import(ss); err != nil {
			return err
		}
		heated = append(heated, h)
		return nil
	}
	return heated, loadCSV(filename, fn)
}

func loadCSV(filename string, fn func([]string) error) error {
	f, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	// TODO: experiment with buffer sizes
	b := bufio.NewReader(f)
	r := io.Reader(b)

	gzipped := filepath.Ext(filename) == ".gz"
	if gzipped {
		gzr, err := gzip.NewReader(b)
		if err != nil {
			return fmt.Errorf("gzip fail: %w", err)
		}
		defer gzr.Close()
		r = gzr
	}
	cr := csv.NewReader(r)
	for {
		records, err := cr.Read()
		if len(records) > 0 {
			if err := fn(records); err != nil {
				return err
			}
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
	}
	return nil
}

func TestRightmost(t *testing.T) {
	list := []int{1, 2, 2, 2, 2, 2, 3, 4, 7, 7, 7, 7, 7, 9, 9, 9, 9, 11, 15, 23, 42, 69}
	val := 9
	T := sort.SearchInts(list, val)
	t.Logf("original value is %d at idx %d", val, T)

	T--
	prior := list[T]
	t.Logf("prior value last is %d at %d", prior, T)
	// T := 7
	idx := Rightmost(T, func(i int) bool {
		fnd := list[i]
		return fnd > T
	})
	// t.Logf("for target idx %d (%d) found idx:%d with value %d", T, list[T], idx, list[idx])
	t.Logf("found idx:%d with value %d", idx, list[idx])

}

func TestPrior(t *testing.T) {
	list := []int{1, 2, 2, 2, 2, 2, 3, 4, 7, 7, 7, 7, 7, 9, 9, 9, 9, 11, 15, 23, 42, 69}
	val := 9
	T := sort.SearchInts(list, val)
	t.Logf("original value is %d at idx %d", val, T)

	T--
	prior := list[T]
	t.Logf("prior value last is %d at %d", prior, T)

	idx := sort.SearchInts(list[:T], prior)
	// t.Logf("original value is %d at idx %d", val, T)

	t.Logf("found idx:%d with value %d", idx, list[idx])

}
