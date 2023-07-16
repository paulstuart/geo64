package geo

import (
	"bufio"
	"compress/gzip"
	"io"
	"math"
	"os"
	"path/filepath"
)

type LineInfo struct {
	Index    int
	Line     string
	Distance float64
}

/*
// LinePoint extracts coordinates from the beginning
// of a line of csv text
func LinePoint(s string, latLon bool) (Point, error) {
	c1 := strings.Index(s, ",")
	if c1 < 1 {
		return Point{}, fmt.Errorf("no commas")
	}
	ss := s[c1+1:]
	c2 := strings.Index(ss, ",")
	if c2 < 1 {
		return Point{}, fmt.Errorf("no terminator")
	}
	s1 := s[:c1]
	s2 := ss[:c2]
	f1, err := strconv.ParseFloat(s1, 32)
	if err != nil {
		return Point{}, fmt.Errorf("bad coord: %w", err)
	}
	f2, err := strconv.ParseFloat(s2, 32)
	if err != nil {
		return Point{}, fmt.Errorf("moar bad: %w", err)
	}
	if latLon {
		return Point{GeoType(f1), GeoType(f2)}, nil
	}
	return Point{GeoType(f2), GeoType(f1)}, nil
}
*/

// Nearest scans a csv file with lon,lat coordinates
// and returns the line that is closest to the given point
func Nearest(filename string, pt Point, latLon bool) (LineInfo, error) {
	info := LineInfo{
		Distance: math.MaxFloat64,
	}
	var idx int
	fn := func(s string) error {
		idx++
		there, err := QueryPoint(s)
		if err != nil {
			return nil // should we log it?
		}
		if !latLon {
			pt.Lat, pt.Lon = pt.Lon, pt.Lat
		}
		dist := pt.Distance(there)
		if dist < info.Distance {
			info.Distance = dist
			info.Line = s
			info.Index = idx
		}
		return nil
	}
	return info, LoadLines(filename, fn)
}

func LoadLines(filename string, fn func(string) error) error {
	f, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	// TODO: experiment with buffer sizes
	b := bufio.NewReader(f)
	r := io.Reader(b)

	// TODO: add lzma and others
	gzipped := filepath.Ext(filename) == ".gz"
	if gzipped {
		gzr, err := gzip.NewReader(b)
		if err != nil {
			return err
		}
		defer gzr.Close()
		r = gzr
	}

	scan := bufio.NewScanner(r)
	for scan.Scan() {
		line := scan.Text()
		if err = fn(line); err != nil {
			return err
		}
	}
	return nil
}
