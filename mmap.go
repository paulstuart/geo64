package geo

import (
	"errors"
	"io"
	"sort"

	"github.com/tidwall/mmap"
)

var ErrNotFound = errors.New("not found")

type Decoder interface {
	Decode([]byte) error
	Size() int // size of struct
	Point() Point
	Less(Point) bool
	JSON(w io.Writer) error
}

type MFile struct {
	B []byte
}

type Iter struct {
	m *MFile
	d Decoder
}

func (m *MFile) Close() error {
	return mmap.Close(m.B)
}

func (m *Iter) Len() int {
	return len(m.m.B) / m.d.Size()
}

func (m *Iter) IndexPoint(i int) Point {
	off := m.d.Size() * i
	end := off + m.d.Size()
	if err := m.d.Decode(m.m.B[off:end]); err != nil {
		panic(err)
	}
	return m.d.Point()
}

func (m *Iter) Load(i int) {
	off := m.d.Size() * i
	end := off + m.d.Size()
	if err := m.d.Decode(m.m.B[off:end]); err != nil {
		panic(err)
	}
}

func (m *Iter) Less(pt Point) bool {
	return m.d.Less(pt)
}

func (m *Iter) JSON(w io.Writer) {
	m.d.JSON(w)
}

func Mmap(filename string) (*MFile, error) {
	b, err := mmap.Open(filename, false)
	if err != nil {
		return nil, err
	}
	return &MFile{b}, err
}

func (m *MFile) ReadAt(p []byte, i int64) (int, error) {
	if i > int64(len(m.B)) {
		return 0, errors.New("index exceeds file size")
	}
	return copy(p, m.B[i:]), nil
}

func (m *MFile) NewIter(d Decoder) *Iter {
	return &Iter{
		m: m,
		d: d,
	}
}

func (m *Iter) Get(i int) interface{} {
	off := m.d.Size() * i
	end := off + m.d.Size()
	if err := m.d.Decode(m.m.B[off:end]); err != nil {
		panic(err)
	}
	return m.d
}

type Container interface {
	ContainsPoint(Point) bool
}

func (m *Iter) Ranger(from, to Point, fn func(interface{}), ctr Container) error {
	size := m.Len()
	idx := sort.Search(size, func(i int) bool {
		return from.Less(m.IndexPoint(i))
	})
	if idx == size {
		return ErrNotFound
	}
	for {
		m.Load(idx)
		if !m.Less(to) {
			break
		}
		pt := m.d.Point()
		if between(pt.Lon, from.Lon, to.Lon) {
			if ctr == nil || ctr.ContainsPoint(m.d.Point()) {
				fn(m.d)
			}
		}
		idx++
	}
	return nil
}
