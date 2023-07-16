//go:build debug

package geo

import (
	"log"
)

func debugf(text string, args ...interface{}) {
	log.Printf("DBG: "+text, args...)
}
