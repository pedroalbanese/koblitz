// Parameters for the Koblitz (secp256k1) Elliptic curve
package secp256k1

import (
	"crypto/ecdsa"
	"crypto/elliptic"
	"encoding/asn1"
	"errors"
	"math/big"
	"sync"
)

var (
	oidS256 = asn1.ObjectIdentifier{1, 3, 132, 0, 10}
)

var initonce sync.Once
var s256 *Curve

func initS256() {
	s256 = new(Curve)
	s256.P, _ = new(big.Int).SetString("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16)
	s256.N, _ = new(big.Int).SetString("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
	s256.B, _ = new(big.Int).SetString("0000000000000000000000000000000000000000000000000000000000000007", 16)
	s256.Gx, _ = new(big.Int).SetString("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16)
	s256.Gy, _ = new(big.Int).SetString("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16)
	s256.BitSize = 256
	s256.Name = "secp256k1"
}

func S256() elliptic.Curve {
	initonce.Do(initS256)
	return s256
}

type Curve struct {
	P       *big.Int // the order of the underlying field
	N       *big.Int // the order of the base point
	B       *big.Int // the constant of the Curve equation
	Gx, Gy  *big.Int // (x,y) of the base point
	BitSize int      // the size of the underlying field
	Name    string   // name of the curve
}

func (curve *Curve) Params() *elliptic.CurveParams {
	return &elliptic.CurveParams{
		P:       curve.P,
		N:       curve.N,
		B:       curve.B,
		Gx:      curve.Gx,
		Gy:      curve.Gy,
		BitSize: curve.BitSize,
	}
}

// IsOnCurve checks if the point (x, y) is on the curve defined by y^2 = x^3 + B mod P.
func (curve *Curve) IsOnCurve(x, y *big.Int) bool {
	y2 := new(big.Int).Mul(y, y) 
	y2.Mod(y2, curve.P)          

	x3 := new(big.Int).Mul(x, x) 
	x3.Mul(x3, x)                

	x3.Add(x3, curve.B)
	x3.Mod(x3, curve.P)

	return x3.Cmp(y2) == 0
}

// affineFromJacobian converts Jacobian coordinates (x, y, z) to affine coordinates (xOut, yOut).
func (curve *Curve) affineFromJacobian(x, y, z *big.Int) (xOut, yOut *big.Int) {
	zinv := new(big.Int).ModInverse(z, curve.P)
	zinvsq := new(big.Int).Mul(zinv, zinv)

	xOut = new(big.Int).Mul(x, zinvsq)
	xOut.Mod(xOut, curve.P)
	zinvsq.Mul(zinvsq, zinv)
	yOut = new(big.Int).Mul(y, zinvsq)
	yOut.Mod(yOut, curve.P)
	return
}

// Add adds two points (x1, y1) and (x2, y2) on the curve and returns the result in affine coordinates.
func (curve *Curve) Add(x1, y1, x2, y2 *big.Int) (*big.Int, *big.Int) {
	z := new(big.Int).SetInt64(1)
	return curve.affineFromJacobian(curve.addJacobian(x1, y1, z, x2, y2, z))
}

// addJacobian performs point addition in Jacobian coordinates and returns the result.
func (curve *Curve) addJacobian(x1, y1, z1, x2, y2, z2 *big.Int) (*big.Int, *big.Int, *big.Int) {
	z1z1 := new(big.Int).Mul(z1, z1)
	z1z1.Mod(z1z1, curve.P)
	z2z2 := new(big.Int).Mul(z2, z2)
	z2z2.Mod(z2z2, curve.P)

	u1 := new(big.Int).Mul(x1, z2z2)
	u1.Mod(u1, curve.P)
	u2 := new(big.Int).Mul(x2, z1z1)
	u2.Mod(u2, curve.P)
	h := new(big.Int).Sub(u2, u1)
	if h.Sign() == -1 {
		h.Add(h, curve.P)
	}
	i := new(big.Int).Lsh(h, 1)
	i.Mul(i, i)
	j := new(big.Int).Mul(h, i)

	s1 := new(big.Int).Mul(y1, z2)
	s1.Mul(s1, z2z2)
	s1.Mod(s1, curve.P)
	s2 := new(big.Int).Mul(y2, z1)
	s2.Mul(s2, z1z1)
	s2.Mod(s2, curve.P)
	r := new(big.Int).Sub(s2, s1)
	if r.Sign() == -1 {
		r.Add(r, curve.P)
	}
	r.Lsh(r, 1)
	v := new(big.Int).Mul(u1, i)

	x3 := new(big.Int).Set(r)
	x3.Mul(x3, x3)
	x3.Sub(x3, j)
	x3.Sub(x3, v)
	x3.Sub(x3, v)
	x3.Mod(x3, curve.P)

	y3 := new(big.Int).Set(r)
	v.Sub(v, x3)
	y3.Mul(y3, v)
	s1.Mul(s1, j)
	s1.Lsh(s1, 1)
	y3.Sub(y3, s1)
	y3.Mod(y3, curve.P)

	z3 := new(big.Int).Add(z1, z2)
	z3.Mul(z3, z3)
	z3.Sub(z3, z1z1)
	if z3.Sign() == -1 {
		z3.Add(z3, curve.P)
	}
	z3.Sub(z3, z2z2)
	if z3.Sign() == -1 {
		z3.Add(z3, curve.P)
	}
	z3.Mul(z3, h)
	z3.Mod(z3, curve.P)

	return x3, y3, z3
}

// Double doubles the point (x1, y1) on the curve and returns the result in affine coordinates.
func (curve *Curve) Double(x1, y1 *big.Int) (*big.Int, *big.Int) {
	z1 := new(big.Int).SetInt64(1)
	return curve.affineFromJacobian(curve.doubleJacobian(x1, y1, z1))
}

// doubleJacobian performs point doubling in Jacobian coordinates and returns the result.
func (curve *Curve) doubleJacobian(x, y, z *big.Int) (*big.Int, *big.Int, *big.Int) {
	a := new(big.Int).Mul(x, x) 
	b := new(big.Int).Mul(y, y) 
	c := new(big.Int).Mul(b, b) 

	d := new(big.Int).Add(x, b) 
	d.Mul(d, d)                 
	d.Sub(d, a)                 
	d.Sub(d, c)                 
	d.Mul(d, big.NewInt(2))     

	e := new(big.Int).Mul(big.NewInt(3), a) 
	f := new(big.Int).Mul(e, e)             

	x3 := new(big.Int).Mul(big.NewInt(2), d) 
	x3.Sub(f, x3)                            
	x3.Mod(x3, curve.P)

	y3 := new(big.Int).Sub(d, x3)                  
	y3.Mul(e, y3)                                  
	y3.Sub(y3, new(big.Int).Mul(big.NewInt(8), c)) 
	y3.Mod(y3, curve.P)

	z3 := new(big.Int).Mul(y, z) 
	z3.Mul(big.NewInt(2), z3)    
	z3.Mod(z3, curve.P)

	return x3, y3, z3
}

// ScalarMult computes the scalar multiplication of a point (Bx, By) by a scalar k.
func (curve *Curve) ScalarMult(Bx, By *big.Int, k []byte) (*big.Int, *big.Int) {
	Bz := new(big.Int).SetInt64(1)
	x := Bx
	y := By
	z := Bz

	seenFirstTrue := false
	for _, byte := range k {
		for bitNum := 0; bitNum < 8; bitNum++ {
			if seenFirstTrue {
				x, y, z = curve.doubleJacobian(x, y, z)
			}
			if byte&0x80 == 0x80 {
				if !seenFirstTrue {
					seenFirstTrue = true
				} else {
					x, y, z = curve.addJacobian(Bx, By, Bz, x, y, z)
				}
			}
			byte <<= 1
		}
	}

	if !seenFirstTrue {
		return nil, nil
	}

	return curve.affineFromJacobian(x, y, z)
}

// ScalarBaseMult computes the scalar multiplication of the curve's base by a scalar k.
func (curve *Curve) ScalarBaseMult(k []byte) (*big.Int, *big.Int) {
	return curve.ScalarMult(curve.Gx, curve.Gy, k)
}

// CompressPoint compresses the point (X, Y) into a compact representation.
func CompressPoint(curve *Curve, X, Y *big.Int) (cp []byte) {
    return curve.CompressPoint(X, Y)
}

// CompressPoint compresses the point (X, Y) into a compact representation.
func (curve *Curve) CompressPoint(X, Y *big.Int) (cp []byte) {
	by := new(big.Int).And(Y, big.NewInt(1)).Int64()
	bx := X.Bytes()
	cp = make([]byte, len(bx)+1)
	if by == 1 {
		cp[0] = byte(3)
	} else {
		cp[0] = byte(2)
	}
	copy(cp[1:], bx)

	return
}

// DecompressPoint decompresses a point from its compact representation.
func (curve *Curve) DecompressPoint(cp []byte) (X, Y *big.Int, err error) {
	var c int64

	switch cp[0] {
	case byte(0x03):
		c = 1
		break
	case byte(0x02):
		c = 0
		break
	case byte(0x04):
		X, Y = elliptic.Unmarshal(curve, cp)
		return
	default:
		return nil, nil, errors.New("Not a compressed point. (Invalid Header)")
	}

	byteLen := (curve.Params().BitSize + 7) >> 3
	if len(cp) != 1+byteLen {
		return nil, nil, errors.New("Not a compressed point. (Require 1 + key size)")
	}

	X = new(big.Int).SetBytes(cp[1:])
	Y = new(big.Int)

	Y.Mod(Y.Mul(X, X), curve.P)
	Y.Mod(Y.Mul(Y, X), curve.P)
	Y.Mod(Y.Add(Y, curve.B), curve.P)

	Y = curve.Sqrt(Y)

	if Y.Cmp(big.NewInt(0)) == 0 {
		return nil, nil, errors.New("Not a compressed point. (Not on curve)")
	}

	if c != new(big.Int).And(Y, big.NewInt(1)).Int64() {
		Y.Sub(curve.P, Y)
	}

	return
}

// Sqrt calculates the modular square root of a, returning the root if it exists.
func (curve *Curve) Sqrt(a *big.Int) *big.Int {
	ZERO := big.NewInt(0)
	ONE := big.NewInt(1)
	TWO := big.NewInt(2)
	THREE := big.NewInt(3)
	FOUR := big.NewInt(4)

	p := curve.P
	c := new(big.Int)
	
	if a.Cmp(ZERO) == 0 {
		return ZERO
	} else if p.Cmp(TWO) == 0 {
		return a.Mod(a,p)
	} else if LegendreSymbol(a, p) != 1 {
		return ZERO
	} else if c.Mod(p, FOUR).Cmp(THREE) == 0 {
		c.Add(p, ONE)
		c.Div(c, FOUR)
		c.Exp(a, c, p)
		return c
	}

	s := new(big.Int)
	s.Sub(p, ONE)

	e := new(big.Int)
	e.Set(ZERO)
	for c.Mod(s, TWO).Cmp(ZERO) == 0 {
		s.Div(s, TWO)
		e.Add(e, ONE)
	}
	
	n := new(big.Int)
	n.Set(TWO)
	for LegendreSymbol(n, p) != -1 {
		n.Add(n, ONE)
	}
	
	x := new(big.Int)
	x.Add(s, ONE)
	x.Div(x, TWO)
	x.Exp(a, x, p)
	
	b := new(big.Int)
	b.Exp(a, s, p)

	g := new(big.Int)
	g.Exp(n, s, p)

	r := new(big.Int)
	r.Set(e)

	t := new(big.Int)
	m := new(big.Int)
	gs := new(big.Int)

	for {
		t.Set(b)
		m.Set(ZERO)

		for ; m.Cmp(r) < 0; m.Add(m, ONE) {
			if t.Cmp(ONE) == 0 {
				break
			}
			t.Exp(t, TWO, p)
		}

		if m.Cmp(ZERO) == 0 {
			return x
		}

		gs.Sub(r, m)
		gs.Sub(gs, ONE)
		gs.Exp(TWO, gs, nil)
		gs.Exp(g, gs, p)

		g.Mod(g.Mul(gs, gs), p)
		x.Mod(x.Mul(x, gs), p)
		b.Mod(b.Mul(b, g), p)
		r.Set(m)
	}
}

// LegendreSymbol calculates the Legendre symbol for a and p.
func LegendreSymbol(a, p *big.Int) int {
	ZERO := big.NewInt(0)
	ONE := big.NewInt(1)
	TWO := big.NewInt(2)

	ls := new(big.Int).Mod(a, p)

	if ls.Cmp(ZERO) == 0 {
		return 0 
	}

	ps := new(big.Int).Sub(p, ONE)

	ls.Div(ps, TWO)
	ls.Exp(a, ls, p)

	if c := ls.Cmp(ps); c == 0 {
		return -1 
	}

	return 1 
}

// Define pkAlgorithmIdentifier to avoid undefined identifier
type pkAlgorithmIdentifier struct {
	Algorithm  asn1.ObjectIdentifier
	Parameters asn1.RawValue
}

type PublicKey struct {
	X, Y  *big.Int
	Curve elliptic.Curve
}

type PrivateKey struct {
	PublicKey PublicKey
	D         *big.Int
}

func (pk *PublicKey) MarshalPKCS8PublicKey(curve elliptic.Curve) ([]byte, error) {
	// Marshal the public key coordinates
	derBytes := elliptic.Marshal(curve, pk.X, pk.Y)

	// Determine the OID based on the curve
	var oid asn1.ObjectIdentifier
	switch curve {
	case S256():
		oid = oidS256
	default:
		return nil, errors.New("unsupported curve")
	}

	// Create a SubjectPublicKeyInfo structure
	subjectPublicKeyInfo := struct {
		Algorithm pkAlgorithmIdentifier
		PublicKey asn1.BitString
	}{
		Algorithm: pkAlgorithmIdentifier{
			Algorithm:  oid,
			Parameters: asn1.RawValue{Tag: asn1.TagOID, Bytes: []byte(oid.String())},
		},
		PublicKey: asn1.BitString{Bytes: derBytes, BitLength: len(derBytes) * 8},
	}

	// Marshal the SubjectPublicKeyInfo structure
	derBytes, err := asn1.Marshal(subjectPublicKeyInfo)
	if err != nil {
		return nil, err
	}

	return derBytes, nil
}

func ParsePublicKey(der []byte) (*PublicKey, error) {
	var publicKeyInfo struct {
		Algorithm pkAlgorithmIdentifier
		PublicKey asn1.BitString
	}

	_, err := asn1.Unmarshal(der, &publicKeyInfo)
	if err != nil {
		return nil, err
	}

	var curve elliptic.Curve
	switch {
	case publicKeyInfo.Algorithm.Algorithm.Equal(oidS256):
		curve = S256()
	default:
		return nil, errors.New("unsupported curve OID")
	}

	// Check if the public key bytes are empty
	if len(publicKeyInfo.PublicKey.Bytes) == 0 {
		return nil, errors.New("public key bytes are empty")
	}

	// Unmarshal the public key coordinates
	X, Y := elliptic.Unmarshal(curve, publicKeyInfo.PublicKey.Bytes)
	if X == nil || Y == nil {
		return nil, errors.New("failed to unmarshal public key")
	}

	// Return the parsed public key with the determined curve
	return &PublicKey{X: X, Y: Y, Curve: curve}, nil
}

func (pk *PrivateKey) MarshalPKCS8PrivateKey(curve elliptic.Curve) ([]byte, error) {
	if !curve.IsOnCurve(pk.PublicKey.X, pk.PublicKey.Y) {
		return nil, errors.New("Public key is not on the curve")
	}

	// Convert the private key D to bytes
	dBytes := pk.D.Bytes()

	curveSize := (curve.Params().BitSize + 7) / 8
	if len(dBytes) < curveSize {
		padding := make([]byte, curveSize-len(dBytes))
		dBytes = append(padding, dBytes...)
	}

	// Determine the OID based on the curve
	var oid asn1.ObjectIdentifier
	switch curve {
	case S256():
		oid = oidS256
	default:
		return nil, errors.New("unsupported curve")
	}

	// Create a PrivateKeyInfo structure
	privateKeyInfo := struct {
		Version             int
		PrivateKeyAlgorithm pkAlgorithmIdentifier
		PublicKey           struct {
			X *big.Int
			Y *big.Int
		}
		PrivateKey []byte
	}{
		Version: 0,
		PrivateKeyAlgorithm: pkAlgorithmIdentifier{
			Algorithm:  oid,
			Parameters: asn1.RawValue{Tag: asn1.TagOID, Bytes: []byte(oid.String())},
		},
		PublicKey: struct {
			X *big.Int
			Y *big.Int
		}{
			X: new(big.Int).SetBytes(pk.PublicKey.X.Bytes()),
			Y: new(big.Int).SetBytes(pk.PublicKey.Y.Bytes()),
		},
		PrivateKey: dBytes,
	}

	// Marshal the PrivateKeyInfo structure
	derBytes, err := asn1.Marshal(privateKeyInfo)
	if err != nil {
		return nil, err
	}

	return derBytes, nil
}

func ParsePrivateKey(der []byte) (*PrivateKey, error) {
	var privateKeyInfo struct {
		Version             int
		PrivateKeyAlgorithm pkAlgorithmIdentifier
		PublicKey           struct {
			X *big.Int
			Y *big.Int
		}
		PrivateKey []byte
	}
	_, err := asn1.Unmarshal(der, &privateKeyInfo)
	if err != nil {
		return nil, err
	}

	// Determine the curve based on the OID
	var curve elliptic.Curve
	switch {
	case privateKeyInfo.PrivateKeyAlgorithm.Algorithm.Equal(oidS256):
		curve = S256()
	default:
		return nil, errors.New("unsupported curve OID")
	}

	X := privateKeyInfo.PublicKey.X
	Y := privateKeyInfo.PublicKey.Y
	D := new(big.Int).SetBytes(privateKeyInfo.PrivateKey)

	if !curve.IsOnCurve(X, Y) {
		return nil, errors.New("Public key is not on the curve")
	}

	// Create and return the private key with the determined curve
	privateKey := &PrivateKey{
		PublicKey: PublicKey{
			X:     X,
			Y:     Y,
			Curve: curve,
		},
		D: D,
	}

	return privateKey, nil
}

func (pk *PublicKey) ToECDSA() *ecdsa.PublicKey {
	return &ecdsa.PublicKey{
		Curve: pk.Curve,
		X:     pk.X,
		Y:     pk.Y,
	}
}

func (pk *PrivateKey) ToECDSAPrivateKey() *ecdsa.PrivateKey {
	return &ecdsa.PrivateKey{
		PublicKey: ecdsa.PublicKey{
			Curve: pk.PublicKey.Curve,
			X:     pk.PublicKey.X,
			Y:     pk.PublicKey.Y,
		},
		D: pk.D,
	}
}

func NewPrivateKey(privateKey *ecdsa.PrivateKey) *PrivateKey {
	return &PrivateKey{
		PublicKey: PublicKey{
			Curve: privateKey.PublicKey.Curve,
			X:     privateKey.PublicKey.X,
			Y:     privateKey.PublicKey.Y,
		},
		D: privateKey.D,
	}
}

func ECDH(privateKey *ecdsa.PrivateKey, publicKey *ecdsa.PublicKey) ([]byte, error) {
	// Compute shared key
	x, _ := privateKey.Curve.ScalarMult(publicKey.X, publicKey.Y, privateKey.D.Bytes())
	return x.Bytes(), nil
}
