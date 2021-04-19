import hashlib
import sys
import binascii
sys.setrecursionlimit(1000000)

b = 256
q = 2**255 - 19
l = 2**252 + 27742317777372353535851937790883648493

def H(m):
  return hashlib.sha512(m).digest()

def expmod(b,e,m):
  if e == 0: return 1
  t = expmod(b,e//2,m)**2 % m #fix //
  if e & 1: t = (t*b) % m
  return t

def inv(x):
  return expmod(x,q-2,q)

d = -121665 * inv(121666)
I = expmod(2,(q-1)//4,q) #fix //

print("d:        " + str(d % q))
print("I:        " + str(I))

print("inv of d: " + str(inv(d) % q))

def xrecover(y):
  xx = (y*y-1) * inv(d*y*y+1)
  x = expmod(xx,(q+3)//8,q) #fix //
  if (x*x - xx) % q != 0: x = (x*I) % q
  if x % 2 != 0: x = q-x
  return x

By = 4 * inv(5)
Bx = xrecover(By)
B = [Bx % q,By % q]
print("B[0]:     " + hex(B[0]))
print("B[1]:     " + hex(B[1]))

def edwards(P,Q):
  x1 = P[0]
  y1 = P[1]
  x2 = Q[0]
  y2 = Q[1]
  x3 = (x1*y2+x2*y1) * inv(1+d*x1*x2*y1*y2)
  y3 = (y1*y2+x1*x2) * inv(1-d*x1*x2*y1*y2)
  return [x3 % q,y3 % q]

def scalarmult(P,e):
  if e == 0: return [0,1]
  Q = scalarmult(P,e//2) #fix
  Q = edwards(Q,Q)
  if e & 1: Q = edwards(Q,P)
  return Q

def encodeint(y):
  bits = [(y >> i) & 1 for i in range(b)]
  return b''.join([sum([bits[i * 8 + j] << j for j in range(8)]).to_bytes(1, 'little') for i in range(b//8)]) #fix

def encodepoint(P):
  x = P[0] % q
  y = P[1] % q
#  print("x (" + str(type(x)) + ") = " + str(x))
#  print("y (" + str(type(y)) + ") = " + str(y))
  bits = [(y >> i) & 1 for i in range(b - 1)] + [x & 1]
#  for i in range(b):
#    print("i = " + str(i) + " - " + str(bits[i]))
#  print()

#  print("buf = " + binascii.hexlify(b''.join([sum([bits[i * 8 + j] << j for j in range(8)]).to_bytes(1, 'little') for i in range(b//8)])).decode("utf-8"))
  return b''.join([sum([bits[i * 8 + j] << j for j in range(8)]).to_bytes(1, 'little') for i in range(b//8)]) #fix

def bit(h,i):
  return (h[i//8] >> (i%8)) & 1 #fix

def publickey(sk):
#  print("PUBLIC KEY:")
  h = H(sk)
#  print("hash(sk) = " + binascii.hexlify(h).decode("utf-8"))
#  for i in range(3, b-2):
#    print(bit(h, i), end="")
#  print()
  a = 2**(b-2) + sum(2**i * bit(h,i) for i in range(3,b-2))
#  print("a = " + hex(a))
  A = scalarmult(B,a)
#  print("A[0] = " + hex(A[0]))
#  print("A[1] = " + hex(A[1]))
#  print("Public key = " + binascii.hexlify(encodepoint(A)).decode("utf-8"))
  return encodepoint(A)

def Hint(m):
  h = H(m)
#  print("hash(bits) = " + binascii.hexlify(h).decode("utf-8"))
#  for i in range(2*b):
#    print(bit(h, i), end="")
#  print()
#  print(hex(sum(2**i * bit(h,i) for i in range(2*b)) % q))
  return sum(2**i * bit(h,i) for i in range(2*b))

def signature(m,sk,pk):
  print("SIGNATURE")
#  print("m = " + binascii.hexlify(m).decode("utf-8"))
#  print("sk = " + binascii.hexlify(sk).decode("utf-8"))
  h = H(sk)
#  print("hash(sk) = " + binascii.hexlify(h).decode("utf-8"))
  a = 2**(b-2) + sum(2**i * bit(h,i) for i in range(3,b-2))
  print("a = " + hex(a))
  print("bits = " + binascii.hexlify(b''.join([h[i].to_bytes(1, 'little') for i in range(b//8,b//4)]) + m).decode("utf-8"))
  r = Hint(b''.join([h[i].to_bytes(1, 'little') for i in range(b//8,b//4)]) + m) #fix
  print("r = " + hex(r))
  print("r % q = " + hex(r % q))
  R = scalarmult(B,r)
  print("R = " + binascii.hexlify(encodepoint(R)).decode("utf-8"))
  S = (r + Hint(encodepoint(R) + pk + m) * a) % l
  print("S = " + hex(S))
  print("Signature = " + binascii.hexlify(encodepoint(R) + encodeint(S)).decode("utf-8"))
  return encodepoint(R) + encodeint(S)

def isoncurve(P):
  x = P[0]
  y = P[1]
  return (-x*x + y*y - 1 - d*x*x*y*y) % q == 0

def decodeint(s):
  return sum(2**i * bit(s,i) for i in range(0,b))

def decodepoint(s):
  y = sum(2**i * bit(s,i) for i in range(0,b-1))
  print(hex(y % q));
  x = xrecover(y)
  if x & 1 != bit(s,b-1): x = q-x
  P = [x,y]
  if not isoncurve(P): raise Exception("decoding point that is not on curve")
  return P

def checkvalid(s,m,pk):
  if len(s) != b/4: raise Exception("signature length is wrong")
  if len(pk) != b/8: raise Exception("public-key length is wrong")
  R = decodepoint(s[0:b//8])
  A = decodepoint(pk)
  S = decodeint(s[b//8:b//4])
  h = Hint(encodepoint(R) + pk + m)
  if scalarmult(B,S) != edwards(R,scalarmult(A,h)):
    raise Exception("signature does not pass verification")
