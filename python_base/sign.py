import sys
import binascii
import ed25519
import hashlib


# examples of inputs: see sign.input
# should produce no output: python sign.py < sign.input

# warning: currently 37 seconds/line on a fast machine

# fields on each input line: sk, pk, m, sm
# each field hex
# each field colon-terminated
# sk includes pk at end
# sm includes m at end

run = 0
while not (run == 1):
  line = sys.stdin.readline()
  if not line: break
  x = line.split(':')
  print("PARAMETERS:")
  print("sk = " + x[0][0:64])
  print("pk = " + x[1])
  print("m = " + x[2])
  print("sig = " + x[3][0:len(x[3])-len(x[2])])
  sk = binascii.unhexlify(x[0][0:64])
#  print("SHA512(skey) = " + binascii.hexlify(hashlib.sha512(sk).digest()).decode("utf-8"))
  print("-----------------------")
#  print("Testing the publickey() function... ", end="")
  pk = ed25519.publickey(sk)
#  print("ok")
  print("Checking the public key... ", end="")
  if binascii.hexlify(pk).decode("utf-8") == x[1]:
    print("same")
  else:
    print("different")
  print("-----------------------")
  
  
  m = binascii.unhexlify(x[2])
  s = ed25519.signature(m,sk,pk)
  ed25519.checkvalid(s,m,pk)
#  forgedsuccess = 0
#  try:
#    if len(m) == 0:
#      forgedm = "x"
#    else:
#      forgedmlen = len(m)
#      forgedm = ''.join([chr(ord(m[i])+(i==forgedmlen-1)) for i in range(forgedmlen)])
#    ed25519.checkvalid(s,forgedm,pk)
#    forgedsuccess = 1
#  except:
#    pass
#  assert not forgedsuccess
#  assert x[0] == binascii.hexlify(sk + pk).decode('utf-8')
#  assert x[1] == binascii.hexlify(pk).decode('utf-8')
#  assert x[3] == binascii.hexlify(s + m).decode('utf-8')
  print("-----------------------")
  run = 1
