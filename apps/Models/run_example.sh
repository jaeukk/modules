#!/usr/bin/env bash
# run_example.sh -- minimal usage examples for ./models.out
# Generates a small RSA ConfigPack and a Poisson SpherePacking into ./ .
#
# stdin token order (legacy protocol, see README.md / main.cpp header):
#   RSA:     dim N rho phi seed threads   nConfig prefix type
#   Poisson: dim N rho mode               nConfig prefix type
set -euo pipefail
cd "$(dirname "$0")"

make >/dev/null

echo "== list models =="
./models.out -h

echo
echo "== RSA: 2D, N=200, rho=1, phi=0.3, seed=7, 4 threads, 3 configs -> rsa_example.ConfigPack =="
./models.out models <<'EOF'
RSA
2
200
1
0.3
7
4
3
rsa_example
ConfigPack
EOF

echo
echo "== Poisson: 2D, N=100, rho=1, canonical, 2 configs -> poisson_example__%05d.txt =="
./models.out models <<'EOF'
Poisson
2
100
1
0
2
poisson_example
SpherePacking
EOF

echo
echo "Outputs:"
ls -1 rsa_example.ConfigPack rsa_example.ConfigPackIndex poisson_example__*.txt
