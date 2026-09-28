Pairing clarified: controls must be compared with their matched-history direct runs, not default-history `direct`. All label/radius and `initial.tsv`/`modes.tsv` byte-identity assertions passed.

| Pair | meanCartesian | maxFractional |
|---|---:|---:|
| direct vs finufft14 | .4237308 | .0429480 |
| direct vs finufft15 | .4137843 | .0357714 |
| history1_direct vs history1_finufft | .1824786 | .0266998 |
| history5_direct vs history5_finufft | .2174896 | .0240911 |
| history10_direct vs history10_finufft | .2467948 | .0323664 |

Controls aggregate is endpointFAIL: all maxFractional values exceed `1e-7`.