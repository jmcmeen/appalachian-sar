# Data Dictionary

Both files come from the appendices of Stout, Jessee & McMeen (2025): amphibian species counts and land areas for the southern Appalachian ecoregion.

## `nested.csv`

The nested model (Appendix 1): seven concentric areas of increasing size, each containing all the smaller ones.

| Column | Type | Description |
|---|---|---|
| `Locality` | string | Area ID, `Area 1` to `Area 7` (see the hierarchy below) |
| `Area` | float | Land area in km² |
| `AmphibianSpecies` | int | Amphibian species count (frogs + salamanders) |
| `FrogSpecies` | int | Frog species count |
| `SalamanderSpecies` | int | Salamander species count |

### Spatial hierarchy

| Area | Locality | Area (km²) |
|---|---|---:|
| 1 | Wetland pond at Steele Creek Park, Bristol, Tennessee | 0.0018 |
| 2 | Steele Creek Park, Bristol, Tennessee | 9.3 |
| 3 | Sullivan County, Tennessee | 1,114 |
| 4 | Northeastern Tennessee (five counties, as defined in Jessee et al. 2022) | 4,137 |
| 5 | Eastern Tennessee (as defined in Jessee et al. 2022) | 37,438 |
| 6 | Area 5 plus the western North Carolina high and low mountains and the Blue Ridge, Ridge and Valley, and Appalachian Plateau of western Virginia | 101,271 |
| 7 | Area 6 plus West Virginia and the Blue Ridge, Ridge and Valley, and Appalachian Plateau of Georgia: southern Appalachia as defined in the study | 178,758 |

Because each area contains the last, the rows are not independent observations.

## `island.csv`

The island model (Appendix 2): 26 sites within the study area, mostly public land, each with species occurrence data from exhaustive sampling or many collecting events over years or decades.

| Column | Type | Description |
|---|---|---|
| `Locality` | string | Site name (state abbreviations as in the paper) |
| `Area` | float | Land area in km² |
| `AmphibianSpecies` | int | Amphibian species count (frogs + salamanders) |
| `FrogSpecies` | int | Frog species count |
| `SalamanderSpecies` | int | Salamander species count |

The sites range from John's Bog (0.0081 km²) to Great Smokies NP (2,114 km²) and include parks, national recreation areas, counties and other well-surveyed localities in Tennessee, North Carolina, West Virginia, Virginia and Georgia. They are separate localities, but not all disjoint: UT Arboretum lies within Anderson Co. TN, Elk Knob SP largely within Watauga Co. NC, and part of Great Smokies NP within Sevier Co. TN.

The file starts with a UTF-8 byte-order mark; `pandas.read_csv` handles it.

## Provenance

Land areas are county or state land areas, or polygons digitized in Google Earth Pro (western Virginia, western North Carolina and northern Georgia), as described in the paper. Species counts and their sources are given in the paper's appendices:

> Stout, J.B., Jessee, L.D. & McMeen, J.N. (2025). Nested and island models for determining the species-area relationship of southern Appalachian amphibians. *Journal of North American Herpetology*, 2025(1), 1–7. https://doi.org/10.17161/jnah.v2025i1.21974
